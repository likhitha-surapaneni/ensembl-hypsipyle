"""
Common base class for variant models.

This module defines `BaseVariant` which encapsulates behaviour shared
between short variants and structural variants.  Concrete subclasses
are responsible for any additional attributes (e.g. `alts` for
`Variant`) or specialised methods.
"""

from typing import Any, Mapping
import re
import os
import json
from common.file_model.utils import decode_population_name


class BaseVariant:
    variant_sources = {}  # cached per-genome source information

    def __init__(self, record: Any, header: Any, genome_uuid: str) -> None:
        # core attributes present in both VCF record types
        self.genome_uuid = genome_uuid
        self.name = self._normalize_variant_name(record.ID[0])
        self.record = record
        self.header = header
        self.chromosome = record.CHROM
        self.position = record.POS
        self.info = record.INFO
        # some structural VCFs may not have ALT/REF but most do; set if
        # available for convenience
        self.ref = getattr(record, "REF", None)
        self.alts = getattr(record, "ALT", None)
        # default type may be overridden by subclass
        self.type = self.__class__.__name__
        # obtain vep version from header; errors propagate so that tests
        # using dummy headers can short-circuit this by providing a
        # minimal implementation
        self.vep_version = re.search("v\d+", header.get_lines("VEP")[0].value).group()
        self.population_map = {}

    # ------------------------------------------------------------------
    # shared helpers
    # ------------------------------------------------------------------

    def _normalize_variant_name(self, name: str | None) -> str | None:
        if not name:
            return name

        if not hasattr(self, "_normalize_spdi_name"):
            return name
        return self._normalize_spdi_name(name)

    def parse_source_from_header(self) -> Mapping:
        """Parse and cache source metadata from the VCF header."""
        genome_uuid = self.genome_uuid
        if genome_uuid not in self.variant_sources:
            self.variant_sources[genome_uuid] = {}

        source_header_lines = self.header.get_lines("source")
        for source_header_line in source_header_lines:
            source, source_info_line = source_header_line.value.split('" ', 1)
            source = source.strip('"').replace(" ", "_")
            source_info = dict(re.findall('(.+?)="(.+?)"\s*', source_info_line))
            # overwrite allowed
            self.variant_sources[genome_uuid][source] = source_info
        return self.variant_sources[genome_uuid]

    def get_primary_source(self) -> Mapping:
        """Return the primary source for this variant."""
        try:
            if "SOURCE" in self.info:
                source = self.info["SOURCE"]
            else:
                # fall back to header default
                srcs = self.header.get_lines("source")
                if not srcs:
                    return None
                source = srcs[0].value.split()[0].strip('"')

            genome_uuid = self.genome_uuid
            if (
                genome_uuid not in self.variant_sources
                or source not in self.variant_sources[genome_uuid]
            ):
                self.parse_source_from_header()

            variant_sources = self.variant_sources[genome_uuid]
            if source in variant_sources:
                info = variant_sources[source]
                return {
                    "id": info.get("ID", source),
                    "name": source,
                    "description": info.get("description"),
                    "url": info.get("url"),
                    "release": info.get("version"),
                }
            return None
        except Exception:
            return None

    def set_allele_type(
        self, alt_one_bp: bool, ref_one_bp: bool, ref_alt_equal_bp: bool
    ):
        match [alt_one_bp, ref_one_bp, ref_alt_equal_bp]:
            case [True, True, True]:
                allele_type = "SNV"
                SO_term = "SO:0001483"
            case [True, False, False]:
                allele_type = "deletion"
                SO_term = "SO:0000159"
            case [False, True, False]:
                allele_type = "insertion"
                SO_term = "SO:0000667"
            case [False, False, False]:
                allele_type = "indel"
                SO_term = "SO:1000032"
            case [False, False, True]:
                allele_type = "substitution"
                SO_term = "SO:1000002"
        return allele_type, SO_term

    def get_allele_type(self, allele) -> Mapping:
        """Determine type of the supplied allele.

        The method is generic so that callers (including the GraphQL layer)
        need not know which subclass they are dealing with.
        """
        if isinstance(allele, str):
            if allele == self.ref:
                allele_type = "biological_region"
                SO_term = "SO:0001411"
            else:
                allele_type, SO_term = self.set_allele_type(
                    len(allele) == 1,
                    len(self.ref) == 1 if self.ref is not None else False,
                    len(allele) == len(self.ref) if self.ref is not None else False,
                )
        else:
            # allele is a list of objects with a `.value` attribute
            alt_length = max(len(a.value) for a in allele)
            allele_type, SO_term = self.set_allele_type(
                alt_length < 2, len(self.ref) < 2, alt_length == len(self.ref)
            )
        return {
            "accession_id": allele_type,
            "value": allele_type,
            "url": f"http://sequenceontology.org/browser/current_release/term/{SO_term}",
            "source": {
                "id": "",
                "name": "Sequence Ontology",
                "url": "www.sequenceontology.org",
                "description": "The Sequence Ontology...",
            },
        }

    def get_slice(self, allele) -> Mapping:
        """Return a location slice for the variant given an allele value."""
        start = self.position
        length = len(self.ref) if self.ref is not None else 0
        end = start + length - 1
        if allele != self.ref:
            allele_type = self.get_allele_type(allele)
            if allele_type["accession_id"] == "insertion":
                end = start
                length = 0
        return {
            "location": {"start": start, "end": end, "length": length},
            "region": {
                "name": self.chromosome,
                "code": "chromosome",
                "topology": "linear",
                "so_term": "SO:0001217",
            },
            "strand": {"code": "forward", "value": 1},
        }

    def get_most_severe_consequence(self) -> Mapping:
        """Return the most severe known consequence from the VEP CSQ data."""
        consequence_map = {}
        csq_records = self.info.get("CSQ", [])
        if csq_records:
            consequence_index = self.get_info_key_index("Consequence")
            if consequence_index is not None:
                directory = os.path.dirname(__file__)
                with open(
                    os.path.join(directory, "variation_consequence_rank.json")
                ) as rank_file:
                    consequence_rank = json.load(rank_file)

                for csq_record in csq_records:
                    csq_record_list = csq_record.split("|")
                    if consequence_index >= len(csq_record_list):
                        continue

                    for consequence in csq_record_list[consequence_index].split("&"):
                        rank = consequence_rank.get(consequence)
                        if rank is None:
                            continue
                        try:
                            consequence_map[int(rank)] = consequence
                        except (TypeError, ValueError):
                            continue

        return {
            "result": consequence_map[min(consequence_map.keys())]
            if consequence_map
            else None,
            "analysis_method": {
                "tool": "Ensembl VEP",
                "qualifier": {"result_type": "most severe consequence", "modes": []},
            },
        }

    def get_gerp_score(self) -> Mapping:
        """Return the first valid GERP conservation score from VEP CSQ data."""
        csq_records = self.info.get("CSQ", [])
        if not csq_records:
            return {}

        gerp_index = self.get_info_key_index("Conservation")
        if gerp_index is None:
            return {}

        for csq_record in csq_records:
            csq_record_list = csq_record.split("|")
            if gerp_index >= len(csq_record_list):
                continue

            raw_score = csq_record_list[gerp_index].strip()
            if not raw_score or raw_score == ".":
                continue

            try:
                score = float(raw_score)
            except (TypeError, ValueError):
                continue

            return {
                "score": score,
                "analysis_method": {
                    "tool": "GERP",
                    "qualifier": {"result_type": "GERP score", "modes": []},
                },
            }

        return {}

    def get_ancestral_allele(self) -> Mapping:
        """Retrieves the ancestral allele from the variant's CSQ record.

        Returns:
            Mapping: The ancestral allele prediction, or an empty mapping when
            it is not available.
        """
        csq = self.info.get("CSQ", [])
        if not csq:
            return {}

        aa_index = self.get_info_key_index("AA")
        if aa_index is None:
            return {}

        csq_record_list = csq[0].split("|")
        ancestral_allele = csq_record_list[aa_index]
        if not ancestral_allele or ancestral_allele == ".":
            return {}

        return {
            "result": ancestral_allele,
            "analysis_method": {
                "tool": "AncestralAllele",
                "qualifier": {"result_type": "Ancestral Allele", "modes": []},
                "version": "110",  # self.vep_version
            },
        }

    def get_info_key_index(self, key: str, info_id: str = "CSQ") -> int:
        info_field = self.header.get_info_field_info(info_id).description
        csq_list = info_field.split("Format: ")[1].split("|")
        for index, value in enumerate(csq_list):
            if value.lower() == key.lower():
                return index
        return None

    def get_csq_field_indices(self, keys, info_id: str = "CSQ") -> Mapping:
        prediction_index_map = {}
        for key in keys:
            index = self.get_info_key_index(key, info_id)
            if index is not None:
                prediction_index_map[key.lower()] = index
        return prediction_index_map

    def traverse_population_info(self) -> Mapping:
        """Return VEP population frequencies keyed by allele and population."""
        csq_records = self.info.get("CSQ", [])
        if not csq_records:
            return {}

        pop_mapping = self.parse_population_file()
        population_frequency_map = {}
        allele_index = self.get_info_key_index("Allele")
        if allele_index is None:
            return population_frequency_map

        for csq_record in csq_records:
            csq_record_list = csq_record.split("|")
            if allele_index >= len(csq_record_list):
                continue

            allele = csq_record_list[allele_index]
            if allele is None or allele in population_frequency_map:
                continue

            population_frequency_map[allele] = {}
            for pop in pop_mapping.values():
                for sub_pop in pop:
                    if sub_pop["name"] in population_frequency_map[allele]:
                        continue

                    allele_count = allele_number = allele_frequency = None
                    for freq_key, freq_val in sub_pop["fields"].items():
                        col_index = self.get_info_key_index(freq_val)
                        if (
                            col_index is not None
                            and col_index < len(csq_record_list)
                            and csq_record_list[col_index] is not None
                        ):
                            value = csq_record_list[col_index].split("&")[0] or None
                            if freq_key == "af":
                                allele_frequency = value
                            elif freq_key == "an":
                                allele_number = value
                            elif freq_key == "ac":
                                allele_count = value
                            else:
                                raise Exception("Frequency metric is not recognised")

                    if allele_frequency is None:
                        try:
                            allele_frequency = int(allele_count) / int(allele_number)
                        except Exception:
                            print(
                                "Cannot calculate AF using expression - "
                                f"{allele_count}/{allele_number}"
                            )

                    if allele_frequency is not None:
                        population_frequency_map[allele][sub_pop["name"]] = {
                            "population_name": decode_population_name(sub_pop["name"]),
                            "allele_frequency": float(allele_frequency),
                            "allele_count": allele_count,
                            "allele_number": allele_number,
                            "is_minor_allele": False,
                            "is_hpmaf": False,
                        }

        return population_frequency_map

    def parse_population_file(self) -> dict:
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory, "populations.json")) as pop_file:
            population_mappings = json.load(pop_file)
        return population_mappings.get(self.genome_uuid, {})

    def get_population_reference_allele(self) -> str:
        """Return the allele key used for reference-frequency calculations."""
        return self.ref

    def set_frequency_flags(self) -> Mapping:
        """Apply minor-allele and HPMAF flags to population frequencies."""
        pop_mapping = self.parse_population_file()
        pop_names = []
        for pop in pop_mapping.values():
            pop_names.extend([sub_pop["name"] for sub_pop in pop])

        hpmaf = []
        pop_frequency_map = self.traverse_population_info()
        if not pop_frequency_map:
            return pop_frequency_map

        pop_frequency_map_transpose = {
            pop_name: {
                pop_allele: pop_frequency_map[pop_allele][pop_name]
                for pop_allele in pop_frequency_map
                if pop_name in pop_frequency_map[pop_allele]
            }
            for pop_name in pop_names
        }

        for pop_name in pop_frequency_map_transpose:
            by_population = []
            for pop_allele, pop_allele_freq in pop_frequency_map_transpose[
                pop_name
            ].items():
                by_population.append(
                    [float(pop_allele_freq["allele_frequency"]), pop_allele, pop_name]
                )
            if not by_population:
                continue

            reference_allele = self.get_population_reference_allele()
            reference_frequency = 1 - float(sum(list(zip(*by_population))[0]))
            if 0 <= reference_frequency <= 1:
                population_frequency_ref = {
                    "population_name": pop_name,
                    "allele_frequency": reference_frequency,
                    "allele_count": None,
                    "allele_number": None,
                    "is_minor_allele": False,
                    "is_hpmaf": False,
                }
                if reference_allele not in pop_frequency_map:
                    pop_frequency_map[reference_allele] = {}
                pop_frequency_map[reference_allele][
                    pop_name
                ] = population_frequency_ref
                by_population.append(
                    [reference_frequency, reference_allele, pop_name]
                )

            by_population_sorted = sorted(by_population, key=lambda item: item[0])
            if len(by_population_sorted) >= 2:
                highest_frequency = by_population_sorted[-1][0]
                maf_frequency = None
                for pop in reversed(by_population_sorted[:-1]):
                    if pop[0] == highest_frequency:
                        continue
                    if pop[0] < highest_frequency and not maf_frequency:
                        maf_frequency, maf_allele, maf_population = pop
                        pop_frequency_map[maf_allele][maf_population][
                            "is_minor_allele"
                        ] = True
                        hpmaf.append([maf_frequency, maf_allele, maf_population])
                    elif (
                        maf_frequency
                        and pop[0] == maf_frequency
                        and maf_allele != reference_allele
                    ):
                        pop_frequency_map[maf_allele][maf_population][
                            "is_minor_allele"
                        ] = True
                        hpmaf.append([maf_frequency, maf_allele, maf_population])
                    elif maf_frequency and pop[0] < maf_frequency:
                        break

        if hpmaf:
            hpmaf_sorted = sorted(hpmaf, key=lambda item: item[0])
            hpmaf_frequency, hpmaf_allele, hpmaf_population = hpmaf_sorted[-1]
            pop_frequency_map[hpmaf_allele][hpmaf_population]["is_hpmaf"] = True
            for hpmaf_pop in reversed(hpmaf_sorted[:-1]):
                if hpmaf_pop[0] == hpmaf_frequency:
                    hpmaf_frequency, hpmaf_allele, hpmaf_population = hpmaf_pop
                    pop_frequency_map[hpmaf_allele][hpmaf_population][
                        "is_hpmaf"
                    ] = True
                elif hpmaf_pop[0] < hpmaf_frequency:
                    break

        return pop_frequency_map

    def get_web_display_data(self) -> Mapping:
        n_citations = self.info.get("NCITE", 0)
        return {"count_citations": n_citations}
