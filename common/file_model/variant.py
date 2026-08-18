"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from typing import Any, Mapping, List, Union
import re
from common.file_model.base_variant import BaseVariant
from common.file_model.variant_allele import VariantAllele
from common.file_model.utils import minimise_allele


def reduce_allele_length(allele_list: List):
    """Returns the maximum length of allele values in the list.

    Args:
        allele_list (List): A list of allele objects with a 'value' attribute.

    Returns:
        int: The maximum length of allele.value, or -1 if no allele has a value longer than -1.
    """
    allele_length = -1
    for allele in allele_list:
        if len(allele.value) > allele_length:
            allele_length = len(allele.value)
    return allele_length


class Variant(BaseVariant):
    def __init__(self, record: Any, header: Any, genome_uuid: str) -> None:
        """Initialise a Variant and delegate shared setup to BaseVariant.

        Variant-specific attributes (like `alts`) are set here.
        """
        super().__init__(record, header, genome_uuid)
        # Variant-specific
        self.alts = getattr(record, "ALT", [])
        self.ref = getattr(record, "REF", None)
        self.type = "Variant"

    def get_alternative_names(self) -> List:
        """Returns alternative names for the variant.

        Returns:
            List: A list of alternative names.
        """
        return []

    def get_primary_source(self) -> Mapping:
        """Fetches the primary source information for the variant.

        It attempts to retrieve the source from the variant INFO columns and falls back
        to header information if necessary.

        Returns:
            Mapping: A mapping containing source details, or None if an error occurs.
        """

        try:
            if "SOURCE" in self.info:
                source = self.info["SOURCE"]
            else:
                source = self.header.get_lines("source")[0].value

            # Get source information from data file header header
            genome_uuid = self.genome_uuid
            if (
                genome_uuid not in self.variant_sources
                or source not in self.variant_sources[genome_uuid]
            ):
                self.parse_source_from_header()
            variant_sources = self.variant_sources[genome_uuid]

            if source in variant_sources:
                source_info = variant_sources[source]

                source_id = source
                source_name = source.replace("_", " ")
                source_description = (
                    source_info["description"] if "description" in source_info else ""
                )
                source_url = source_info["url"] if "url" in source_info else ""
                source_release = (
                    source_info["version"] if "version" in source_info else ""
                )

                if "accession_url" in source_info:
                    source_url_id = source_info["accession_url"]
                    if re.search("^Ensembl", source):
                        variant_id = f"{self.chromosome}:{self.position}:{self.name}"
                    else:
                        variant_id = self.name
                else:
                    source_url_id = source_url
                    variant_id = ""

            # If source information not found in data file try using default value for main accessioning sources
            elif re.search("^dbSNP", source):
                source_id = "dbSNP"
                source_name = "dbSNP"
                source_description = "NCBI db of human variants"
                source_url = "https://www.ncbi.nlm.nih.gov/snp/"
                source_url_id = source_url
                source_release = 156
                variant_id = self.name

            elif re.search("^EVA", source):
                source_id = "EVA"
                source_name = "EVA"
                source_description = "European Variation Archive"
                source_url = "https://www.ebi.ac.uk/eva"
                source_url_id = "https://www.ebi.ac.uk/eva/?variant&accessionID="
                source_release = "release_6"
                variant_id = self.name

            elif re.search("^Ensembl", source):
                source_id = "Ensembl"
                source_name = "Ensembl"
                source_description = "Ensembl"
                source_url = "https://beta.ensembl.org"
                source_url_id = "https://beta.ensembl.org/"
                source_release = "110"  # to be fetched from the file
                variant_id = f"{self.chromosome}:{self.position}:{self.name}"

        except Exception:
            return None

        return {
            "accession_id": self.name,
            "name": self.name,
            "description": f"{source_description}",
            "assignment_method": {
                "type": "DIRECT",
                "description": "A reference made by an external resource of annotation to an Ensembl feature that Ensembl imports without modification",
            },
            "url": f"{source_url_id}{variant_id}",
            "source": {
                "id": f"{source_id}",
                "name": f"{source_name}",
                "description": f"{source_description}",
                "url": f"{source_url}",
                "release": f"{source_release}",
            },
        }

    def get_alleles(self) -> List:
        """Generates a list of VariantAllele instances for the variant.

        Returns:
            List: A list of VariantAllele objects, including both alternate and reference alleles.
        """
        variant_allele_list = []

        for index, alt in enumerate(self.alts):
            if index + 1 <= len(self.alts):
                variant_allele = VariantAllele(index + 1, alt.value, self)
                variant_allele_list.append(variant_allele)
        reference_allele = VariantAllele(0, self.ref, self)
        variant_allele_list.append(reference_allele)
        return variant_allele_list

    def get_population_reference_allele(self) -> str:
        """Return the minimised reference key used by short-variant CSQ data."""
        return minimise_allele(self.ref, self.ref)

    def get_statistics_info(self) -> Mapping:
        """Collects statistical information for each allele of the variant.

        Returns:
            Mapping: A mapping where each allele is associated with various counts and frequency data.
        """
        alleles = [i.value for i in self.alts]
        statistics_info = {}
        RAF_exists = True
        if not self.parse_population_file():
            RAF_exists = False
            print(f"No representative allele frequency for - {self.genome_uuid}")

        for index, allele in enumerate(alleles):
            statistics_info[allele] = {
                "count_transcript_consequences": self.info["NTCSQ"][index]
                if "NTCSQ" in self.info
                else 0,
                "count_overlapped_genes": self.info["NGENE"][index]
                if "NGENE" in self.info
                else 0,
                "count_regulatory_consequences": self.info["NRCSQ"][index]
                if "NRCSQ" in self.info
                else 0,
                "count_variant_phenotypes": self.info["NVPHN"][index]
                if "NVPHN" in self.info
                else 0,
                "count_gene_phenotypes": self.info["NGPHN"][index]
                if "NGPHN" in self.info
                else 0,
                "representative_population_allele_frequency": self.info["RAF"][index]
                if "RAF" in self.info and RAF_exists
                else None,
            }

        statistics_info[self.ref] = {
            "count_transcript_consequences": 0,
            "count_overlapped_genes": 0,
            "count_regulatory_consequences": 0,
            "count_variant_phenotypes": 0,
            "count_gene_phenotypes": 0,
            "representative_population_allele_frequency": 1
            - float(sum(filter(None, self.info["RAF"])))
            if "RAF" in self.info and RAF_exists
            else None,
        }

        return statistics_info
