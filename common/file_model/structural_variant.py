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

from typing import Any, List

from common.file_model.base_variant import BaseVariant
from common.file_model.structural_variant_allele import StructuralVariantAllele


class StructuralVariant(BaseVariant):
    """StructuralVariant model – inherits shared behaviour from BaseVariant."""

    def __init__(
        self,
        record: Any,
        header: Any,
        genome_uuid: str,
        synonyms_by_allele_id: dict | None = None,
    ) -> None:
        """Initialise SV-specific attributes and delegate shared setup.

        The `length` attribute is derived from the VCF INFO `SVLEN` where
        available; otherwise zero.
        """
        super().__init__(record, header, genome_uuid)
        self.type = "StructuralVariant"
        self.synonyms_by_allele_id = synonyms_by_allele_id or {}
        # SVLEN may be a list or integer depending on the generator; coerce to int
        svlen = self.info.get("SVLEN") if isinstance(self.info, dict) else None
        try:
            if isinstance(svlen, (list, tuple)) and svlen:
                self.length = int(svlen[0])
            elif svlen is not None:
                self.length = int(svlen)
            else:
                self.length = 0
        except Exception:
            self.length = 0

    def get_primary_source(self) -> dict:
        return super().get_primary_source()

    def set_synonyms_by_allele_id(self, synonyms_by_allele_id: dict | None) -> None:
        self.synonyms_by_allele_id = synonyms_by_allele_id or {}

    def get_alternative_names(self, allele_id: str | None = None) -> list:
        if allele_id:
            return [
                self._build_synonym(synonym["synonym"], synonym.get("source"))
                for synonym in self.synonyms_by_allele_id.get(allele_id, [])
            ]

        synonyms = []
        seen = set()
        for allele_synonyms in self.synonyms_by_allele_id.values():
            for synonym in allele_synonyms:
                synonym_id = synonym["synonym"]
                source = synonym.get("source")
                key = (
                    synonym_id,
                    source,
                )
                if key not in seen:
                    seen.add(key)
                    synonyms.append(self._build_synonym(synonym_id, source))

        return synonyms

    @staticmethod
    def _build_synonym(synonym: str, source: str) -> dict:
        source = source or "dbVar"
        return {
            "accession_id": synonym,
            "name": synonym,
            "description": f"Structural variant synonym from {source}",
            "assignment_method": {
                "type": "DIRECT",
                "description": (
                    "A reference made by an external resource of annotation to an "
                    "Ensembl feature that Ensembl imports without modification"
                ),
            },
            "url": None,
            "source": {
                "id": source,
                "name": source,
                "description": "",
                "url": None,
                "release": None,
            },
        }

    @staticmethod
    def _build_allele_type_payload(allele_type: str, so_term: str) -> dict:
        return {
            "accession_id": allele_type,
            "value": allele_type,
            "url": f"http://sequenceontology.org/browser/current_release/term/{so_term}",
            "source": {
                "id": "",
                "name": "Sequence Ontology",
                "url": "www.sequenceontology.org",
                "description": "The Sequence Ontology...",
            },
        }

    def get_slice(self, allele=None) -> dict:
        target_allele = self.alts if allele is None else allele
        return super().get_slice(target_allele)

    def get_allele_type(self, allele: Any | None = None) -> dict:
        svtype = self.info.get("SVTYPE") if isinstance(self.info, dict) else None

        svtype_to_term = {
            "DEL": ("deletion", "SO:0000159"),
            "INS": ("insertion", "SO:0000667"),
            "DUP": ("duplication", "SO:1000035"),
            "INV": ("inversion", "SO:1000036"),
            "CNV": ("copy_number_variation", "SO:0001019"),
            "BND": ("translocation", "SO:0000199"),
        }

        if isinstance(svtype, str):
            normalized_svtype = svtype.upper()
            if normalized_svtype in svtype_to_term:
                allele_type, so_term = svtype_to_term[normalized_svtype]
                return self._build_allele_type_payload(allele_type, so_term)

        try:
            return super().get_allele_type(self.alts if allele is None else allele)
        except Exception:
            # Keep GraphQL non-null contract even when ALT/SVTYPE is malformed.
            return self._build_allele_type_payload("structural_variant", "SO:0001537")

    def get_length(self) -> int:
        return self.length

    def get_alleles(self) -> List:
        variant_allele_list = []
        alts = self.alts or []
        for index, alt in enumerate(alts):
            alt_value = alt.value if hasattr(alt, "value") else str(alt)
            variant_allele_list.append(StructuralVariantAllele(index + 1, alt_value, self))
        if self.ref:
            variant_allele_list.append(StructuralVariantAllele(0, self.ref, self))
        return variant_allele_list

    def get_prediction_results(self) -> list:
        return []

    def get_ensembl_website_display_data(self) -> dict:
        return {}

    def get_web_display_data(self) -> dict:
        return {}
