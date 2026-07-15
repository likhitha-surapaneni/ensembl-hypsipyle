"""
Structural variant allele model.
"""

from typing import Any, List, Mapping, Optional


class StructuralVariantAllele:
    def __init__(self, allele_index: int, alt: str, variant: dict) -> None:
        self.ref_len = variant.length if hasattr(variant, "length") else 0
        self.variant = variant
        self.allele_index = allele_index
        self.alt = alt
        self.copy_number = self._get_copy_number_from_info(allele_index, variant)
        self.type = "StructuralVariantAllele"

    def get_name(self) -> str:
        allele_name = None
        if hasattr(self.variant, "info") and isinstance(self.variant.info, dict):
            raw_allele_name = self.variant.info.get("ALLELE_NAME")
            if isinstance(raw_allele_name, (list, tuple)):
                # allele_index of 1 means first ALT; 0 refers to ref
                if self.allele_index > 0 and self.allele_index - 1 < len(raw_allele_name):
                    allele_name = raw_allele_name[self.allele_index - 1]
                elif self.allele_index == 0 and self.variant.ref is not None:
                    # no ref-specific name in list, keep fallback
                    self.alt_len = self.ref_len
                    allele_name = f"{self.variant.chromosome}:{self.variant.position}:{self.ref_len}:{self.ref_len}"
            elif isinstance(raw_allele_name, str):
                allele_name = raw_allele_name

            if allele_name:
                self.name = allele_name
            else:
                if self.allele_index == 0 and self.variant.ref is not None:
                    self.alt_len = self.ref_len
                else:
                    allele_type = self.get_allele_type() 
                    self.alt_len = 0 if allele_type["value"] == "deletion" else self.get_length()
                self.name = f"{self.variant.chromosome}:{self.variant.position}:{self.ref_len}:{self.alt_len}"
        return self.name
    
    
    def get_allele_type(self) -> Mapping:
        return self.variant.get_allele_type(self.alt)

    def get_length(self) -> int:
        is_symbolic_alt = (
            isinstance(self.alt, str) and self.alt in ("DEL", "DUP", "INV", "INS", "CNV")
        )

        # Prefer SVLEN entry from INFO when available. Only non-symbolic ALTs carry sequence.
        self.alt_len = (
            len(self.alt)
            if self.allele_index > 0 and self.alt is not None and not is_symbolic_alt
            else 0
        )
        if hasattr(self.variant, "info") and isinstance(self.variant.info, dict):
            raw_svlen = self.variant.info.get("SVLEN") or self.variant.info.get("END") - self.variant.position + 1
            svlen_value = None
            if isinstance(raw_svlen, (list, tuple)):
                if self.allele_index > 0 and self.allele_index - 1 < len(raw_svlen):
                    svlen_value = raw_svlen[self.allele_index - 1]
            elif self.allele_index > 0:
                svlen_value = raw_svlen

            if svlen_value is not None:
                try:
                    self.alt_len = abs(int(svlen_value))
                except (TypeError, ValueError):
                    pass
        if "CN" in self.variant.info:
            if self.get_copy_number():
                return self.alt_len * self.get_copy_number()
        return self.alt_len

    def _get_copy_number_from_info(self, allele_index: int, variant: Any) -> Optional[int]:
        if self.get_allele_type()["value"] == "biological_region":
                return 1

        raw_copy_number = variant.info.get("CN")
        if raw_copy_number is None:
            return None

        if isinstance(raw_copy_number, (list, tuple)):
            copy_number_values = raw_copy_number
            if (
                len(copy_number_values) == 1
                and isinstance(copy_number_values[0], str)
                and "," in copy_number_values[0]
            ):
                copy_number_values = copy_number_values[0].rstrip(";").split(",")
        elif isinstance(raw_copy_number, str):
            copy_number_values = raw_copy_number.rstrip(";").split(",")
        else:
            copy_number_values = [raw_copy_number]

        copy_number_index = allele_index - 1
        if copy_number_index >= len(copy_number_values):
            return None

        copy_number = copy_number_values[copy_number_index]
        if isinstance(copy_number, str):
            copy_number = copy_number.strip().rstrip(";")
            if copy_number in ("", "."):
                return None

        try:
            return int(copy_number)
        except (TypeError, ValueError):
            return None

    def get_copy_number(self) -> Optional[int]:
        return self.copy_number

    def get_alternative_names(self) -> list:
        return self.variant.get_alternative_names(self.name)

    def get_slice(self) -> Mapping:
        return self.variant.get_slice(self.alt)

    def get_phenotype_assertions(self) -> list:
        return []

    def get_predicted_molecular_consequences(self) -> list:
        return []

    def get_prediction_results(self) -> list:
        return []

    def get_population_allele_frequencies(self) -> list:
        return []

    def get_web_display_data(self) -> Mapping:
        return {}

    def create_allele_prediction_results(
        self,
        current_prediction_results: Mapping,
        csq_record: List,
        prediction_index_map: dict,
    ) -> list:
        """Creates prediction results for the allele based on a CSQ record.

        Args:
            current_prediction_results (Mapping): Existing prediction results.
            csq_record (List): The CSQ record split into fields.
            prediction_index_map (dict): A mapping from annotation keys to their indices.

        Returns:
            list: A list of new prediction results.
        """
        prediction_results = []
        if "cadd_phred" in prediction_index_map.keys():
            if not self.prediction_result_already_exists(
                current_prediction_results, "CADD"
            ):
                cadd_prediction_result = (
                    {
                        "score": csq_record[prediction_index_map["cadd_phred"]],
                        "analysis_method": {
                            "tool": "CADD",
                            "qualifier": {"result_type": "CADD Phred score"},
                            "reference_data": [],
                        },
                    }
                    if csq_record[prediction_index_map["cadd_phred"]]
                    else None
                )
                if cadd_prediction_result:
                    prediction_results.append(cadd_prediction_result)

        return prediction_results

    def prediction_result_already_exists(
        self, current_prediction_results: Mapping, tool: str
    ) -> bool:
        """Checks if a prediction result for a specific tool already exists.

        Args:
            current_prediction_results (Mapping): The current prediction results.
            tool (str): The analysis tool name.

        Returns:
            bool: True if a result exists, False otherwise.
        """
        for prediction_result in current_prediction_results:
            if prediction_result["analysis_method"]["tool"] == tool:
                return True

        return False
