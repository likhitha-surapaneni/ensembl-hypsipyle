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

from typing import Dict, Optional
from ariadne import ObjectType
from graphql import GraphQLResolveInfo

from graphql_service.resolver.exceptions import VariantNotFoundError
from graphql_service.resolver.variant_model import QUERY_TYPE
# Define Query types for GraphQL
# Don't forget to import these into ariadne_app.py if you add a new type


STRUCTURAL_VARIANT_TYPE = ObjectType("StructuralVariant")
STRUCTURAL_VARIANT_ALLELE_TYPE = ObjectType("StructuralVariantAllele")


@QUERY_TYPE.field("structural_variant")
async def resolve_structural_variant(
    _,
    info: GraphQLResolveInfo,
    by_id: Dict[str, str] = None,
) -> Dict:
    "Load structural variants via variant id"

    file_client = info.context["file_client"]
    result = file_client.get_structural_variant_record(
        by_id["genome_id"],
        by_id["variant_id"],
        by_id.get("source_name"),  # source_name is optional
    )
    if not result:
        raise VariantNotFoundError(by_id["variant_id"])
    return result


@STRUCTURAL_VARIANT_TYPE.field("primary_source")
def primary_source(structural_variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load primary source for variant
    """
    return structural_variant.get_primary_source()


@STRUCTURAL_VARIANT_TYPE.field("allele_type")
def allele_type(structural_variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load allele_type for variant
    """
    return structural_variant.get_allele_type()


@STRUCTURAL_VARIANT_TYPE.field("alternative_names")
def alternative_names(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load alternative names for variant
    """
    return variant.get_alternative_names()


@STRUCTURAL_VARIANT_TYPE.field("slice")
def slice(structural_variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load slice for variant
    """
    return structural_variant.get_slice()


@STRUCTURAL_VARIANT_TYPE.field("prediction_results")
def prediction_results(structural_variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load prediction result for variant
    """
    prediction_results = []
    prediction_results.append(structural_variant.get_most_severe_consequence())
    if structural_variant.get_gerp_score():
        prediction_results.append(structural_variant.get_gerp_score())
    if structural_variant.get_ancestral_allele():
        prediction_results.append(structural_variant.get_ancestral_allele())
    return prediction_results


# @STRUCTURAL_VARIANT_TYPE.field("ensembl_website_display_data")
# def ensembl_website_display_data(structural_variant: Dict, info: GraphQLResolveInfo) -> Dict:
#     """
#     Load ensembl website display data for variant
#     """
#     return structural_variant.get_web_display_data()


@STRUCTURAL_VARIANT_TYPE.field("alleles")
def resolve_alleles_from_structural_variant(structural_variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load alleles for variant
    """
    return structural_variant.get_alleles()    


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("name")
def resolve_name_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load name for variant allele
    """
    return structural_variant_allele.get_name()


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("alternative_names")
def resolve_alternative_names_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load alternative names for variant allele
    """
    return structural_variant_allele.get_alternative_names()


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("length")
def resolve_length_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> int:
    """
    Load length for variant allele
    """
    return int(structural_variant_allele.get_length())


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("copy_number")
def resolve_copy_number_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Optional[int]:
    """
    Load copy number for variant allele
    """
    copy_number = structural_variant_allele.get_copy_number()
    return int(copy_number) if copy_number is not None else None


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("slice")
def resolve_slice_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load slice for variant allele
    """
    return structural_variant_allele.get_slice()


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("allele_type")
def resolve_allele_type_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load allele type for variant allele
    """
    return structural_variant_allele.get_allele_type()


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("phenotype_assertions")
def resolve_phenotype_assertions_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load phenotype assertions for variant allele
    """
    return structural_variant_allele.get_phenotype_assertions()


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("prediction_results")
def resolve_prediction_results_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load prediction results for variant allele
    """
    return structural_variant_allele.get_prediction_results()


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("predicted_molecular_consequences")
def resolve_predicted_molecular_consequences_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load predicted molecular consequences for variant allele
    """
    return structural_variant_allele.get_predicted_molecular_consequences()


@STRUCTURAL_VARIANT_ALLELE_TYPE.field("population_frequencies")
def resolve_population_frequencies_from_structural_variant_allele(
    structural_variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load population frequencies for variant allele
    """
    return structural_variant_allele.get_population_allele_frequencies()


# @STRUCTURAL_VARIANT_ALLELE_TYPE.field("ensembl_website_display_data")
# def resolve_ensembl_website_display_data_from_variant_allele(
#     variant_allele: Dict, info: GraphQLResolveInfo
# ) -> Dict:
#     """
#     Load ensembl website display data for variant allele
#     """
#     return variant_allele.get_web_display_data()
