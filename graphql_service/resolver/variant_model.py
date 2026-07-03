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

from typing import Dict
from ariadne import QueryType, ObjectType
from graphql import GraphQLResolveInfo

from graphql_service.resolver.exceptions import VariantNotFoundError

# Define Query types for GraphQL
# Don't forget to import these into ariadne_app.py if you add a new type that use schema

QUERY_TYPE = QueryType()
VARIANT_TYPE = ObjectType("Variant")
VARIANT_ALLELE_TYPE = ObjectType("VariantAllele")


@QUERY_TYPE.field("variant")
async def resolve_variant(
    _,
    info: GraphQLResolveInfo,
    by_id: Dict[str, str] = None,
) -> Dict:
    "Load variants via variant id"

    {
        "type": "Variant",
        "variant_id": by_id["variant_id"],
        "genome_id": by_id["genome_id"],
    }
    file_client = info.context["file_client"]
    result = file_client.get_variant_record(
        by_id["genome_id"],
        by_id["variant_id"],
        by_id.get("source_name"),  # source_name is optional
    )
    if not result:
        raise VariantNotFoundError(by_id["variant_id"])
    return result


@VARIANT_TYPE.field("primary_source")
def primary_source(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load primary source for variant
    """
    return variant.get_primary_source()


@VARIANT_TYPE.field("allele_type")
def allele_type(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load allele_type for variant
    """
    return variant.get_allele_type(variant.alts)


@VARIANT_TYPE.field("alternative_names")
def alternative_names(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load alternative names for variant
    """
    return variant.get_alternative_names()


@VARIANT_TYPE.field("slice")
def slice(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load slice for variant
    """
    return variant.get_slice(variant.alts)


@VARIANT_TYPE.field("prediction_results")
def prediction_results(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load prediction result for variant
    """
    prediction_results = []
    prediction_results.append(variant.get_most_severe_consequence())
    if variant.get_gerp_score():
        prediction_results.append(variant.get_gerp_score())
    if variant.get_ancestral_allele():
        prediction_results.append(variant.get_ancestral_allele())
    return prediction_results


@VARIANT_TYPE.field("ensembl_website_display_data")
def ensembl_website_display_data(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load ensembl website display data for variant
    """
    return variant.get_web_display_data()


@VARIANT_TYPE.field("alleles")
def resolve_alleles_from_variant(variant: Dict, info: GraphQLResolveInfo) -> Dict:
    """
    Load alleles for variant
    """
    return variant.get_alleles()


@VARIANT_ALLELE_TYPE.field("name")
def resolve_name_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load name for variant allele
    """
    return variant_allele.name


@VARIANT_ALLELE_TYPE.field("alternative_names")
def resolve_alternative_names_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load alternative names for variant allele
    """
    return variant_allele.get_alternative_names()


@VARIANT_ALLELE_TYPE.field("slice")
def resolve_slice_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load slice for variant allele
    """
    return variant_allele.get_slice()


@VARIANT_ALLELE_TYPE.field("allele_type")
def resolve_allele_type_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load allele type for variant allele
    """
    return variant_allele.get_allele_type()


@VARIANT_ALLELE_TYPE.field("phenotype_assertions")
def resolve_phenotype_assertions_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load phenotype assertions for variant allele
    """
    return variant_allele.get_phenotype_assertions()


@VARIANT_ALLELE_TYPE.field("predicted_molecular_consequences")
def resolve_predicted_molecular_consequences_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load predicted molecular consequences for variant allele
    """
    return variant_allele.get_predicted_molecular_consequences()


@VARIANT_ALLELE_TYPE.field("prediction_results")
def resolve_prediction_results_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load prediction results for variant allele
    """
    return variant_allele.get_prediction_results()


@VARIANT_ALLELE_TYPE.field("population_frequencies")
def resolve_population_frequencies_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load population frequencies for variant allele
    """
    return variant_allele.get_population_allele_frequencies()


@VARIANT_ALLELE_TYPE.field("ensembl_website_display_data")
def resolve_ensmebl_website_display_data_from_variant_allele(
    variant_allele: Dict, info: GraphQLResolveInfo
) -> Dict:
    """
    Load ensembl website display data for variant allele
    """
    return variant_allele.get_web_display_data()


@QUERY_TYPE.field("version")
def resolve_api(
    _: None, info: GraphQLResolveInfo
) -> Dict:  # the second argument must be named `info` to avoid a NameError
    return {"api": {"major": "0", "minor": "1", "patch": "0-beta"}}

