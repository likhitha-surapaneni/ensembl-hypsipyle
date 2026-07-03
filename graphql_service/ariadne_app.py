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

from typing import Dict, Callable

import ariadne
from graphql import GraphQLSchema
from starlette.requests import Request

from graphql_service.resolver.variant_model import (
    QUERY_TYPE,
    VARIANT_TYPE,
    VARIANT_ALLELE_TYPE,
)
from graphql_service.resolver.structural_variant_model import (
    QUERY_TYPE as STRUCTURAL_QUERY_TYPE,
    STRUCTURAL_VARIANT_TYPE,
    STRUCTURAL_VARIANT_ALLELE_TYPE,
)
from graphql_service.resolver.population_model import (
    QUERY_TYPE as POPULATION_QUERY_TYPE,
    POPULATION_TYPE,
)


def prepare_executable_schema() -> GraphQLSchema:
    """Combines schema definitions with the corresponding resolvers to produce an executable schema.

    Loads the GraphQL schema from the "common/schemas" directory and integrates it with the
    query and variant resolvers.

    Returns:
        GraphQLSchema: The executable GraphQL schema.
    """
    schema = ariadne.load_schema_from_path("common/schemas")
    return ariadne.make_executable_schema(
        schema,
        QUERY_TYPE,
        STRUCTURAL_QUERY_TYPE,
        POPULATION_QUERY_TYPE,
        VARIANT_TYPE,
        VARIANT_ALLELE_TYPE,
        POPULATION_TYPE,
        STRUCTURAL_VARIANT_TYPE,
        STRUCTURAL_VARIANT_ALLELE_TYPE,
    )


def prepare_context_provider(context: Dict) -> Callable[[Request], Dict]:
    """Creates a context provider function for GraphQL executions.

    Returns a closure that injects a fresh context for each request. The context will include
    the incoming request and the shared file client, ensuring that request-specific data does
    not leak between executions.

    Args:
        context (Dict): A dictionary containing shared objects (e.g. file_client) for the application.

    Returns:
        Callable[[Request], Dict]: A function that takes a Request and returns a context dictionary.
    """
    file_client = context["file_client"]

    def context_provider(request: Request) -> Dict:
        """Provides a fresh context for each GraphQL execution.

        Args:
            request (Request): The incoming HTTP request.

        Returns:
            Dict: A dictionary containing the request and the shared file_client.
        """
        return {
            "request": request,
            "file_client": file_client,
        }

    return context_provider
