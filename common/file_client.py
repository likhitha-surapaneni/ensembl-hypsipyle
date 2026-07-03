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

import os
from common.base_file_client import BaseFileClient
from common.file_model.variant import Variant
from common.file_model.structural_variant import StructuralVariant
from common.structural_variant_index_client import StructuralVariantIndexClient


class FileClient(BaseFileClient):
    """Client to load file into memory.

    This class provides methods to retrieve and process variant data from VCF files.
    """

    record_class = Variant

    def __init__(self, config) -> None:
        """Initialises a FileClient instance.

        Args:
            config (dict): A configuration dictionary containing at least the "data_root" key.
        """
        self.data_root = config.get("data_root")

    def get_variant_record(
        self, genome_uuid: str, variant_id: str, source_name: str = None
    ) -> Variant | None:
        """Retrieves a variant entry using the specified variant identifier.

        This method loads the VCF file from the configured data root and genome UUID,
        then attempts to fetch the variant record that matches the given variant_id.
        The variant_id should be formatted as 'contig:position:identifier'.

        Args:
            genome_uuid (str): The UUID for the genome.
            variant_id (str): The variant identifier in the format 'contig:position:identifier'.
            source_name (str, optional): The name of the source for the variant.
                                         If None, searches the default file first and
                                         then any available source files.

        Returns:
            Variant or None: A Variant instance if a matching record is found; otherwise, None.
        """
        try:
            [contig, pos, id] = self.split_variant_id(variant_id)
            pos = int(pos)
        except (TypeError, ValueError) as e:
            # TODO: This needs to go to thoas logger
            print(
                f"Invalid variant_id format '{variant_id}': {str(e)}. Expected format: contig:position:identifier"
            )
            return None

        genome_dir = os.path.join(self.data_root, genome_uuid)
        if source_name:
            datafile = os.path.join(
                genome_dir,
                f"variation/variation_{source_name}.vcf.gz",
            )
            variant = self.search_in_file(datafile, contig, pos, id, genome_uuid)
            if variant:
                return variant
            print(f"Variant not found in source '{source_name}'")
            return None

        default_datafile = os.path.join(genome_dir, "variation.vcf.gz")
        if os.path.exists(default_datafile):
            variant = self.search_in_file(default_datafile, contig, pos, id, genome_uuid)
            if variant:
                return variant

        variation_dir = os.path.join(genome_dir, "variation")
        if not os.path.isdir(variation_dir):
            print(f"Variation directory not found: {variation_dir}")
            return None

        vcf_files = self.get_vcf_files(variation_dir)
        if not vcf_files:
            print(f"No variation VCF files found in {variation_dir}")
            return None

        return self.search_all_files(
            variation_dir, vcf_files, contig, pos, id, genome_uuid
        )

    def get_structural_variant_record(
        self, genome_uuid: str, variant_id: str, source_name: str = None
    ) -> StructuralVariant | None:
        """Retrieves a structural variant entry using the specified variant identifier.

        This method loads the VCF file from the configured data root, genome UUID and source name,
        then attempts to fetch the structural variant record that matches the given structural_variant_id.
        The structural_variant_id should be formatted as 'contig:position:identifier'.

        If source_name is not provided, searches across all available structural variation source files at once.

        Args:
            genome_uuid (str): The UUID for the genome.
            variant_id (str): The variant identifier in the format 'contig:position:identifier'.
            source_name (str, optional): The name of the source for the structural variant.
                                        If None, searches all available sources.

        Returns:
            StructuralVariant or None: A StructuralVariant instance if a matching record is found; otherwise, None.
        """
        try:
            [contig, pos, id] = self.split_variant_id(variant_id)
            pos = int(pos)
        except (TypeError, ValueError) as e:
            print(
                f"Invalid variant_id format '{variant_id}': {str(e)}. Expected format: contig:position:identifier"
            )
            return None

        # If source_name is provided, search that specific file
        if source_name:
            datafile = os.path.join(
                self.data_root,
                genome_uuid,
                f"structural-variation/variation_{source_name}.vcf.gz",
            )
            variant = self.search_in_file(
                datafile,
                contig,
                pos,
                id,
                genome_uuid,
                source_name,
                record_class=StructuralVariant,
                use_bcftools=True,
            )
            if variant:
                index_client = StructuralVariantIndexClient(os.path.dirname(datafile))
                index_client.load_structural_variant_synonyms(variant)
                return variant
            print(f"Structural variant not found in source '{source_name}'")
            return None
        else:
            # Search across all available structural variation sources at once using bcftools
            sv_dir = os.path.join(self.data_root, genome_uuid, "structural-variation")
            if not os.path.isdir(sv_dir):
                print(f"Structural variation directory not found: {sv_dir}")
                return None

            # Get all valid VCF files
            vcf_files = self.get_vcf_files(sv_dir)

            if not vcf_files:
                print(f"No structural variation VCF files found in {sv_dir}")
                return None

            index_client = StructuralVariantIndexClient(sv_dir)

            # If there is just one source file, prefer the single-file search implementation
            # instead of invoking the more general multi-file scan.  This keeps behaviour
            # consistent with consumers that pass a ``source_name`` and avoids building
            # a bcftools command with a single path for no reason.
            if len(vcf_files) == 1:
                single = vcf_files[0]
                datafile = os.path.join(sv_dir, single)
                print(
                    f"Only one structural variant file ({single}) found; using single-file lookup"
                )
                variant = self.search_in_file(
                    datafile,
                    contig,
                    pos,
                    id,
                    genome_uuid,
                    record_class=StructuralVariant,
                    use_bcftools=True,
                )
                if variant:
                    index_client.load_structural_variant_synonyms(variant)
                return variant

            variant = None
            if index_client.exists():
                variant = self._get_structural_variant_from_index(
                    index_client, sv_dir, contig, pos, id, genome_uuid
                )
                if not variant:
                    print("DuckDB index search returned no matching record")
            else:
                print(f"DuckDB index not found in {sv_dir}")
                variant = self.search_all_files(
                    sv_dir,
                    vcf_files,
                    contig,
                    pos,
                    id,
                    genome_uuid,
                    record_class=StructuralVariant,
                    use_bcftools=True,
                )

            if variant:
                index_client.load_structural_variant_synonyms(variant)
                return variant

        return None

    def _get_structural_variant_from_index(
        self,
        index_client: StructuralVariantIndexClient,
        sv_dir: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
    ) -> StructuralVariant | None:
        results = index_client.get_variant_locations(id)
        if not results:
            print(f"Variant {id} not found in DuckDB index")
            return None

        for result in results:
            if result["chr"] == contig and result["start"] == pos:
                datafile = os.path.join(sv_dir, result["source_file"])
                return self.search_in_file(
                    datafile,
                    contig,
                    pos,
                    id,
                    genome_uuid,
                    record_class=StructuralVariant,
                    use_bcftools=True,
                )

        print(f"Variant {id} found in index but no coordinates matched")
        return None
