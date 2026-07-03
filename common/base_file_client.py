"""
Shared file and index helpers for variant clients.
"""

import io
import os
import subprocess
from common.file_model.variant import Variant
from common.file_model.structural_variant import StructuralVariant

import vcfpy


class BaseFileClient:
    """Common helpers for loading VCF records and locating related files."""

    INDEX_FILENAMES = ("variants.duckdb", "variants.db")
    record_class = None
    use_bcftools = False

    def split_variant_id(self, variant_id: str):
        """Splits a variant identifier in the form contig:position:identifier."""
        return variant_id.split(":")

    def search_in_file(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        source_name: str = None,
        record_class=None,
        use_bcftools: bool = None,
    ) -> None | Variant | StructuralVariant:
        """Search for a variant in a specific VCF file."""
        record_class = record_class or self.record_class
        use_bcftools = self.use_bcftools if use_bcftools is None else use_bcftools

        if not os.path.exists(datafile) or record_class is None:
            return None

        try:
            if use_bcftools:
                return self.get_record_with_bcftools(
                    datafile,
                    contig,
                    pos,
                    id,
                    genome_uuid,
                    record_class,
                    source_name,
                )

            return self.get_record_from_file(
                datafile, contig, pos, id, genome_uuid, record_class
            )
        except Exception as e:
            backend = "bcftools" if use_bcftools else "VCF"
            print(f"{backend} search failed for {source_name or datafile}: {str(e)}")
            return None

    def get_record_from_file(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        record_class,
        fallback_to_iteration: bool = False,
        error_prefix: str = None,
    ) -> None | Variant | StructuralVariant:
        """Read a matching record from a VCF file using vcfpy."""
        if not os.path.exists(datafile):
            return None

        reader = None
        try:
            reader = vcfpy.Reader.from_path(datafile)
            try:
                for rec in reader.fetch(contig, pos - 1, pos):
                    if self.record_matches(rec, contig, pos, id):
                        return record_class(rec, reader.header, genome_uuid)
            except (NotImplementedError, TypeError):
                if not fallback_to_iteration:
                    return None
                self.close_reader(reader)
                reader = vcfpy.Reader.from_path(datafile)
                for rec in reader:
                    if self.record_matches(rec, contig, pos, id):
                        return record_class(rec, reader.header, genome_uuid)
        except Exception as e:
            if error_prefix:
                print(f"{error_prefix}: {str(e)}")
            return None
        finally:
            self.close_reader(reader)

        return None

    def get_record_with_bcftools(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        record_class,
        source_name: str = None,
    ) -> None | Variant | StructuralVariant:
        """Read a matching record from a VCF file using bcftools view.
        bcftools gives us raw VCF text, but the rest of this codebase
        expects parsed VCF record objects
        """
        region = f"{contig}:{pos}-{pos}"
        try:
            result = subprocess.run(
                ["bcftools", "view", "-r", region, datafile],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
        except subprocess.TimeoutExpired:
            print(f"bcftools timeout for {source_name or datafile}")
            return None

        if result.returncode != 0:
            return None

        reader = vcfpy.Reader(io.StringIO(result.stdout))
        for rec in reader:
            if self.record_matches(rec, contig, pos, id):
                return record_class(rec, reader.header, genome_uuid)

        return None

    def get_full_record(
        self,
        datafile: str,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        record_class,
    ) -> None | Variant | StructuralVariant:
        """Read a full VCF record using vcfpy."""
        try:
            return self.get_record_from_file(
                datafile,
                contig,
                pos,
                id,
                genome_uuid,
                record_class,
                fallback_to_iteration=True,
            )
        except Exception as e:
            print(f"Error fetching full record from {datafile}: {str(e)}")
            return None

    def get_vcf_files(
        self, directory: str, prefix: str = "variation_", suffix: str = ".vcf.gz"
    ) -> list[str]:
        if not os.path.isdir(directory):
            return []

        return sorted(
            filename
            for filename in os.listdir(directory)
            if filename.startswith(prefix) and filename.endswith(suffix)
        )

    def search_all_files(
        self,
        directory: str,
        vcf_files: list,
        contig: str,
        pos: int,
        id: str,
        genome_uuid: str,
        record_class=None,
        use_bcftools: bool = None,
    ) -> None | Variant | StructuralVariant:
        """Search for a variant across all supplied VCF files."""
        try:
            for vcf_file in vcf_files:
                datafile = os.path.join(directory, vcf_file)
                variant = self.search_in_file(
                    datafile,
                    contig,
                    pos,
                    id,
                    genome_uuid,
                    record_class=record_class,
                    use_bcftools=use_bcftools,
                )
                if variant:
                    return variant

            return None
        except Exception as e:
            print(f"Error searching all variation files: {str(e)}")
            return None

    def get_index_db_path(self, directory: str) -> str | None:
        for filename in self.INDEX_FILENAMES:
            db_path = os.path.join(directory, filename)
            if os.path.exists(db_path):
                return db_path

        return None

    def fetch_duckdb_rows(self, db_path: str, query: str, parameters: list) -> list:
        if not db_path:
            return []

        try:
            import duckdb

            conn = duckdb.connect(db_path, read_only=True)
            try:
                return conn.execute(query, parameters).fetchall()
            finally:
                conn.close()
        except Exception:
            return []

    def record_matches(self, rec, contig: str, pos: int, id: str) -> bool:
        return rec.CHROM == contig and rec.POS == pos and rec.ID and rec.ID[0] == id

    def close_reader(self, reader) -> None:
        if reader is None:
            return

        close = getattr(reader, "close", None)
        if close:
            close()
