#!/usr/bin/env bash
# This finds all the vcf files in DATA_DIR of the genome_uuid and creates a row for each individual allele id
# with columns CHR POS ID ALLELE_INDEX ALLELE ALLELE_ID SOURCE_FILE
# This is then piped to duckdb to produce the index tables
set -euo pipefail

VCF_DIR="${1:-.}"
DB_FILE="${2:-variants.duckdb}"
SYNONYMS_FILE="${3:-$VCF_DIR/nstd-syn.txt}"

if [ ! -d "$VCF_DIR" ]; then
  echo "Error: VCF directory not found: $VCF_DIR" >&2
  exit 1
fi

VCF_FILE_COUNT=$(find "$VCF_DIR" -maxdepth 1 -type f -name "*.vcf.gz" | wc -l | tr -d '[:space:]')
if [ "$VCF_FILE_COUNT" -eq 0 ]; then
  echo "Error: no .vcf.gz files found in: $VCF_DIR" >&2
  echo "Usage: bash scripts/source/create_index.sh <vcf_dir> [db_file] [synonyms_file]" >&2
  exit 1
fi

if [ -f "$SYNONYMS_FILE" ]; then
  SYNONYMS_SQL_FILE="${SYNONYMS_FILE//\'/\'\'}"
  SYNONYMS_SQL="
CREATE OR REPLACE TABLE allele_synonyms AS
SELECT
    variation_name AS allele_id,
    alias,
    name AS source
FROM read_csv('$SYNONYMS_SQL_FILE',
    delim='\t',
    header=true,
    columns={
        'variation_name': 'VARCHAR',
        'alias': 'VARCHAR',
        'name': 'VARCHAR'
    }
)
WHERE variation_name IS NOT NULL
  AND alias IS NOT NULL
  AND alias != '';
"
else
  SYNONYMS_SQL="
CREATE OR REPLACE TABLE allele_synonyms (
    allele_id VARCHAR,
    alias VARCHAR,
    source VARCHAR
);
"
fi

find "$VCF_DIR" -maxdepth 1 -type f -name "*.vcf.gz" | parallel -j 8 \
  "bcftools query -u -f '%CHROM\t%POS\t%ID\t%ALT\t%INFO/ALLELE_NAME\n' {} | awk -v f='{}' 'BEGIN{OFS=\"\t\"; n=split(f,p,\"/\"); base=p[n]} {delete alts; delete names; alt_count=split(\$4,alts,\",\"); name_count=(\$5 == \".\" || \$5 == \"\" ? 0 : split(\$5,names,\",\")); max_count=alt_count; if (name_count > max_count) max_count=name_count; for (i=1; i<=max_count; i++) {allele=(i in alts ? alts[i] : \"\"); allele_id=(i in names ? names[i] : \"\"); print \$1, \$2, \$3, i, allele, allele_id, base}}'" \
| duckdb "$DB_FILE" -c "
CREATE OR REPLACE TABLE variant_alleles AS
SELECT *
FROM read_csv('/dev/stdin',
    delim='\t',
    columns={
        'chr': 'VARCHAR',
        'start': 'INTEGER',
        'id': 'VARCHAR',
        'allele_index': 'INTEGER',
        'allele': 'VARCHAR',
        'allele_id': 'VARCHAR',
        'source_file': 'VARCHAR'
    }
)
ORDER BY chr, start, id, allele_index;

CREATE OR REPLACE TABLE variants AS
SELECT
    chr,
    start,
    id,
    string_agg(allele, ',' ORDER BY allele_index) AS allele,
    string_agg(allele_id, ',' ORDER BY allele_index) FILTER (WHERE allele_id IS NOT NULL AND allele_id != '') AS allele_ids,
    source_file
FROM variant_alleles
GROUP BY chr, start, id, source_file
ORDER BY chr, start, id;

$SYNONYMS_SQL

CREATE OR REPLACE TABLE variant_aliases AS
SELECT
    va.chr,
    va.start,
    va.id,
    va.allele_index,
    va.allele,
    va.allele_id,
    s.alias,
    s.source,
    va.source_file
FROM variant_alleles va
JOIN allele_synonyms s ON s.allele_id = va.allele_id
ORDER BY va.chr, va.start, va.id, va.allele_index, s.alias;

CREATE INDEX IF NOT EXISTS idx_variants_id ON variants(id);
CREATE INDEX IF NOT EXISTS idx_variants_chr_start ON variants(chr, start);
CREATE INDEX IF NOT EXISTS idx_variant_alleles_allele_id ON variant_alleles(allele_id);
CREATE INDEX IF NOT EXISTS idx_allele_synonyms_allele_id ON allele_synonyms(allele_id);
CREATE INDEX IF NOT EXISTS idx_variant_aliases_id ON variant_aliases(id);
CREATE INDEX IF NOT EXISTS idx_variant_aliases_alias ON variant_aliases(alias);
"
