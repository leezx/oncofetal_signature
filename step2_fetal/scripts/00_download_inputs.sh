#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
step_dir="$(cd "$script_dir/.." && pwd)"
source "$step_dir/config/step2.env"

hgca_raw="$STEP2_DATA_ROOT/scRNAseq/HGCA_Elmentaite2021/raw"
gao_raw="$STEP2_DATA_ROOT/scRNAseq/GSE103239_Gao2018/raw"
mouse_raw="$STEP2_DATA_ROOT/bulkRNAseq/GSE230581/raw"
orth_raw="$STEP2_DATA_ROOT/1.Databases/Ensembl_orthologues/release_116/raw"
mouse_reuse="$STEP2_DATA_ROOT/2.PROJECTS/1.TWEAKR-oncoFetal/data/external/Perturbation/GSE230581_CA-YAP-S127A_organoid/GSE230581_in_vivo_counts.tsv.gz"

mkdir -p "$hgca_raw" "$gao_raw" "$mouse_raw" "$orth_raw"

download_if_missing() {
  local url="$1"
  local output="$2"
  if [[ ! -s "$output" ]]; then
    curl -L --fail --retry 3 -o "$output" "$url"
  fi
}

download_if_missing \
  "https://cellgeni.cog.sanger.ac.uk/gutcellatlas/epi_raw_counts02_v2.h5ad" \
  "$hgca_raw/epi_raw_counts02_v2.h5ad"
download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95630/suppl/GSE95630_Digestion_TPM_new.txt.gz" \
  "$gao_raw/GSE95630_Digestion_TPM_new.txt.gz"
download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95630/suppl/GSE95630_README_barcodes.xlsx" \
  "$gao_raw/GSE95630_README_barcodes.xlsx"
download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103154/suppl/GSE103154_All_Merge_umi_tpm_gene.txt.gz" \
  "$gao_raw/GSE103154_All_Merge_umi_tpm_gene.txt.gz"
download_if_missing \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103154/suppl/GSE103154_README.xlsx" \
  "$gao_raw/GSE103154_README.xlsx"
download_if_missing \
  "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41556-018-0105-4/MediaObjects/41556_2018_105_MOESM4_ESM.xlsx" \
  "$gao_raw/41556_2018_105_MOESM4_ESM.xlsx"

if [[ ! -s "$MOUSE_COUNTS" ]]; then
  cp "$mouse_reuse" "$MOUSE_COUNTS"
fi

if [[ ! -s "$ORTHOLOGUES" ]]; then
  curl -G -L --fail --retry 3 \
    --data-urlencode "query@$step_dir/config/ensembl_release116_query.xml" \
    -o "$ORTHOLOGUES" \
    "https://www.ensembl.org/biomart/martservice"
fi

printf 'All Step 2 inputs are present under %s\n' "$STEP2_DATA_ROOT"
