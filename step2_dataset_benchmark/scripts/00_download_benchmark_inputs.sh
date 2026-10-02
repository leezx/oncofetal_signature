#!/usr/bin/env bash
set -euo pipefail

data_root="${BENCHMARK_DATA_ROOT:-/Volumes/Stelligen_SSD/Stelligen/DATA/2.PROJECTS/1.TWEAKR-oncoFetal/dataset_qualification_v0.1}"
mkdir -p "$data_root/raw/GSE101531" "$data_root/raw/GSE158702" \
  "$data_root/raw/GSE160449" "$data_root/raw/GSE44433" "$data_root/processed"

fetch() {
  local url="$1"
  local dest="$2"
  if [[ ! -s "$dest" ]]; then
    curl -fL "$url" -o "$dest"
  fi
}

fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE101nnn/GSE101531/suppl/GSE101531_RAW.tar" \
  "$data_root/raw/GSE101531/GSE101531_RAW.tar"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160449/suppl/GSE160449_org_ln_rnaseq_counts.txt.gz" \
  "$data_root/raw/GSE160449/GSE160449_org_ln_rnaseq_counts.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44433/matrix/GSE44433_series_matrix.txt.gz" \
  "$data_root/raw/GSE44433/GSE44433_series_matrix.txt.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz" \
  "$data_root/raw/Mus_musculus.gene_info.gz"

fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4808nnn/GSM4808339/suppl/GSM4808339_EPI1_RUN3.tar.gz" \
  "$data_root/raw/GSE158702/GSM4808339_EPI1_RUN3.tar.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4808nnn/GSM4808340/suppl/GSM4808340_EPI2_RUN3.tar.gz" \
  "$data_root/raw/GSE158702/GSM4808340_EPI2_RUN3.tar.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4808nnn/GSM4808341/suppl/GSM4808341_EPI3_RUN3.tar.gz" \
  "$data_root/raw/GSE158702/GSM4808341_EPI3_RUN3.tar.gz"
fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4808nnn/GSM4808345/suppl/GSM4808345_EPI_run2.tar.gz" \
  "$data_root/raw/GSE158702/GSM4808345_EPI_run2.tar.gz"

find "$data_root/raw" -type f -exec shasum -a 256 {} \; | sort > "$data_root/raw/checksums.sha256"
