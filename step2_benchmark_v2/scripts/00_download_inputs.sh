#!/usr/bin/env bash
# Acquire the three new human fetal/adult resources for the Step 2 human
# benchmark v2 into their managed DATA dataset directories. Idempotent: an
# existing non-empty file is never re-downloaded or overwritten.
set -euo pipefail

DATA_ROOT="${STEP2V2_DATA_ROOT:-/Volumes/Stelligen_SSD/Stelligen/DATA}"
geo="https://ftp.ncbi.nlm.nih.gov/geo/series"

fetch() {  # url dest
  local url="$1" dest="$2"
  mkdir -p "$(dirname "$dest")"
  if [[ -s "$dest" ]]; then echo "exists  $dest"; return; fi
  curl -fsSL --retry 3 -o "$dest.part" "$url" && mv "$dest.part" "$dest"
  echo "fetched $dest"
}

# H-new1: Fawkner-Corbett 2021 spatial transcriptomics (Visium), fetal + adult.
fawk_st="$DATA_ROOT/SpatialTranscriptomics/GSE158328_FawknerCorbett2021/raw"
fetch "$geo/GSE158nnn/GSE158328/suppl/GSE158328_RAW.tar" "$fawk_st/GSE158328_RAW.tar"
fetch "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158328&targ=gsm&form=text&view=brief" \
  "$fawk_st/GSE158328_samples_brief.txt"

# H-new2 fetal: Fawkner-Corbett 2021 epithelial pools (GEX) + hashtag (HTO) libraries.
fawk_sc="$DATA_ROOT/scRNAseq/GSE158702_FawknerCorbett2021/raw"
for f in GSM4808339_EPI1_RUN3 GSM4808340_EPI2_RUN3 GSM4808341_EPI3_RUN3 GSM4808345_EPI_run2 \
         GSM4808349_HTO1 GSM4808350_HTO2 GSM4808351_HTO3 GSM4808352_HTO4 GSM4808353_HTO5 \
         GSM4808354_HTO6 GSM4808355_HTO_epi_4 GSM4808356_HTO_stromal_4 GSM4808357_Pool5_HTO \
         GSM4808358_Pool6_HTO; do
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4808nnn/${f%%_*}/suppl/$f.tar.gz" "$fawk_sc/$f.tar.gz"
done
fetch "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158702&targ=gsm&form=text&view=brief" \
  "$fawk_sc/GSE158702_samples_brief.txt"

# H-new2 adult: Burclaff 2022 healthy adult epithelium, 3 donors.
burc="$DATA_ROOT/scRNAseq/GSE185224_Burclaff2022/raw"
for f in GSE185224_Donor1_filtered_feature_bc_matrix.h5 GSE185224_Donor2_filtered_feature_bc_matrix.h5 \
         GSE185224_Donor3_filtered_feature_bc_matrix.h5 GSE185224_Donors_IntestinalRegions_Hashtags.txt.gz \
         GSE185224_clustered_annotated_adata_k10_lr0.92_v1.7.h5ad.gz; do
  fetch "$geo/GSE185nnn/GSE185224/suppl/$f" "$burc/$f"
done

# H-new3 adult: Wang 2020 adult ileum/colon/rectum epithelium.
wang="$DATA_ROOT/scRNAseq/GSE125970_Wang2020/raw"
for f in GSE125970_raw_UMIcounts.txt.gz GSE125970_cell_info.txt.gz; do
  fetch "$geo/GSE125nnn/GSE125970/suppl/$f" "$wang/$f"
done

for d in "$fawk_st" "$fawk_sc" "$burc" "$wang"; do
  (cd "$d" && shasum -a 256 ./*.gz ./*.tar ./*.h5 2>/dev/null > checksums.sha256 || true)
done
echo "Downloads complete."
