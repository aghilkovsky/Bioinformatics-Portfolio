#!/usr/bin/env bash
# QIIME 2 workflow — 16S rRNA analysis for Propolis
# Author: Anastasia Ghilkovsky
# Usage: bash microbiome_16S/scripts/qiime2_commands.sh
# Requires: a QIIME 2 conda environment (e.g., 2024.10)

set -euo pipefail

############################
# 0) Activate QIIME 2 env  #
############################
# If you’re on your own laptop:
#   conda activate qiime2-amplicon-2024.10
# If you’re on a cluster with modules, uncomment and adjust:
#   module load miniconda3
#   conda activate qiime2-amplicon-2024.10

############################
# 1) Paths (relative)      #
############################
# Project root = repo folder
ROOT="microbiome_16S"
META="${ROOT}/sample-metadata.tsv"

# Inputs expected to already exist from your prior steps:
# - 02_dada2/table.qza
# - 03_taxonomy/taxonomy-reclassified.qza
# - 03_taxonomy/rep-seqs-v4.qza
# If you need to run DADA2/classifier first, keep those steps in a separate script.

TABLE_QZA="${ROOT}/02_dada2/table.qza"
REPK_QZA="${ROOT}/03_taxonomy/rep-seqs-v4.qza"
TAXA_QZA="${ROOT}/03_taxonomy/taxonomy-reclassified.qza"

# Outputs
TABLE_SUMMARY_QZV="${ROOT}/02_dada2/table-summary.qzv"
TABLE_NOCHLORO_QZA="${ROOT}/02_dada2/table-no-mitochondria-chloroplast.qza"
REPS_NOCHLORO_QZA="${ROOT}/03_taxonomy/rep-seqs-no-mitochondria-chloroplast.qza"

ALIGN_QZA="${ROOT}/04_tree/aligned-rep-seqs.qza"
MASK_ALIGN_QZA="${ROOT}/04_tree/masked-aligned-rep-seqs.qza"
TREE_UNROOT_QZA="${ROOT}/04_tree/unrooted-tree.qza"
TREE_ROOT_QZA="${ROOT}/04_tree/rooted-tree.qza"

CORE_DIR="${ROOT}/05_core_metrics"
SHANNON_QZV="${CORE_DIR}/shannon_group_significance.qzv"
BRAY_PERM_QZV="${CORE_DIR}/bray_curtis-permanova-Group.qzv"

TAXA_BAR_QZV="${ROOT}/03_taxonomy/taxa-bar-plots.qzv"

DIFF_DIR="${ROOT}/06_differential"
COMP_TABLE_QZA="${DIFF_DIR}/composition.qza"
ANCOM_QZV="${DIFF_DIR}/ancom-Group.qzv"

UNASSIGNED_IDS="${ROOT}/exported-taxonomy/unassigned-ids.txt"
UNASSIGNED_REP_QZA="${ROOT}/unassigned-rep-seqs.qza"
UNASSIGNED_DIR="${ROOT}/exported-unassigned"

############################
# 2) Create folders        #
############################
mkdir -p "${ROOT}/02_dada2" "${ROOT}/03_taxonomy" "${ROOT}/04_tree" "${CORE_DIR}" "${DIFF_DIR}" "${ROOT}/exported-taxonomy" "${UNASSIGNED_DIR}"

############################
# 3) Sanity checks         #
############################
if [[ ! -f "${META}" ]]; then
  echo "ERROR: Missing metadata file: ${META}"
  echo "Create it with a column named 'SampleID' and a grouping column 'Group'."
  exit 1
fi
for f in "${TABLE_QZA}" "${REPK_QZA}" "${TAXA_QZA}"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing input artifact: $f"; exit 1; }
done

#############################################
# 4) Summarize the feature table            #
#############################################
qiime feature-table summarize \
  --i-table "${TABLE_QZA}" \
  --o-visualization "${TABLE_SUMMARY_QZV}" \
  --m-sample-metadata-file "${META}"

echo "[OK] Feature table summary -> ${TABLE_SUMMARY_QZV}"

########################################################
# 5) Remove mitochondrial & chloroplast assignments     #
########################################################
qiime taxa filter-table \
  --i-table "${TABLE_QZA}" \
  --i-taxonomy "${TAXA_QZA}" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table "${TABLE_NOCHLORO_QZA}"

qiime taxa filter-seqs \
  --i-sequences "${REPK_QZA}" \
  --i-taxonomy "${TAXA_QZA}" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences "${REPS_NOCHLORO_QZA}"

echo "[OK] Filtered mito/chloro."

############################################
# 6) Phylogenetic tree (MAFFT + FastTree)  #
############################################
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${REPS_NOCHLORO_QZA}" \
  --o-alignment "${ALIGN_QZA}" \
  --o-masked-alignment "${MASK_ALIGN_QZA}" \
  --o-tree "${TREE_UNROOT_QZA}" \
  --o-rooted-tree "${TREE_ROOT_QZA}"

echo "[OK] Rooted tree -> ${TREE_ROOT_QZA}"

############################################
# 7) Core diversity metrics                #
############################################
# IMPORTANT: choose a sampling depth from your table summary rarefaction curve.
# For portfolio/demo, keep it small. Edit the value below as needed.
SAMPLING_DEPTH=465

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${TREE_ROOT_QZA}" \
  --i-table "${TABLE_NOCHLORO_QZA}" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --m-metadata-file "${META}" \
  --output-dir "${CORE_DIR}"

# Alpha diversity – Shannon
qiime diversity alpha-group-significance \
  --i-alpha-diversity "${CORE_DIR}/shannon_vector.qza" \
  --m-metadata-file "${META}" \
  --o-visualization "${SHANNON_QZV}"

# Beta diversity – Bray-Curtis PERMANOVA by 'Group' column in metadata
qiime diversity beta-group-significance \
  --i-distance-matrix "${CORE_DIR}/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "${META}" \
  --m-metadata-column Group \
  --o-visualization "${BRAY_PERM_QZV}" \
  --p-pairwise

echo "[OK] Core diversity + stats."

############################################
# 8) Taxa bar plot                         #
############################################
qiime taxa barplot \
  --i-table "${TABLE_NOCHLORO_QZA}" \
  --i-taxonomy "${TAXA_QZA}" \
  --m-metadata-file "${META}" \
  --o-visualization "${TAXA_BAR_QZV}"

echo "[OK] Taxa bar plot."

############################################
# 9) Differential abundance (ANCOM)        #
############################################
qiime composition add-pseudocount \
  --i-table "${TABLE_NOCHLORO_QZA}" \
  --o-composition-table "${COMP_TABLE_QZA}"

qiime composition ancom \
  --i-table "${COMP_TABLE_QZA}" \
  --m-metadata-file "${META}" \
  --m-metadata-column Group \
  --o-visualization "${ANCOM_QZV}"

echo "[OK] ANCOM completed."

############################################
# 10) Unassigned sequences → FASTA         #
############################################
# First, export taxonomy so we can grep unassigned IDs
qiime tools export \
  --input-path "${TAXA_QZA}" \
  --output-path "${ROOT}/exported-taxonomy"

# Collect feature IDs annotated as Unassigned
grep -P '\tUnassigned' "${ROOT}/exported-taxonomy/taxonomy.tsv" | cut -f1 > "${UNASSIGNED_IDS}"
sed -i.bak '1i#OTU ID' "${UNASSIGNED_IDS}" && rm -f "${UNASSIGNED_IDS}.bak"

# Filter unassigned sequences from the filtered rep-seqs
qiime feature-table filter-seqs \
  --i-data "${REPS_NOCHLORO_QZA}" \
  --m-metadata-file "${UNASSIGNED_IDS}" \
  --o-filtered-data "${UNASSIGNED_REP_QZA}"

# Export to FASTA (for BLAST)
qiime tools export \
  --input-path "${UNASSIGNED_REP_QZA}" \
  --output-path "${UNASSIGNED_DIR}"

echo "[OK] Unassigned sequences exported to ${UNASSIGNED_DIR} (for BLAST)."

echo "DONE. Open the .qzv files at https://view.qiime2.org/"
