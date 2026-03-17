#!/bin/bash
# =============================================================================
# ANCOM-BC Differential Abundance Analysis
# ISS Chile Pepper Rhizosphere Microbiome Study
#
# Identifies taxa that are statistically differentially abundant between
# Space Flight and Terrestrial Soil (q < 0.05, Holm-Bonferroni corrected).
#
# Prerequisites:
#   - version-2_integrated/exported_table_clean/feature-table.tsv
#   - version-2_integrated/integrated_metadata.tsv
#
# Run:
#   conda activate qiime2-amplicon
#   bash 18_ancombc.sh
#
# Output (version-2_integrated/ancombc_out/):
#   - differentials.qza / differentials.qzv
#   - ancombc_differentials.tsv    (exported)
# Then visualized by 18_ancombc_plot.py
# =============================================================================

set -euo pipefail
WORKDIR="version-2_integrated"
OUTDIR="${WORKDIR}/ancombc_out"
mkdir -p "${OUTDIR}"

echo "============================================================"
echo "ANCOM-BC Differential Abundance Analysis"
echo "============================================================"

# ── Step 1: Filter metadata to Space_Flight and Terrestrial_Soil ──
echo "[1] Filtering metadata..."
python - <<'PYEOF'
import pandas as pd
meta = pd.read_csv('version-2_integrated/integrated_metadata.tsv', sep='\t', index_col=0)
meta = meta[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])]
meta.to_csv('version-2_integrated/ancombc_out/metadata_filtered.tsv', sep='\t')
print(f"  Kept {len(meta)} samples: Space={sum(meta.Study_Group=='Space_Flight')}, Terr={sum(meta.Study_Group=='Terrestrial_Soil')}")
PYEOF

# ── Step 2: Filter feature table and convert to BIOM ──
echo "[2] Converting feature table to BIOM..."
python - <<'PYEOF'
import pandas as pd, subprocess, os
meta = pd.read_csv('version-2_integrated/ancombc_out/metadata_filtered.tsv', sep='\t', index_col=0)
table = pd.read_csv('version-2_integrated/exported_table_clean/feature-table.tsv',
                    sep='\t', skiprows=1, index_col=0)
keep = [c for c in table.columns if c in meta.index]
table = table[keep]
# Remove ASVs with zero total counts after filtering
table = table[table.sum(axis=1) > 0]
print(f"  Table: {table.shape[0]} ASVs × {table.shape[1]} samples")

# Write TSV with BIOM header
out = 'version-2_integrated/ancombc_out/feature-table-filtered.tsv'
with open(out, 'w') as f:
    f.write('# Constructed from biom file\n')
    table.index.name = '#OTU ID'
    table.to_csv(f, sep='\t')
print(f"  Saved TSV: {out}")
PYEOF

# Convert TSV to BIOM HDF5
biom convert \
  -i "${OUTDIR}/feature-table-filtered.tsv" \
  -o "${OUTDIR}/feature-table-filtered.biom" \
  --table-type="OTU table" \
  --to-hdf5
echo "  Converted to BIOM HDF5"

# ── Step 3: Import BIOM to QIIME2 ──
echo "[3] Importing to QIIME2..."
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path  "${OUTDIR}/feature-table-filtered.biom" \
  --output-path "${OUTDIR}/feature-table-filtered.qza"
echo "  Imported: feature-table-filtered.qza"

# ── Step 4: Run ANCOM-BC ──
echo "[4] Running ANCOM-BC (formula: Study_Group)..."
qiime composition ancombc \
  --i-table    "${OUTDIR}/feature-table-filtered.qza" \
  --m-metadata-file "${OUTDIR}/metadata_filtered.tsv" \
  --p-formula  "Study_Group" \
  --o-differentials "${OUTDIR}/ancombc-differentials.qza" \
  --verbose
echo "  Done: ancombc-differentials.qza"

# ── Step 5: Tabulate and export ──
echo "[5] Exporting results..."
qiime composition tabulate \
  --i-data "${OUTDIR}/ancombc-differentials.qza" \
  --o-visualization "${OUTDIR}/ancombc-differentials.qzv"

qiime tools export \
  --input-path  "${OUTDIR}/ancombc-differentials.qza" \
  --output-path "${OUTDIR}/ancombc-exported"
echo "  Exported to: ${OUTDIR}/ancombc-exported/"

echo ""
echo "============================================================"
echo "ANCOM-BC complete. Run: python 18_ancombc_plot.py"
echo "============================================================"
