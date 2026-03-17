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
echo "ANCOM-BC Differential Abundance Analysis (Genus Level)"
echo "============================================================"

# ── Step 1: Filter metadata ──
echo "[1] Filtering metadata..."
python - <<'PYEOF'
import pandas as pd
meta = pd.read_csv('version-2_integrated/integrated_metadata.tsv', sep='\t', index_col=0)
meta = meta[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])]
meta.to_csv('version-2_integrated/ancombc_out/metadata_filtered.tsv', sep='\t')
print(f"  Kept {len(meta)} samples: Space={sum(meta.Study_Group=='Space_Flight')}, Terr={sum(meta.Study_Group=='Terrestrial_Soil')}")
PYEOF

# ── Step 2: Collapse to genus level in Python ──
echo "[2] Collapsing to genus level..."
python - <<'PYEOF'
import pandas as pd

meta  = pd.read_csv('version-2_integrated/ancombc_out/metadata_filtered.tsv', sep='\t', index_col=0)
table = pd.read_csv('version-2_integrated/exported_table_clean/feature-table.tsv',
                    sep='\t', skiprows=1, index_col=0)
tax   = pd.read_csv('version-2_integrated/exported_taxonomy/taxonomy.tsv', sep='\t', index_col=0)

# Filter samples
keep  = [c for c in table.columns if c in meta.index]
table = table[keep]

# Extract genus from SILVA taxonomy string
def get_genus(t):
    if pd.isna(t): return 'Unassigned'
    for p in str(t).split(';'):
        p = p.strip()
        if p.startswith('g__') and len(p) > 3:
            return p[3:]
    return 'Unassigned'

tax['Genus'] = tax['Taxon'].apply(get_genus)

# Collapse to genus
shared = table.index.intersection(tax.index)
table  = table.loc[shared].copy()
table['Genus'] = tax.loc[shared, 'Genus']
genus_table = table.groupby('Genus').sum()
genus_table = genus_table[genus_table.sum(axis=1) > 0]

# Remove Unassigned
genus_table = genus_table[genus_table.index != 'Unassigned']

print(f"  Genus table: {genus_table.shape[0]} genera × {genus_table.shape[1]} samples")

# Write with BIOM header
out = 'version-2_integrated/ancombc_out/genus-table.tsv'
with open(out, 'w') as f:
    f.write('# Constructed from biom file\n')
    genus_table.index.name = '#OTU ID'
    genus_table.to_csv(f, sep='\t')
print(f"  Saved: {out}")
PYEOF

# ── Step 3: Convert to BIOM and import ──
echo "[3] Converting to BIOM and importing to QIIME2..."
biom convert \
  -i "${OUTDIR}/genus-table.tsv" \
  -o "${OUTDIR}/genus-table.biom" \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path  "${OUTDIR}/genus-table.biom" \
  --output-path "${OUTDIR}/genus-table.qza"
echo "  Imported: genus-table.qza"

# ── Step 4: Run ANCOM-BC ──
echo "[4] Running ANCOM-BC (formula: Study_Group)..."
qiime composition ancombc \
  --i-table         "${OUTDIR}/genus-table.qza" \
  --m-metadata-file "${OUTDIR}/metadata_filtered.tsv" \
  --p-formula       "Study_Group" \
  --o-differentials "${OUTDIR}/ancombc-genus-differentials.qza" \
  --verbose
echo "  Done: ancombc-genus-differentials.qza"

# ── Step 5: Export ──
echo "[5] Exporting results..."
qiime composition tabulate \
  --i-data "${OUTDIR}/ancombc-genus-differentials.qza" \
  --o-visualization "${OUTDIR}/ancombc-genus-differentials.qzv"

rm -rf "${OUTDIR}/ancombc-genus-exported"
qiime tools export \
  --input-path  "${OUTDIR}/ancombc-genus-differentials.qza" \
  --output-path "${OUTDIR}/ancombc-genus-exported"
echo "  Exported to: ${OUTDIR}/ancombc-genus-exported/"

echo ""
echo "============================================================"
echo "ANCOM-BC (genus) complete. Run: python 18_ancombc_plot.py"
echo "============================================================"
