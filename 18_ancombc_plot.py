"""
ANCOM-BC Results Visualization
================================
Reads exported ANCOM-BC differentials and creates a volcano-style plot
and a ranked bar chart of top differentially abundant taxa.

Run after: bash 18_ancombc.sh

Input:  version-2_integrated/ancombc_out/ancombc-exported/
Output: version-2_integrated/Main_Fig_ANCOMBC.png/.pdf
        version-2_integrated/ANCOMBC_significant.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import config

OUTPUT_DIR  = config.VERSION_DIR
ANCOMBC_DIR = OUTPUT_DIR / 'ancombc_out' / 'ancombc-exported'

# ─────────────────────────────────────────────
# Load taxonomy for labeling
# ─────────────────────────────────────────────

def get_genus(taxon_str):
    if pd.isna(taxon_str):
        return 'Unassigned'
    for part in str(taxon_str).split(';'):
        part = part.strip()
        if part.startswith('g__') and len(part) > 3:
            return part[3:]
    return 'Unknown'

# ─────────────────────────────────────────────
# Load ANCOM-BC output
# ─────────────────────────────────────────────

def load_ancombc():
    # QIIME2 exports split CSVs: lfc_slice.csv, q_val_slice.csv, etc.
    lfc_file = ANCOMBC_DIR / 'lfc_slice.csv'
    q_file   = ANCOMBC_DIR / 'q_val_slice.csv'
    if not lfc_file.exists():
        raise FileNotFoundError(f'{lfc_file} not found. Run 18_ancombc.sh first.')

    lfc = pd.read_csv(lfc_file, index_col=0)
    q   = pd.read_csv(q_file,   index_col=0)

    # Column of interest: Study_GroupTerrestrial_Soil
    # positive lfc = Terrestrial > Space; negative = Space > Terrestrial
    term_col = [c for c in lfc.columns if 'Terrestrial' in c or 'terrestrial' in c]
    if not term_col:
        term_col = [c for c in lfc.columns if c != '(Intercept)']
    term_col = term_col[0]
    print(f'  Term column: {term_col}')

    df = pd.DataFrame({
        'lfc':  lfc[term_col],
        'qval': q[term_col],
    })
    # Flip sign: positive = Space-enriched, negative = Terrestrial-enriched
    df['lfc'] = -df['lfc']
    df['sig'] = df['qval'] < 0.05
    print(f'  Loaded {len(df)} features, {df["sig"].sum()} significant (q<0.05)')
    return df

# ─────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────

def plot_ancombc(df, tax):
    # Add genus labels
    df = df.copy()
    if 'id' in df.columns:
        df.index = df['id']

    # Merge taxonomy
    df['genus'] = tax.reindex(df.index)['Taxon'].apply(
        lambda x: get_genus(x) if pd.notna(x) else df.index[df.index == x][0][:8]
    ) if 'Taxon' in tax.columns else df.index.map(lambda x: x[:12])

    # Identify column names (QIIME2 uses Study_Group[T.Terrestrial_Soil] style)
    # load_ancombc() already provides 'lfc', 'qval', 'sig' columns
    if 'lfc' not in df.columns:
        raise ValueError(f'Expected lfc column. Available: {list(df.columns)}')

    print(f'  Using pre-processed lfc/qval columns')

    # Save significant taxa
    sig = df[df['sig']].sort_values('lfc', ascending=False)
    sig.to_csv(OUTPUT_DIR / 'ANCOMBC_significant.csv')
    print(f'  Significant (q<0.05): {len(sig)} features')
    print(f'    Space-enriched (lfc>0): {sum(sig["lfc"]>0)}')
    print(f'    Terrestrial-enriched (lfc<0): {sum(sig["lfc"]<0)}')

    # ── Plot: top significant taxa ──
    top_n = 30
    top = pd.concat([
        sig[sig['lfc'] > 0].nlargest(top_n // 2, 'lfc'),
        sig[sig['lfc'] < 0].nsmallest(top_n // 2, 'lfc'),
    ]).sort_values('lfc')

    if top.empty:
        print('  No significant taxa to plot.')
        return

    fig, ax = plt.subplots(figsize=(10, max(6, len(top) * 0.42)))
    colors = ['#2980b9' if v > 0 else '#c0392b' for v in top['lfc']]
    ax.barh(range(len(top)), top['lfc'], color=colors, edgecolor='white', height=0.7)

    labels = top['genus'].tolist()
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=9, fontstyle='italic')
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('ANCOM-BC Log Fold-Change (Space Flight / Terrestrial Soil)', fontsize=10)
    ax.set_title(
        'ANCOM-BC Differential Abundance\nSpace Flight vs. Terrestrial Soil (q < 0.05)',
        fontsize=12, fontweight='bold'
    )

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor='#2980b9', label=f'Space-enriched (n={sum(sig["lfc"]>0)})'),
        Patch(facecolor='#c0392b', label=f'Terrestrial-enriched (n={sum(sig["lfc"]<0)})'),
    ], fontsize=9)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'Main_Fig_ANCOMBC.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'Main_Fig_ANCOMBC.pdf', bbox_inches='tight')
    plt.close()
    print(f'  Saved: Main_Fig_ANCOMBC.png')
    print(f'  Saved: Main_Fig_ANCOMBC.pdf')


if __name__ == '__main__':
    print('=' * 60)
    print('ANCOM-BC Results Visualization')
    print('=' * 60)

    tax = pd.read_csv(config.TAXONOMY_FILE, sep='\t', index_col=0)

    print('\n[1] Loading ANCOM-BC output...')
    df = load_ancombc()

    print('\n[2] Plotting...')
    plot_ancombc(df, tax)

    print('\nDone.')
