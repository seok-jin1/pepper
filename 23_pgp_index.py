"""
Plant Growth-Promoting (PGP) Functional Index
===============================================
Calculates a composite index of plant growth-promoting (PGP) bacterial taxa
per sample, comparing Space Flight vs. Terrestrial Soil.

PGP index directly links microbial compositional shifts (ANCOM-BC, FAPROTAX)
to plant-relevant functional consequences, supporting the narrative for a
Plant Physiology audience.

PGP categories and genera (literature-curated):
  N-fixation    : Bradyrhizobium, Mesorhizobium, Rhizobium, Sinorhizobium,
                  Azorhizobium, Azospirillum, Herbaspirillum, Gluconacetobacter,
                  Azoarcus, Frankia, Cupriavidus, Burkholderia-Caballeronia-Paraburkholderia
  P-solubilization: Pseudomonas, Bacillus, Paenibacillus, Burkholderia-Caballeronia-Paraburkholderia,
                  Serratia, Enterobacter, Chryseobacterium, Alcaligenes, Rahnella
  IAA/phytohormone: Azospirillum, Pseudomonas, Bacillus, Sphingomonas,
                  Herbaspirillum, Streptomyces, Bradyrhizobium, Mesorhizobium
  Biocontrol    : Bacillus, Pseudomonas, Streptomyces, Lysobacter,
                  Stenotrophomonas, Serratia, Paenibacillus
  ACC-deaminase : Pseudomonas, Mesorhizobium, Bradyrhizobium, Burkholderia-Caballeronia-Paraburkholderia,
                  Alcaligenes, Variovorax, Cupriavidus

References:
  - Compant et al. (2019) Soil Biol Biochem
  - Verbon & Liberman (2016) Front Plant Sci
  - Glick (2014) Scientifica

Prerequisites:
  - version-2_integrated/exported_table_clean/feature-table.tsv
  - version-2_integrated/exported_taxonomy/taxonomy.tsv
  - version-2_integrated/integrated_metadata.tsv
  - version-2_integrated/ANCOMBC_significant.csv   (optional, for cross-ref)

Run:
  python 23_pgp_index.py

Output (version-2_integrated/):
  - PGP_Index_Results.csv       : per-sample PGP indices
  - PGP_Index_Summary.csv       : group statistics + Wilcoxon tests
  - Main_Fig_PGP_Index.png/.pdf
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

OUTPUT_DIR   = config.VERSION_DIR
RANDOM_SEED  = config.RANDOM_SEED

# ── PGP genera by functional category ────────────────────────────────────────
PGP_CATEGORIES = {
    'N-Fixation': [
        'Bradyrhizobium', 'Mesorhizobium', 'Rhizobium', 'Sinorhizobium',
        'Azorhizobium', 'Azospirillum', 'Herbaspirillum', 'Gluconacetobacter',
        'Azoarcus', 'Frankia', 'Cupriavidus',
        'Burkholderia-Caballeronia-Paraburkholderia',
        'Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium',
    ],
    'P-Solubilization': [
        'Pseudomonas', 'Bacillus', 'Paenibacillus',
        'Burkholderia-Caballeronia-Paraburkholderia',
        'Serratia', 'Enterobacter', 'Chryseobacterium',
        'Alcaligenes', 'Rahnella', 'Gluconobacter',
        'Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium',
    ],
    'IAA/Phytohormone': [
        'Azospirillum', 'Pseudomonas', 'Bacillus', 'Sphingomonas',
        'Herbaspirillum', 'Streptomyces', 'Bradyrhizobium', 'Mesorhizobium',
        'Gluconobacter', 'Arthrobacter',
    ],
    'Biocontrol': [
        'Bacillus', 'Pseudomonas', 'Streptomyces', 'Lysobacter',
        'Stenotrophomonas', 'Serratia', 'Paenibacillus',
        'Burkholderia-Caballeronia-Paraburkholderia',
    ],
    'ACC-Deaminase': [
        'Pseudomonas', 'Mesorhizobium', 'Bradyrhizobium',
        'Burkholderia-Caballeronia-Paraburkholderia',
        'Variovorax', 'Cupriavidus', 'Alcaligenes',
        'Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium',
    ],
}

# All unique PGP genera
ALL_PGP = sorted(set(g for gs in PGP_CATEGORIES.values() for g in gs))

CATEGORY_COLORS = {
    'N-Fixation':      '#27ae60',
    'P-Solubilization': '#f39c12',
    'IAA/Phytohormone': '#8e44ad',
    'Biocontrol':       '#e74c3c',
    'ACC-Deaminase':    '#2980b9',
}

GROUP_COLORS = {'Space_Flight': '#3498db', 'Terrestrial_Soil': '#2ecc71'}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Build genus-level feature table
# ─────────────────────────────────────────────────────────────────────────────

def aggregate_to_genus(feature_table_path, taxonomy_path):
    table = pd.read_csv(feature_table_path, sep='\t', skiprows=1, index_col=0)
    tax   = pd.read_csv(taxonomy_path, sep='\t', index_col=0)

    def get_genus(s):
        if pd.isna(s):
            return 'Unassigned'
        for part in reversed(s.split(';')):
            part = part.strip()
            if part.startswith('g__'):
                name = part.split('__', 1)[-1].strip()
                return name if name else 'Unassigned'
        return 'Unassigned'

    tax['Genus'] = tax['Taxon'].apply(get_genus)
    shared = table.index.intersection(tax.index)
    table  = table.loc[shared]
    table['Genus'] = tax.loc[shared, 'Genus']
    genus_table = table.groupby('Genus').sum()
    genus_table = genus_table.drop(index='Unassigned', errors='ignore')
    return genus_table.div(genus_table.sum(axis=0), axis=1) * 100  # rel. abund (%)


# ─────────────────────────────────────────────────────────────────────────────
# 2. Compute PGP indices
# ─────────────────────────────────────────────────────────────────────────────

def compute_pgp_indices(genus_rel, meta):
    """
    Returns DataFrame: index=samples, columns=PGP categories + 'Total_PGP'
    Value = sum of relative abundances of PGP genera in that category (%).
    """
    samples = genus_rel.columns.intersection(meta.index)
    genus_rel = genus_rel[samples]
    records   = {}

    # Category sub-indices
    for cat, genera in PGP_CATEGORIES.items():
        present = genus_rel.index.intersection(genera)
        if len(present) > 0:
            records[cat] = genus_rel.loc[present].sum(axis=0)
        else:
            records[cat] = pd.Series(0.0, index=genus_rel.columns)

    # Total PGP index (union, no double-counting)
    all_present = genus_rel.index.intersection(ALL_PGP)
    records['Total_PGP'] = genus_rel.loc[all_present].sum(axis=0) if len(all_present) > 0 \
                            else pd.Series(0.0, index=genus_rel.columns)

    # Individual PGP genus abundances (for stacked bar)
    for genus in ALL_PGP:
        if genus in genus_rel.index:
            records[f'genus_{genus}'] = genus_rel.loc[genus]
        else:
            records[f'genus_{genus}'] = pd.Series(0.0, index=genus_rel.columns)

    df = pd.DataFrame(records).T   # (categories × samples)
    return df.T                    # (samples × categories)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Statistical comparison
# ─────────────────────────────────────────────────────────────────────────────

def compare_pgp(pgp_df, meta):
    sp_idx = meta.index[meta['Study_Group'] == 'Space_Flight']
    tr_idx = meta.index[meta['Study_Group'] == 'Terrestrial_Soil']
    categories = list(PGP_CATEGORIES.keys()) + ['Total_PGP']

    rows   = []
    p_vals = []
    for cat in categories:
        sp = pgp_df.loc[pgp_df.index.intersection(sp_idx), cat].dropna()
        tr = pgp_df.loc[pgp_df.index.intersection(tr_idx), cat].dropna()
        stat, p = mannwhitneyu(sp, tr, alternative='two-sided')
        rows.append({
            'Category':  cat,
            'Space_Mean':  round(sp.mean(), 4),
            'Space_SD':    round(sp.std(), 4),
            'Terr_Mean':   round(tr.mean(), 4),
            'Terr_SD':     round(tr.std(), 4),
            'log2FC':      round(np.log2((sp.mean() + 1e-6) / (tr.mean() + 1e-6)), 4),
            'MWU_stat':    round(stat, 1),
            'p_raw':       round(p, 6),
        })
        p_vals.append(p)

    df = pd.DataFrame(rows)
    _, q_vals, _, _ = multipletests(p_vals, method='fdr_bh')
    df['q_BH'] = np.round(q_vals, 6)
    df['sig']  = df['q_BH'].apply(
        lambda q: '***' if q < 0.001 else ('**' if q < 0.01
                   else ('*' if q < 0.05 else 'n.s.')))
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 4. Plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_pgp(pgp_df, meta, summary_df):
    sp_idx = meta.index[meta['Study_Group'] == 'Space_Flight']
    tr_idx = meta.index[meta['Study_Group'] == 'Terrestrial_Soil']
    cats   = list(PGP_CATEGORIES.keys())

    fig = plt.figure(figsize=(18, 12), constrained_layout=True)
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

    # ── Panel A: Total PGP Index boxplot ─────────────────────────────────────
    ax_total = fig.add_subplot(gs[0, 0])
    plot_data_total = pd.DataFrame({
        'PGP_Index': pd.concat([
            pgp_df.loc[pgp_df.index.intersection(sp_idx), 'Total_PGP'],
            pgp_df.loc[pgp_df.index.intersection(tr_idx), 'Total_PGP'],
        ]),
        'Group': (
            ['Space Flight'] * len(sp_idx.intersection(pgp_df.index)) +
            ['Terrestrial Soil'] * len(tr_idx.intersection(pgp_df.index))
        )
    })
    tot_row = summary_df[summary_df['Category'] == 'Total_PGP'].iloc[0]

    for grp, col in [('Space Flight', '#3498db'), ('Terrestrial Soil', '#2ecc71')]:
        vals = plot_data_total.loc[plot_data_total['Group'] == grp, 'PGP_Index']
        ax_total.boxplot(vals, positions=[['Space Flight', 'Terrestrial Soil'].index(grp)],
                         widths=0.4, patch_artist=True, notch=False,
                         boxprops=dict(facecolor=col, alpha=0.7),
                         medianprops=dict(color='black', linewidth=2),
                         flierprops=dict(marker='o', markersize=3, alpha=0.4))

    y_max = plot_data_total['PGP_Index'].max()
    ax_total.annotate(
        f"p = {tot_row['p_raw']:.3f} {tot_row['sig']}",
        xy=(0.5, 0.92), xycoords='axes fraction', ha='center', fontsize=10
    )
    ax_total.set_xticks([0, 1])
    ax_total.set_xticklabels(['Space\nFlight', 'Terrestrial\nSoil'], fontsize=10)
    ax_total.set_ylabel('Total PGP Index (%)', fontsize=11)
    ax_total.set_title('Total PGP Index', fontsize=11, fontweight='bold')
    ax_total.grid(axis='y', alpha=0.3)

    # ── Panel B: Category sub-indices bar chart ───────────────────────────────
    ax_bar = fig.add_subplot(gs[0, 1:])
    x      = np.arange(len(cats))
    width  = 0.35

    sp_means = [pgp_df.loc[pgp_df.index.intersection(sp_idx), c].mean() for c in cats]
    tr_means = [pgp_df.loc[pgp_df.index.intersection(tr_idx), c].mean() for c in cats]
    sp_se    = [pgp_df.loc[pgp_df.index.intersection(sp_idx), c].sem()  for c in cats]
    tr_se    = [pgp_df.loc[pgp_df.index.intersection(tr_idx), c].sem()  for c in cats]

    ax_bar.bar(x - width/2, sp_means, width, yerr=sp_se,
               color='#3498db', alpha=0.8, label='Space Flight',
               capsize=4, error_kw={'elinewidth': 1.2})
    ax_bar.bar(x + width/2, tr_means, width, yerr=tr_se,
               color='#2ecc71', alpha=0.8, label='Terrestrial Soil',
               capsize=4, error_kw={'elinewidth': 1.2})

    for i, cat in enumerate(cats):
        row = summary_df[summary_df['Category'] == cat]
        if len(row) == 0:
            continue
        sig = row.iloc[0]['sig']
        if sig != 'n.s.':
            y_top = max(sp_means[i] + sp_se[i], tr_means[i] + tr_se[i]) + 0.1
            ax_bar.text(x[i], y_top, sig, ha='center', va='bottom',
                        fontsize=12, fontweight='bold')

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(cats, fontsize=10)
    ax_bar.set_ylabel('Mean Relative Abundance (%, mean ± SE)', fontsize=11)
    ax_bar.set_title('PGP Sub-Index by Category', fontsize=11, fontweight='bold')
    ax_bar.legend(fontsize=10)
    ax_bar.grid(axis='y', alpha=0.3)
    ax_bar.set_ylim(bottom=0)

    # ── Panel C: Stacked bar — individual PGP genus contributions ────────────
    ax_stack = fig.add_subplot(gs[1, :])

    genus_cols = [c for c in pgp_df.columns if c.startswith('genus_')]
    genus_names = [c.replace('genus_', '') for c in genus_cols]

    sp_genus_means = pgp_df.loc[pgp_df.index.intersection(sp_idx), genus_cols].mean()
    tr_genus_means = pgp_df.loc[pgp_df.index.intersection(tr_idx), genus_cols].mean()

    # Only show genera with mean RA > 0.01% in at least one group
    mask = (sp_genus_means > 0.01) | (tr_genus_means > 0.01)
    sp_genus_means = sp_genus_means[mask]
    tr_genus_means = tr_genus_means[mask]
    genus_names_filt = [n for n, m in zip(genus_names, mask) if m]

    # Sort by terrestrial mean for readability
    order = tr_genus_means.sort_values(ascending=False).index
    sp_vals = sp_genus_means[order].values
    tr_vals = tr_genus_means[order].values
    gnames  = [c.replace('genus_', '') for c in order]

    # Color by category membership
    def genus_to_color(g):
        # Prioritize first category match
        priority = ['N-Fixation', 'P-Solubilization', 'IAA/Phytohormone',
                    'Biocontrol', 'ACC-Deaminase']
        for cat in priority:
            if g in PGP_CATEGORIES[cat]:
                return CATEGORY_COLORS[cat]
        return 'gray'

    bar_colors = [genus_to_color(g) for g in gnames]

    x2 = np.arange(len(gnames))
    w2 = 0.35
    for xi, (sv, tv, col) in enumerate(zip(sp_vals, tr_vals, bar_colors)):
        ax_stack.bar(xi - w2/2, sv, w2, color=col, alpha=0.65,
                     edgecolor='white', linewidth=0.5)
        ax_stack.bar(xi + w2/2, tv, w2, color=col, alpha=0.95,
                     edgecolor='white', linewidth=0.5)

    # Custom legend for groups
    from matplotlib.patches import Patch
    group_patches = [
        Patch(facecolor='gray', alpha=0.65, label='Space Flight'),
        Patch(facecolor='gray', alpha=0.95, label='Terrestrial Soil'),
    ]
    cat_patches = [
        Patch(facecolor=CATEGORY_COLORS[c], label=c)
        for c in CATEGORY_COLORS
    ]
    ax_stack.legend(handles=group_patches + cat_patches,
                    ncol=3, fontsize=8, loc='upper right',
                    framealpha=0.8)

    ax_stack.set_xticks(x2)
    ax_stack.set_xticklabels(gnames, rotation=40, ha='right', fontsize=9)
    ax_stack.set_ylabel('Mean Relative Abundance (%)', fontsize=11)
    ax_stack.set_title('Individual PGP Genus Contributions\n'
                        '(Left bar = Space Flight, Right bar = Terrestrial Soil; '
                        'color = PGP category)',
                        fontsize=11, fontweight='bold')
    ax_stack.grid(axis='y', alpha=0.3)
    ax_stack.set_ylim(bottom=0)

    fig.suptitle(
        'Plant Growth-Promoting (PGP) Bacterial Index\n'
        'Space Flight vs. Terrestrial Soil Rhizosphere',
        fontsize=14, fontweight='bold'
    )

    for ext in ['png', 'pdf']:
        out = OUTPUT_DIR / f'Main_Fig_PGP_Index.{ext}'
        fig.savefig(out, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f'  Saved: {out.name}')
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('Plant Growth-Promoting (PGP) Functional Index')
    print('=' * 60)

    print('\n[1] Building genus-level feature table...')
    genus_rel = aggregate_to_genus(config.FEATURE_TABLE_CLEAN, config.TAXONOMY_FILE)
    print(f'  {genus_rel.shape[0]} genera × {genus_rel.shape[1]} samples')

    print('\n[2] Loading metadata...')
    meta = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    meta = meta[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])]
    print(f'  {len(meta)} samples '
          f'(Space={sum(meta.Study_Group=="Space_Flight")}, '
          f'Terrestrial={sum(meta.Study_Group=="Terrestrial_Soil")})')

    print('\n[3] PGP genera coverage in dataset:')
    for genus in sorted(ALL_PGP):
        present = genus in genus_rel.index
        if present:
            sp_idx = meta.index[meta.Study_Group == 'Space_Flight']
            tr_idx = meta.index[meta.Study_Group == 'Terrestrial_Soil']
            sp_m = genus_rel.loc[genus, genus_rel.columns.intersection(sp_idx)].mean()
            tr_m = genus_rel.loc[genus, genus_rel.columns.intersection(tr_idx)].mean()
            if sp_m > 0.01 or tr_m > 0.01:
                print(f'  {genus:50s}  Space={sp_m:.3f}%  Terr={tr_m:.3f}%')

    print('\n[4] Computing PGP indices...')
    pgp_df = compute_pgp_indices(genus_rel, meta)
    pgp_df = pgp_df.loc[pgp_df.index.intersection(meta.index)]

    cats = list(PGP_CATEGORIES.keys()) + ['Total_PGP']
    sp_idx = meta.index[meta.Study_Group == 'Space_Flight']
    tr_idx = meta.index[meta.Study_Group == 'Terrestrial_Soil']
    print('\n  PGP index summary:')
    for cat in cats:
        sp_m = pgp_df.loc[pgp_df.index.intersection(sp_idx), cat].mean()
        tr_m = pgp_df.loc[pgp_df.index.intersection(tr_idx), cat].mean()
        print(f'  {cat:20s}  Space={sp_m:.3f}%  Terr={tr_m:.3f}%')

    print('\n[5] Statistical comparison...')
    summary_df = compare_pgp(pgp_df, meta)
    print(summary_df[['Category', 'Space_Mean', 'Terr_Mean',
                       'log2FC', 'p_raw', 'q_BH', 'sig']].to_string(index=False))

    print('\n[6] Cross-reference with ANCOM-BC significant taxa...')
    ancombc_path = OUTPUT_DIR / 'ANCOMBC_significant.csv'
    if ancombc_path.exists():
        ancombc = pd.read_csv(ancombc_path)
        ancombc_genera = set(ancombc['genus'].dropna())
        pgp_ancombc = [g for g in ALL_PGP if g in ancombc_genera]
        print(f'  PGP genera also significant in ANCOM-BC ({len(pgp_ancombc)}):')
        for g in pgp_ancombc:
            row = ancombc[ancombc['genus'] == g].iloc[0]
            direction = 'Space-enriched' if row['lfc'] > 0 else 'Terrestrial-enriched'
            print(f'    {g}: lfc={row["lfc"]:.2f}, q={row["qval"]:.2e} [{direction}]')
    else:
        print('  ANCOMBC_significant.csv not found — skipping cross-reference.')

    print('\n[7] Saving results...')
    out_df = pgp_df[cats].copy()
    out_df.insert(0, 'Study_Group', meta.loc[out_df.index, 'Study_Group'])
    out_df.to_csv(OUTPUT_DIR / 'PGP_Index_Results.csv')
    summary_df.to_csv(OUTPUT_DIR / 'PGP_Index_Summary.csv', index=False)
    print('  Saved: PGP_Index_Results.csv')
    print('  Saved: PGP_Index_Summary.csv')

    print('\n[8] Plotting...')
    plot_pgp(pgp_df, meta, summary_df)

    print('\nDone.')
