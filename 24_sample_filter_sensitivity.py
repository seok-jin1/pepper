"""
Sample Type Filter Sensitivity Analysis
=========================================
Critical check: Space_Flight includes 11 different sample types (fruit, leaf,
stem, seed, arcillite, foam, wick, swab, water, plant roots, rhizosphere).
The terrestrial control is rhizosphere only (n=20).

This script:
  1. Visualizes sample type composition
  2. Defines three filter levels for Space_Flight:
       ALL      : all 106 Space_Flight samples (current approach)
       ROOT_ZONE: Rhizosphere + Plant Roots + Arcillite (n=36)
                  → most defensible rhizosphere-equivalent in hydroponic ISS system
       STRICT   : Rhizosphere only (n=4)
  3. Compares alpha diversity, beta diversity (PCoA), and FAPROTAX nitrification
     across filter levels vs Terrestrial_Soil (n=20)
  4. Tests whether core conclusions hold when restricted to root-zone samples

Interpretation of arcillite:
  Arcillite (calcined clay aggregate) is the primary grow substrate in the
  VEGGIE plant growth hardware aboard ISS — functionally analogous to soil
  in terrestrial rhizosphere studies.

Output (version-2_integrated/):
  - SampleType_Composition.png/.pdf
  - FilterSensitivity_AlphaDiversity.png/.pdf
  - FilterSensitivity_BetaDiversity.png/.pdf
  - FilterSensitivity_Nitrification.png/.pdf
  - FilterSensitivity_Summary.csv
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, kruskal
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from sklearn.decomposition import PCA
import seaborn as sns

OUTPUT_DIR  = config.VERSION_DIR
RANDOM_SEED = config.RANDOM_SEED

# ── Filter definitions ───────────────────────────────────────────────────────
FILTER_LEVELS = {
    'ALL'       : ['Arcillite', 'Foam', 'Plant Roots', 'Rhizosphere',
                   'Swab', 'Water', 'Wick', 'fruit', 'leaf', 'seed', 'stem'],
    'ROOT_ZONE' : ['Rhizosphere', 'Plant Roots', 'Arcillite'],
    'STRICT'    : ['Rhizosphere'],
}
FILTER_COLORS = {
    'ALL'       : '#e74c3c',
    'ROOT_ZONE' : '#f39c12',
    'STRICT'    : '#27ae60',
    'Terrestrial_Soil': '#2ecc71',
}
FILTER_LABELS = {
    'ALL'       : 'All Space\n(n=106)',
    'ROOT_ZONE' : 'Root-Zone Space\n(Rhizosphere+Roots+Arcillite, n=36)',
    'STRICT'    : 'Strict Rhizosphere\nSpace (n=4)',
    'Terrestrial_Soil': 'Terrestrial\nRhizosphere (n=20)',
}

# ─────────────────────────────────────────────────────────────────────────────
def load_feature_table():
    table = pd.read_csv(config.FEATURE_TABLE_CLEAN, sep='\t',
                        skiprows=1, index_col=0)
    return table   # ASVs × samples

def relative_abundance(table):
    rel = table.div(table.sum(axis=0), axis=1) * 100
    return rel

def aggregate_to_genus(table_path, taxonomy_path):
    table = pd.read_csv(table_path, sep='\t', skiprows=1, index_col=0)
    tax   = pd.read_csv(taxonomy_path, sep='\t', index_col=0)
    def get_genus(s):
        if pd.isna(s): return 'Unassigned'
        for part in reversed(s.split(';')):
            part = part.strip()
            if part.startswith('g__'):
                name = part.split('__', 1)[-1].strip()
                return name if name else 'Unassigned'
        return 'Unassigned'
    tax['Genus'] = tax['Taxon'].apply(get_genus)
    shared = table.index.intersection(tax.index)
    table = table.loc[shared]
    table['Genus'] = tax.loc[shared, 'Genus']
    gt = table.groupby('Genus').sum()
    gt = gt.drop(index='Unassigned', errors='ignore')
    return gt.div(gt.sum(axis=0), axis=1) * 100

def bray_curtis(rel_table):
    """Bray-Curtis distance matrix on columns (samples)."""
    arr = rel_table.values.T  # (samples × ASVs)
    arr = arr / arr.sum(axis=1, keepdims=True) * 100
    n = arr.shape[0]
    dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            bc = np.sum(np.abs(arr[i] - arr[j])) / (arr[i].sum() + arr[j].sum())
            dm[i, j] = dm[j, i] = bc
    return dm

def bray_curtis_fast(rel_table):
    """Faster BC using scipy."""
    from scipy.spatial.distance import cdist
    arr = rel_table.values.T.astype(float)
    n = arr.shape[0]
    dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            num = np.sum(np.abs(arr[i] - arr[j]))
            denom = arr[i].sum() + arr[j].sum()
            dm[i, j] = dm[j, i] = num / denom if denom > 0 else 0
    return dm

def shannon_alpha(table):
    """Shannon diversity per sample."""
    rel = table.div(table.sum(axis=0), axis=1)
    return rel.apply(lambda s: -np.sum(s[s > 0] * np.log(s[s > 0])), axis=0)

# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('Sample Type Filter Sensitivity Analysis')
    print('=' * 60)

    # Load data
    meta  = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    table = load_feature_table()
    faprotax_func = pd.read_csv(OUTPUT_DIR / 'FAPROTAX_functional_table.tsv',
                                sep='\t', index_col=0)

    # ── Panel 1: Sample type composition ─────────────────────────────────────
    print('\n[1] Visualizing sample type composition...')

    sp_meta  = meta[meta['Study_Group'] == 'Space_Flight']
    tr_meta  = meta[meta['Study_Group'] == 'Terrestrial_Soil']
    type_counts = sp_meta['sample_type'].value_counts()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

    # Pie chart
    ax = axes[0]
    colors_pie = plt.cm.Set3(np.linspace(0, 1, len(type_counts)))
    wedges, texts, autotexts = ax.pie(
        type_counts.values, labels=type_counts.index,
        autopct='%1.0f%%', colors=colors_pie,
        pctdistance=0.75, startangle=90)
    for t in autotexts: t.set_fontsize(8)
    ax.set_title('Space_Flight Sample Type\nComposition (n=106)',
                 fontweight='bold', fontsize=12)

    # Bar chart with filter grouping
    ax = axes[1]
    categories = {'Root-Zone\n(Rhizosphere+\nRoots+Arcillite)':
                      ['Rhizosphere', 'Plant Roots', 'Arcillite'],
                  'Irrigation\nSystem': ['Foam', 'Wick', 'Water'],
                  'Plant\nTissue':      ['fruit', 'leaf', 'stem', 'seed'],
                  'Other':              ['Swab']}
    cat_counts = {cat: sum(type_counts.get(t, 0) for t in types)
                  for cat, types in categories.items()}
    ax.bar(range(len(cat_counts)), cat_counts.values(),
           color=['#2ecc71', '#3498db', '#e74c3c', '#95a5a6'], alpha=0.8)
    ax.set_xticks(range(len(cat_counts)))
    ax.set_xticklabels(list(cat_counts.keys()), fontsize=9)
    ax.set_ylabel('Number of Samples', fontsize=11)
    ax.set_title('Space_Flight by Ecological Category\n'
                 '(Terrestrial Soil: 20 Rhizosphere samples)',
                 fontweight='bold', fontsize=12)
    for i, (cat, cnt) in enumerate(cat_counts.items()):
        ax.text(i, cnt + 0.3, str(cnt), ha='center', fontsize=10, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    fig.suptitle('OSD-772 Space_Flight Dataset: Sample Type Heterogeneity',
                 fontsize=13, fontweight='bold')
    for ext in ['png', 'pdf']:
        fig.savefig(OUTPUT_DIR / f'SampleType_Composition.{ext}',
                    dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close()
    print('  Saved: SampleType_Composition.png/.pdf')

    print(f'\n  Filter level sample counts:')
    for level, types in FILTER_LEVELS.items():
        n = sum(sp_meta['sample_type'].isin(types))
        print(f'    {level:12s}: Space n={n}, Terrestrial n=20')

    # ── Panel 2: Alpha diversity sensitivity ─────────────────────────────────
    print('\n[2] Alpha diversity sensitivity...')

    shannon = shannon_alpha(table)
    results_alpha = {}

    # Terrestrial
    tr_samples = tr_meta.index.intersection(table.columns)
    results_alpha['Terrestrial_Soil'] = shannon[tr_samples].values

    # Each Space filter level
    for level, types in FILTER_LEVELS.items():
        sp_samples = sp_meta[sp_meta['sample_type'].isin(types)].index
        sp_samples = sp_samples.intersection(table.columns)
        results_alpha[level] = shannon[sp_samples].values

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    positions = list(range(len(results_alpha)))
    labels    = []
    for i, (key, vals) in enumerate(results_alpha.items()):
        col = FILTER_COLORS.get(key, '#95a5a6')
        bp  = ax.boxplot(vals, positions=[i], widths=0.5,
                         patch_artist=True, notch=False,
                         boxprops=dict(facecolor=col, alpha=0.7),
                         medianprops=dict(color='black', lw=2),
                         flierprops=dict(marker='o', ms=3, alpha=0.4))
        # MWU vs Terrestrial
        if key != 'Terrestrial_Soil':
            _, p = mannwhitneyu(vals, results_alpha['Terrestrial_Soil'],
                                alternative='two-sided')
            sig = '***' if p < 0.001 else ('**' if p < 0.01 else
                  ('*' if p < 0.05 else 'n.s.'))
            y_top = np.percentile(vals, 90) + 0.3
            ax.text(i, y_top, sig, ha='center', fontsize=11, fontweight='bold')
        n = len(vals)
        labels.append(FILTER_LABELS.get(key, key).replace('\n', '\n') + f'\nn={n}')

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Shannon Diversity Index', fontsize=11)
    ax.set_title('Alpha Diversity Across Filter Levels\n'
                 '(stars = vs Terrestrial_Soil, MWU p-value)',
                 fontweight='bold', fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    for ext in ['png', 'pdf']:
        fig.savefig(OUTPUT_DIR / f'FilterSensitivity_AlphaDiversity.{ext}',
                    dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close()
    print('  Saved: FilterSensitivity_AlphaDiversity.png/.pdf')

    # ── Panel 3: Beta diversity PCoA sensitivity ─────────────────────────────
    print('\n[3] Beta diversity PCoA sensitivity...')

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    for ax_idx, (level, types) in enumerate(FILTER_LEVELS.items()):
        ax = axes[ax_idx]
        sp_samples = sp_meta[sp_meta['sample_type'].isin(types)].index
        sp_samples = sp_samples.intersection(table.columns)
        tr_samples = tr_meta.index.intersection(table.columns)
        all_samples = list(sp_samples) + list(tr_samples)

        sub_table = table[all_samples]
        sub_rel   = sub_table.div(sub_table.sum(axis=0), axis=1) * 100
        arr = sub_rel.values.T.astype(float)

        # Fast BC PCoA
        bc = np.zeros((len(all_samples), len(all_samples)))
        for i in range(len(all_samples)):
            for j in range(i+1, len(all_samples)):
                s1, s2 = arr[i], arr[j]
                d = np.sum(np.abs(s1-s2)) / (s1.sum()+s2.sum()) if (s1.sum()+s2.sum())>0 else 0
                bc[i,j] = bc[j,i] = d

        # Double-center for PCoA
        n = len(all_samples)
        A = -0.5 * bc**2
        row_mean = A.mean(axis=1, keepdims=True)
        col_mean = A.mean(axis=0, keepdims=True)
        total_mean = A.mean()
        G = A - row_mean - col_mean + total_mean
        eigvals, eigvecs = np.linalg.eigh(G)
        idx = np.argsort(eigvals)[::-1]
        eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
        pos_idx = eigvals > 0
        coords = eigvecs[:, pos_idx] * np.sqrt(eigvals[pos_idx])
        pct = (eigvals[pos_idx] / eigvals[pos_idx].sum() * 100)

        sp_idx_arr = np.arange(len(sp_samples))
        tr_idx_arr = np.arange(len(sp_samples), len(all_samples))

        ax.scatter(coords[sp_idx_arr, 0], coords[sp_idx_arr, 1],
                   c='#3498db', alpha=0.6, s=40, label=f'Space ({len(sp_samples)})')
        ax.scatter(coords[tr_idx_arr, 0], coords[tr_idx_arr, 1],
                   c='#2ecc71', alpha=0.8, s=60, marker='^',
                   label=f'Terrestrial (20)')
        ax.set_xlabel(f'PC1 ({pct[0]:.1f}%)', fontsize=10)
        ax.set_ylabel(f'PC2 ({pct[1]:.1f}%)', fontsize=10)
        ax.set_title(f'{level}\nSpace filter', fontweight='bold', fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.2)

    fig.suptitle('Bray-Curtis PCoA Across Filter Levels',
                 fontsize=13, fontweight='bold')
    for ext in ['png', 'pdf']:
        fig.savefig(OUTPUT_DIR / f'FilterSensitivity_BetaDiversity.{ext}',
                    dpi=300 if ext == 'png' else None, bbox_inches='tight')
    plt.close()
    print('  Saved: FilterSensitivity_BetaDiversity.png/.pdf')

    # ── Panel 4: FAPROTAX nitrification sensitivity ───────────────────────────
    print('\n[4] FAPROTAX nitrification sensitivity...')

    if 'nitrification' in faprotax_func.index:
        nit_vals = {}
        tr_nit = faprotax_func.loc['nitrification', tr_meta.index.intersection(faprotax_func.columns)]
        nit_vals['Terrestrial_Soil'] = tr_nit.values

        for level, types in FILTER_LEVELS.items():
            sp_samples = sp_meta[sp_meta['sample_type'].isin(types)].index
            sp_samples = sp_samples.intersection(faprotax_func.columns)
            nit_vals[level] = faprotax_func.loc['nitrification', sp_samples].values

        fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
        for i, (key, vals) in enumerate(nit_vals.items()):
            col = FILTER_COLORS.get(key, '#95a5a6')
            ax.boxplot(vals, positions=[i], widths=0.5,
                       patch_artist=True, notch=False,
                       boxprops=dict(facecolor=col, alpha=0.7),
                       medianprops=dict(color='black', lw=2),
                       flierprops=dict(marker='o', ms=3, alpha=0.4))
            if key != 'Terrestrial_Soil':
                _, p = mannwhitneyu(vals, nit_vals['Terrestrial_Soil'],
                                    alternative='two-sided')
                sig = '***' if p < 0.001 else ('**' if p < 0.01 else
                      ('*' if p < 0.05 else 'n.s.'))
                ax.text(i, np.percentile(vals, 95) + 0.05,
                        sig, ha='center', fontsize=11, fontweight='bold')
            n = len(vals)
            ax.set_title(f'Nitrification Abundance Across Filter Levels',
                         fontweight='bold', fontsize=12)

        labels_nit = [FILTER_LABELS.get(k, k) + f'\nn={len(v)}'
                      for k, v in nit_vals.items()]
        ax.set_xticks(range(len(nit_vals)))
        ax.set_xticklabels(labels_nit, fontsize=9)
        ax.set_ylabel('Nitrification Relative Abundance (%)', fontsize=11)
        ax.grid(axis='y', alpha=0.3)
        for ext in ['png', 'pdf']:
            fig.savefig(OUTPUT_DIR / f'FilterSensitivity_Nitrification.{ext}',
                        dpi=300 if ext == 'png' else None, bbox_inches='tight')
        plt.close()
        print('  Saved: FilterSensitivity_Nitrification.png/.pdf')

    # ── Summary table ─────────────────────────────────────────────────────────
    print('\n[5] Generating summary table...')
    rows = []
    tr_sh = results_alpha['Terrestrial_Soil']
    tr_nit = nit_vals.get('Terrestrial_Soil', np.array([np.nan]))
    for level in list(FILTER_LEVELS.keys()) + ['Terrestrial_Soil']:
        sp_samples = (sp_meta[sp_meta['sample_type'].isin(FILTER_LEVELS[level])].index
                      if level in FILTER_LEVELS else tr_meta.index)
        sp_samples = sp_samples.intersection(table.columns)
        sh_vals = results_alpha.get(level, np.array([np.nan]))
        nit = nit_vals.get(level, np.array([np.nan]))
        if level != 'Terrestrial_Soil' and len(sh_vals) > 0:
            _, p_sh = mannwhitneyu(sh_vals, tr_sh, alternative='two-sided')
            _, p_nit = mannwhitneyu(nit, tr_nit, alternative='two-sided')
        else:
            p_sh, p_nit = np.nan, np.nan
        rows.append({
            'Filter': level,
            'n_Space': len(sp_samples),
            'Shannon_Space_Mean': round(sh_vals.mean(), 3) if len(sh_vals) > 0 else np.nan,
            'Shannon_Space_SD':   round(sh_vals.std(), 3)  if len(sh_vals) > 0 else np.nan,
            'Shannon_vs_Terrestrial_p': round(p_sh, 4) if not np.isnan(p_sh) else np.nan,
            'Nitrification_Space_Mean': round(nit.mean(), 4) if len(nit) > 0 else np.nan,
            'Nitrification_vs_Terrestrial_p': round(p_nit, 4) if not np.isnan(p_nit) else np.nan,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUTPUT_DIR / 'FilterSensitivity_Summary.csv', index=False)
    print('  Saved: FilterSensitivity_Summary.csv')
    print('\n' + summary.to_string(index=False))
    print('\nDone.')
