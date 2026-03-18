"""
Temporal Stability Analysis (Q1–Q4)
=====================================
Assesses whether the functional collapse signal (nitrification loss, reduced
functional redundancy) is consistent across all four quarters (Q1–Q4) of the
ISS experiment.

Key questions:
  1. Is nitrification abundance consistently low in Space across all quarters?
  2. Is N-cycle completeness consistently lower in Space across all quarters?
  3. Are Shannon diversity differences between Space and Terrestrial consistent
     across quarters? (checks for batch effects)

Method:
  - Extract quarter (Q1–Q4) from metadata column 'time_point' or 'quarter'
  - Per quarter: compute mean nitrification FAPROTAX abundance (Space samples)
  - Per quarter: compute N-cycle completeness (Space samples)
  - Per quarter: compare Shannon diversity (Space vs. Terrestrial)
  - Plot as multi-panel figure: one panel per metric across Q1–Q4

Prerequisites:
  - version-2_integrated/FAPROTAX_report.txt
  - version-2_integrated/exported_table_clean/feature-table.tsv
  - version-2_integrated/exported_taxonomy/taxonomy.tsv
  - version-2_integrated/integrated_metadata.tsv
  - version-2_integrated/exported_diversity/shannon_vector.tsv (or computed)

Run:
  python 28_temporal_stability.py

Output (version-2_integrated/):
  - TemporalStability_Results.csv
  - Main_Fig_TemporalStability.png/.pdf
"""

import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUTPUT_DIR = config.VERSION_DIR

COLORS = {'Space_Flight': '#3498db', 'Terrestrial_Soil': '#2ecc71'}
GROUP_LABELS = {'Space_Flight': 'Space Flight', 'Terrestrial_Soil': 'Terrestrial Soil'}

NCYCLE_STEPS = [
    'nitrogen_fixation', 'nitrification', 'denitrification',
    'nitrate_reduction', 'ureolysis',
]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def parse_faprotax_report(report_path):
    func2genera = {}
    current_func = None
    with open(report_path) as fh:
        for line in fh:
            line = line.rstrip()
            m = re.match(r'^# (\w+) \(\d+ records\):', line)
            if m:
                current_func = m.group(1)
                func2genera[current_func] = set()
                continue
            if current_func and line.startswith('    '):
                for part in line.split(';'):
                    part = part.strip()
                    if part.startswith('g__'):
                        name = part.split('__', 1)[-1].strip()
                        if name and name not in ('uncultured', 'metagenome',
                                                 'unidentified', ''):
                            func2genera[current_func].add(name)
    return func2genera


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
    return genus_table.div(genus_table.sum(axis=0), axis=1) * 100


def compute_shannon(genus_rel):
    """Compute Shannon H' per sample from relative abundance (%)."""
    p = genus_rel / 100.0
    p = p.clip(lower=0)
    h = -(p * np.log(p + 1e-12)).sum(axis=0)
    return h


def get_nitrification_abund(genus_rel, func2genera):
    """Sum relative abundance of nitrification genera per sample."""
    nitri_genera = func2genera.get('nitrification', set())
    present = genus_rel.index.intersection(nitri_genera)
    if len(present) == 0:
        return pd.Series(0.0, index=genus_rel.columns)
    return genus_rel.loc[present].sum(axis=0)


def compute_ncycle_completeness(genus_rel, func2genera):
    """Per sample: fraction of N-cycle steps with ≥1 genus present."""
    rows = {}
    for sample in genus_rel.columns:
        step_scores = []
        for step in NCYCLE_STEPS:
            assigned = func2genera.get(step, set())
            present  = assigned.intersection(genus_rel.index)
            if len(present) == 0:
                step_scores.append(0)
            else:
                step_scores.append(int((genus_rel.loc[list(present), sample] > 0).any()))
        rows[sample] = np.mean(step_scores)
    return pd.Series(rows)


def detect_quarter_column(meta):
    """Try to find a quarter/time_point column in metadata."""
    for col in ['quarter', 'Quarter', 'time_point', 'TimePoint',
                'time', 'Time', 'timepoint']:
        if col in meta.columns:
            return col
    return None


def assign_quarter_from_sampleid(meta):
    """
    Attempt to derive quarter from sample IDs.
    OSD-772 sample IDs typically contain 'Q1', 'Q2', 'Q3', 'Q4'.
    Returns a Series with Q1/Q2/Q3/Q4 or None.
    """
    quarters = {}
    for sid in meta.index:
        sid_str = str(sid).upper()
        found = None
        for q in ['Q1', 'Q2', 'Q3', 'Q4']:
            if q in sid_str:
                found = q
                break
        quarters[sid] = found
    return pd.Series(quarters)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('Temporal Stability Analysis (Q1–Q4)')
    print('=' * 60)

    print('\n[1] Loading data...')
    func2genera = parse_faprotax_report(OUTPUT_DIR / 'FAPROTAX_report.txt')
    genus_rel   = aggregate_to_genus(config.FEATURE_TABLE_CLEAN, config.TAXONOMY_FILE)
    meta        = pd.read_csv(config.INTEGRATED_METADATA, sep='\t', index_col=0)
    meta        = meta[meta['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])]
    genus_rel   = genus_rel[genus_rel.columns.intersection(meta.index)]
    print(f'  {genus_rel.shape[0]} genera × {genus_rel.shape[1]} samples')

    # ── Detect or derive quarter ─────────────────────────────────────────────
    print('\n[2] Detecting temporal grouping...')
    qcol = detect_quarter_column(meta)
    if qcol:
        print(f'  Found column: {qcol}')
        meta['Quarter'] = meta[qcol].astype(str)
    else:
        print('  No quarter column found — deriving from sample IDs...')
        meta['Quarter'] = assign_quarter_from_sampleid(meta)

    # Only Space_Flight has quarterly structure; Terrestrial is a single group
    sp_meta = meta[meta['Study_Group'] == 'Space_Flight'].copy()
    tr_meta = meta[meta['Study_Group'] == 'Terrestrial_Soil'].copy()

    quarters = [q for q in ['Q1', 'Q2', 'Q3', 'Q4']
                if q in sp_meta['Quarter'].values]
    if not quarters:
        # Try to use unique non-null quarter values
        quarters = sorted([q for q in sp_meta['Quarter'].dropna().unique()
                           if q not in ('None', 'nan', 'NaN')])
    print(f'  Quarters found: {quarters}')
    if not quarters:
        print('  WARNING: No quarterly structure detected. '
              'Will use full Space group as single time point.')
        quarters = ['ALL']
        sp_meta['Quarter'] = 'ALL'

    print(f'  Space samples per quarter:')
    for q in quarters:
        n = (sp_meta['Quarter'] == q).sum()
        print(f'    {q}: n={n}')

    # ── Compute per-sample metrics ───────────────────────────────────────────
    print('\n[3] Computing per-sample metrics...')
    nitri_abund = get_nitrification_abund(genus_rel, func2genera)
    ncycle_comp = compute_ncycle_completeness(genus_rel, func2genera)
    shannon_h   = compute_shannon(genus_rel)

    # ── Aggregate per quarter (Space) and overall (Terrestrial) ─────────────
    print('\n[4] Aggregating by quarter...')

    records = []
    for q in quarters:
        q_ids = sp_meta.index[sp_meta['Quarter'] == q]
        q_ids = q_ids.intersection(genus_rel.columns)
        if len(q_ids) == 0:
            continue

        n_q         = len(q_ids)
        nitri_mean  = nitri_abund[q_ids].mean()
        nitri_se    = nitri_abund[q_ids].sem()
        nc_mean     = ncycle_comp[q_ids].mean()
        nc_se       = ncycle_comp[q_ids].sem()
        shannon_mean = shannon_h[q_ids].mean()
        shannon_se   = shannon_h[q_ids].sem()

        records.append({
            'Quarter': q, 'Group': 'Space_Flight', 'n': n_q,
            'Nitri_Mean': round(nitri_mean, 4), 'Nitri_SE': round(nitri_se, 4),
            'NC_Mean': round(nc_mean, 4), 'NC_SE': round(nc_se, 4),
            'Shannon_Mean': round(shannon_mean, 4), 'Shannon_SE': round(shannon_se, 4),
        })
        print(f'  {q}: n={n_q}, Nitri={nitri_mean:.3f}%, NC={nc_mean:.3f}, '
              f'Shannon={shannon_mean:.3f}')

    # Terrestrial (single time point, shown as reference line)
    tr_ids = tr_meta.index.intersection(genus_rel.columns)
    tr_nitri  = nitri_abund[tr_ids]
    tr_nc     = ncycle_comp[tr_ids]
    tr_shannon = shannon_h[tr_ids]
    records.append({
        'Quarter': 'Terrestrial', 'Group': 'Terrestrial_Soil', 'n': len(tr_ids),
        'Nitri_Mean':   round(tr_nitri.mean(), 4),   'Nitri_SE':   round(tr_nitri.sem(), 4),
        'NC_Mean':      round(tr_nc.mean(), 4),       'NC_SE':      round(tr_nc.sem(), 4),
        'Shannon_Mean': round(tr_shannon.mean(), 4),  'Shannon_SE': round(tr_shannon.sem(), 4),
    })

    result_df = pd.DataFrame(records)
    print('\n  Summary:')
    print(result_df.to_string(index=False))

    # ── Kruskal-Wallis: are Space quarters consistent? ───────────────────────
    print('\n[5] Kruskal-Wallis test across quarters (Space only)...')
    kw_results = {}
    for metric, series in [('Nitrification', nitri_abund),
                            ('NC_Completeness', ncycle_comp),
                            ('Shannon', shannon_h)]:
        groups_vals = []
        for q in quarters:
            q_ids = sp_meta.index[sp_meta['Quarter'] == q]
            q_ids = q_ids.intersection(series.index)
            if len(q_ids) >= 2:
                groups_vals.append(series[q_ids].values)
        if len(groups_vals) >= 2:
            try:
                stat, p = kruskal(*groups_vals)
                kw_results[metric] = (round(stat, 3), round(p, 6))
                print(f'  {metric}: H={stat:.3f}, p={p:.4f}')
            except Exception as e:
                print(f'  {metric}: KW failed ({e})')
                kw_results[metric] = (None, None)
        else:
            print(f'  {metric}: insufficient quarters for KW')
            kw_results[metric] = (None, None)

    # ── Save results ─────────────────────────────────────────────────────────
    print('\n[6] Saving results...')
    result_df.to_csv(OUTPUT_DIR / 'TemporalStability_Results.csv', index=False)
    print('  Saved: TemporalStability_Results.csv')

    # ── Plot ─────────────────────────────────────────────────────────────────
    print('\n[7] Plotting...')

    sp_df = result_df[result_df['Group'] == 'Space_Flight']
    tr_row = result_df[result_df['Group'] == 'Terrestrial_Soil'].iloc[0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    metrics = [
        ('Nitri_Mean', 'Nitri_SE', 'Nitrification Abundance (%)',
         'Nitrification\nAcross Quarters', 'Nitrification'),
        ('NC_Mean',    'NC_SE',    'N-Cycle Completeness (0–1)',
         'N-Cycle Completeness\nAcross Quarters', 'NC_Completeness'),
        ('Shannon_Mean', 'Shannon_SE', "Shannon H'",
         "Shannon Diversity\nAcross Quarters", 'Shannon'),
    ]

    for ax, (mean_col, se_col, ylabel, title, kw_key) in zip(axes, metrics):
        if len(sp_df) == 0:
            ax.text(0.5, 0.5, 'No quarterly data', transform=ax.transAxes,
                    ha='center', va='center')
            continue

        x = np.arange(len(sp_df))
        ax.errorbar(x, sp_df[mean_col], yerr=sp_df[se_col],
                    color=COLORS['Space_Flight'], marker='o', ms=8, lw=2,
                    capsize=4, label='Space Flight (per quarter)')

        # Terrestrial reference band
        tr_m = tr_row[mean_col]
        tr_se = tr_row[se_col]
        ax.axhline(tr_m, color=COLORS['Terrestrial_Soil'], lw=2,
                   linestyle='--', label=f'Terrestrial Soil (mean)')
        ax.fill_between([-0.5, len(sp_df) - 0.5],
                        tr_m - tr_se, tr_m + tr_se,
                        color=COLORS['Terrestrial_Soil'], alpha=0.15)

        ax.set_xticks(x)
        ax.set_xticklabels(sp_df['Quarter'].values, fontsize=11)
        ax.set_xlabel('Quarter', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.5, len(sp_df) - 0.5)

        # Annotate KW result
        kw_h, kw_p = kw_results.get(kw_key, (None, None))
        if kw_h is not None:
            ax.text(0.98, 0.02,
                    f'KW across quarters:\nH={kw_h:.2f}, p={kw_p:.4f}',
                    transform=ax.transAxes, ha='right', va='bottom',
                    fontsize=8, color='gray',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                              alpha=0.7, edgecolor='gray'))

    fig.suptitle(
        'Temporal Stability: Space Flight Functional Metrics Across Q1–Q4\n'
        '(Consistency indicates no temporal batch effect)',
        fontsize=12, fontweight='bold'
    )

    for ext in ['png', 'pdf']:
        out = OUTPUT_DIR / f'Main_Fig_TemporalStability.{ext}'
        fig.savefig(out, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f'  Saved: {out.name}')
    plt.close()

    print('\nDone.')
