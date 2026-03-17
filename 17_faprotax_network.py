"""
FAPROTAX × Network Keystone Cross-Reference
============================================
Links network hub genera (Fig 1E) with FAPROTAX functional assignments.

Key question: Do the genera that collapse as network hubs also carry
the nitrogen-cycling functions that collapsed?

Input:
  - version-2_integrated/Keystone_Genus_Comparison.csv
  - version-2_integrated/FAPROTAX_report.txt
  - version-2_integrated/FAPROTAX_group_comparison.csv

Output (version-2_integrated/):
  - FAPROTAX_Network_CrossRef.csv   : genus × function membership table
  - Main_Fig_FAPROTAX_Network.png/.pdf : bubble plot + heatmap
"""

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
import config

OUTPUT_DIR = config.VERSION_DIR

# N cycle functions to track
N_FUNCTIONS = [
    'nitrification', 'aerobic_ammonia_oxidation', 'aerobic_nitrite_oxidation',
    'denitrification', 'nitrogen_fixation', 'nitrate_reduction',
    'cellulolysis', 'chitinolysis',
]

FUNCTION_LABELS = {
    'nitrification':             'Nitrification',
    'aerobic_ammonia_oxidation': 'Aerobic NH₃ Oxidation',
    'aerobic_nitrite_oxidation': 'Aerobic NO₂ Oxidation',
    'denitrification':           'Denitrification',
    'nitrogen_fixation':         'N₂ Fixation',
    'nitrate_reduction':         'Nitrate Reduction',
    'cellulolysis':              'Cellulolysis',
    'chitinolysis':              'Chitinolysis',
}


# ─────────────────────────────────────────────
# 1. Parse FAPROTAX report → genus → functions
# ─────────────────────────────────────────────

def parse_faprotax_report(report_path):
    """Return dict: function → set of genus names."""
    func2genera = {}
    current_func = None

    with open(report_path) as f:
        for line in f:
            line = line.rstrip()
            # Section header: "# nitrification (17 records):"
            m = re.match(r'^# (\w+) \(\d+ records\):', line)
            if m:
                current_func = m.group(1)
                func2genera[current_func] = set()
                continue
            # Taxon line: "    d__Bacteria; ...; g__Nitrosospira"
            if current_func and line.startswith('    '):
                genus = None
                for part in line.split(';'):
                    part = part.strip()
                    if part.startswith('g__'):
                        genus = part[3:].strip()
                if genus and genus not in ('', 'uncultured', 'metagenome'):
                    func2genera[current_func].add(genus)

    return func2genera


# ─────────────────────────────────────────────
# 2. Cross-reference with keystone genera
# ─────────────────────────────────────────────

def build_crossref(keystone_df, func2genera, log2fc_df):
    """Build genus × function binary membership matrix."""
    genera = keystone_df['Genus'].tolist()

    rows = []
    for genus in genera:
        row = {'Genus': genus,
               'Space_Degree':      keystone_df.loc[keystone_df['Genus'] == genus, 'Space_Degree'].values[0],
               'Terrestrial_Degree': keystone_df.loc[keystone_df['Genus'] == genus, 'Terrestrial_Degree'].values[0]}
        for func in N_FUNCTIONS:
            genera_in_func = func2genera.get(func, set())
            # Match: genus name must appear in the set
            row[func] = int(any(genus.lower() in g.lower() or g.lower() in genus.lower()
                                for g in genera_in_func))
        rows.append(row)

    crossref = pd.DataFrame(rows)

    # Also collect which genera ARE nitrifiers (but not in keystone list)
    nitrif_genera = func2genera.get('nitrification', set()) | \
                    func2genera.get('aerobic_ammonia_oxidation', set())
    nitrif_not_keystone = sorted(nitrif_genera - set(genera))

    return crossref, nitrif_not_keystone


# ─────────────────────────────────────────────
# 3. Plot
# ─────────────────────────────────────────────

def plot_crossref(crossref, nitrif_not_keystone, log2fc_df):
    fig = plt.figure(figsize=(16, 8))
    gs  = fig.add_gridspec(1, 2, width_ratios=[1.4, 1], wspace=0.35)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # ── Panel A: bubble plot (Space degree vs Terrestrial degree) ──
    # Color by functional role
    def genus_color(row):
        if row['nitrification'] or row['aerobic_ammonia_oxidation'] or row['aerobic_nitrite_oxidation']:
            return '#27ae60', 'Nitrifier'
        if row['denitrification']:
            return '#e74c3c', 'Denitrifier'
        if row['nitrogen_fixation']:
            return '#8e44ad', 'N₂ Fixer'
        return '#95a5a6', 'No N function'

    colors, roles = zip(*[genus_color(r) for _, r in crossref.iterrows()])

    sc = ax1.scatter(crossref['Terrestrial_Degree'], crossref['Space_Degree'],
                     s=180, c=colors, edgecolors='white', linewidths=0.8, zorder=3)

    for _, row in crossref.iterrows():
        ax1.annotate(row['Genus'],
                     (row['Terrestrial_Degree'], row['Space_Degree']),
                     textcoords='offset points', xytext=(5, 4),
                     fontsize=8, fontstyle='italic')

    ax1.plot([0, 9], [0, 9], '--', color='#bdc3c7', linewidth=1, zorder=1)
    ax1.set_xlabel('Terrestrial Network Degree', fontsize=11)
    ax1.set_ylabel('Space Flight Network Degree', fontsize=11)
    ax1.set_title('A  Network Connectivity vs. N Function', fontsize=11, fontweight='bold', loc='left')
    ax1.set_xlim(-0.5, 9.5); ax1.set_ylim(-0.5, 8)

    legend_elements = [
        mpatches.Patch(facecolor='#27ae60', label='Nitrifier'),
        mpatches.Patch(facecolor='#e74c3c', label='Denitrifier'),
        mpatches.Patch(facecolor='#8e44ad', label='N₂ Fixer'),
        mpatches.Patch(facecolor='#95a5a6', label='No N function'),
    ]
    ax1.legend(handles=legend_elements, fontsize=8, loc='upper right')
    ax1.text(0.02, 0.98, 'Above diagonal: space-enriched connectivity\nBelow diagonal: terrestrial-enriched',
             transform=ax1.transAxes, fontsize=7.5, va='top', color='#666')

    # ── Panel B: log2FC of N functions ──
    funcs_to_plot = [f for f in N_FUNCTIONS if f in log2fc_df.index]
    lfc = log2fc_df.loc[funcs_to_plot, 'Log2FC'].sort_values()
    labels = [FUNCTION_LABELS.get(f, f) for f in lfc.index]
    bar_colors = ['#c0392b' if v < 0 else '#2980b9' for v in lfc.values]

    ax2.barh(range(len(lfc)), lfc.values, color=bar_colors, edgecolor='white', height=0.65)
    ax2.set_yticks(range(len(lfc)))
    ax2.set_yticklabels(labels, fontsize=9)
    ax2.axvline(0, color='black', linewidth=0.8)
    ax2.set_xlabel('Log₂FC (Space / Terrestrial)', fontsize=10)
    ax2.set_title('B  Functional Guild Log₂FC', fontsize=11, fontweight='bold', loc='left')

    for i, v in enumerate(lfc.values):
        ax2.text(v + (0.2 if v >= 0 else -0.2), i, f'{v:+.1f}',
                 va='center', ha='left' if v >= 0 else 'right', fontsize=8)

    # Annotate nitrifier non-keystones
    note = f"Nitrifier genera absent from network:\n" + \
           '\n'.join(f"  • {g}" for g in nitrif_not_keystone[:8])
    ax2.text(1.02, 0.5, note, transform=ax2.transAxes,
             fontsize=7.5, va='center', color='#555',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f9fa', edgecolor='#dee2e6'))

    fig.suptitle('Network Hub Collapse Coincides with Nitrogen Cycle Functional Collapse\n'
                 'FAPROTAX × Co-occurrence Network Cross-reference',
                 fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()

    png = OUTPUT_DIR / 'Main_Fig_FAPROTAX_Network.png'
    pdf = OUTPUT_DIR / 'Main_Fig_FAPROTAX_Network.pdf'
    plt.savefig(png, dpi=300, bbox_inches='tight')
    plt.savefig(pdf, bbox_inches='tight')
    plt.close()
    print(f'  Saved: {png.name}')
    print(f'  Saved: {pdf.name}')


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('FAPROTAX × Network Keystone Cross-Reference')
    print('=' * 60)

    keystone = pd.read_csv(OUTPUT_DIR / 'Keystone_Genus_Comparison.csv')
    log2fc   = pd.read_csv(OUTPUT_DIR / 'FAPROTAX_group_comparison.csv', index_col=0)
    report   = OUTPUT_DIR / 'FAPROTAX_report.txt'

    print('\n[1] Parsing FAPROTAX report...')
    func2genera = parse_faprotax_report(report)
    for f in N_FUNCTIONS:
        print(f'  {f}: {len(func2genera.get(f, set()))} genera  → {sorted(func2genera.get(f, set()))[:5]}')

    print('\n[2] Cross-referencing with keystone genera...')
    crossref, nitrif_not_keystone = build_crossref(keystone, func2genera, log2fc)

    out_csv = OUTPUT_DIR / 'FAPROTAX_Network_CrossRef.csv'
    crossref.to_csv(out_csv, index=False)
    print(f'  Saved: {out_csv.name}')
    print('\n  Nitrifier genera NOT in keystone list:')
    for g in nitrif_not_keystone:
        print(f'    • {g}')

    print('\n[3] Plotting...')
    plot_crossref(crossref, nitrif_not_keystone, log2fc)

    print('\nDone.')
