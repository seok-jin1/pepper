"""
Network Robustness Simulation
==============================
Simulates targeted (highest-degree-first) and random node removal to quantify
the structural resilience of co-occurrence networks under hub loss.

Key question:
  Is the terrestrial network more vulnerable to targeted removal of keystone
  taxa than the space network (which has already lost them)?

Robustness score R = mean(S_q) for q = 0 → 1
  where S_q = size of largest connected component (LCC) after removing
  fraction q of nodes, normalized to the original network size.

High R  → resilient (LCC remains large even as nodes are removed)
Low  R  → fragile   (LCC collapses quickly)

Fragility index = R_random / R_targeted
  High value → targeted hubs matter; losing them is especially destructive

Prerequisites:
  - version-2_integrated/exported_table_space/feature-table.tsv
  - version-2_integrated/exported_table_terrestrial/feature-table.tsv
  - version-2_integrated/exported_taxonomy/taxonomy.tsv

Run:
  python 20_network_robustness.py

Output (version-2_integrated/):
  - Network_Robustness_Results.csv
  - Network_Robustness_Summary.csv
  - Main_Fig_NetworkRobustness.png/.pdf
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import config

import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Parameters (must match 10_network_analysis.py) ──────────────────────────
R_THRESHOLD   = config.NETWORK_PARAMS["r_threshold"]
P_THRESHOLD   = config.NETWORK_PARAMS["p_threshold"]
TOP_N         = 15
N_RANDOM_ITER = 1000
RANDOM_SEED   = config.RANDOM_SEED
OUTPUT_DIR    = config.VERSION_DIR

COLORS = {'Space_Flight': '#3498db', 'Terrestrial_Soil': '#2ecc71'}
LABELS = {'Space_Flight': 'Space Flight', 'Terrestrial_Soil': 'Terrestrial Soil'}


# ─────────────────────────────────────────────────────────────────────────────
# 1. Network reconstruction (mirrors 10_network_analysis.py)
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
                name = part.split('__')[-1].strip()
                return name if name else 'Unassigned'
        return 'Unassigned'

    tax['Genus'] = tax['Taxon'].apply(get_genus)
    shared = table.index.intersection(tax.index)
    table  = table.loc[shared]
    table['Genus'] = tax.loc[shared, 'Genus']
    genus_table = table.groupby('Genus').sum()
    genus_table = genus_table.drop(index='Unassigned', errors='ignore')
    return genus_table.div(genus_table.sum(axis=0), axis=1) * 100


def select_top_shared_genera(sp_rel, tr_rel, n=TOP_N):
    shared = sp_rel.index.intersection(tr_rel.index)
    combined_mean = (sp_rel.loc[shared].mean(axis=1) +
                     tr_rel.loc[shared].mean(axis=1)) / 2
    return combined_mean.sort_values(ascending=False).head(n).index


def build_network(genus_rel, top_genera):
    df   = genus_rel.loc[top_genera]
    n    = len(top_genera)
    data = df.values.T  # (samples × genera)
    corr, pvals = spearmanr(data)
    if n == 2:
        corr  = np.array([[1, corr], [corr, 1]])
        pvals = np.array([[0, pvals], [pvals, 0]])

    G = nx.Graph()
    G.add_nodes_from(top_genera)
    for i in range(n):
        for j in range(i + 1, n):
            r, p = corr[i, j], pvals[i, j]
            if abs(r) > R_THRESHOLD and p < P_THRESHOLD and not np.isnan(r):
                G.add_edge(top_genera[i], top_genera[j], weight=r)
    return G


# ─────────────────────────────────────────────────────────────────────────────
# 2. Attack simulations
# ─────────────────────────────────────────────────────────────────────────────

def targeted_attack(G):
    """Remove the current highest-degree node at each step."""
    N      = len(G.nodes())
    G_copy = G.copy()
    fracs  = [0.0]
    sizes  = [1.0]

    while G_copy.number_of_nodes() > 0:
        # Tie-break by node name for reproducibility
        node = max(G_copy.degree(), key=lambda x: (x[1], x[0]))[0]
        G_copy.remove_node(node)
        n_left = G_copy.number_of_nodes()
        if n_left > 0:
            lcc_size = len(max(nx.connected_components(G_copy), key=len))
            sizes.append(lcc_size / N)
        else:
            sizes.append(0.0)
        fracs.append(1 - n_left / N)

    return fracs, sizes


def random_attack(G, n_iter=N_RANDOM_ITER, seed=RANDOM_SEED):
    """Average robustness curve over n_iter random removal orders."""
    N     = len(G.nodes())
    nodes = list(G.nodes())
    rng   = np.random.default_rng(seed)

    all_sizes = np.zeros((n_iter, N + 1))
    all_sizes[:, 0] = 1.0

    for it in range(n_iter):
        G_copy = G.copy()
        order  = rng.permutation(nodes)
        for i, node in enumerate(order):
            G_copy.remove_node(node)
            n_left = G_copy.number_of_nodes()
            if n_left > 0:
                lcc_size = len(max(nx.connected_components(G_copy), key=len))
                all_sizes[it, i + 1] = lcc_size / N
            else:
                all_sizes[it, i + 1] = 0.0

    fracs = np.linspace(0, 1, N + 1).tolist()
    return fracs, all_sizes.mean(axis=0).tolist(), all_sizes.std(axis=0).tolist()


def robustness_score(sizes):
    """R = area under robustness curve = mean(S_q)."""
    return float(np.mean(sizes))


# ─────────────────────────────────────────────────────────────────────────────
# 3. Plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_robustness(results_dict):
    """
    results_dict keys: 'Space_Flight', 'Terrestrial_Soil'
    Each value: dict with keys targ_f, targ_s, rand_f, rand_s, rand_std,
                R_targ, R_rand, fragility
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    attack_titles = [
        ('Targeted Attack\n(highest-degree node first)', 'targ_f', 'targ_s', None, 'R_targ'),
        (f'Random Attack\n(mean ± SD, n={N_RANDOM_ITER} iterations)', 'rand_f', 'rand_s', 'rand_std', 'R_rand'),
    ]

    for ax, (title, fkey, skey, sdkey, rkey) in zip(axes, attack_titles):
        for group, res in results_dict.items():
            col   = COLORS[group]
            label = f"{LABELS[group]}  (R = {res[rkey]:.3f})"
            x     = res[fkey]
            y     = res[skey]

            if sdkey:
                ax.plot(x, y, color=col, lw=2.5, label=label)
                y_arr  = np.array(y)
                sd_arr = np.array(res[sdkey])
                ax.fill_between(x, y_arr - sd_arr, y_arr + sd_arr,
                                color=col, alpha=0.15)
            else:
                ax.step(x, y, where='post', color=col, lw=2.5, label=label)

        ax.set_xlabel('Fraction of Nodes Removed', fontsize=12)
        ax.set_ylabel('Relative Size of Largest Connected Component', fontsize=12)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=10, loc='upper right')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        ax.axhline(0.5, color='gray', lw=1, linestyle=':', alpha=0.7)

    fig.suptitle(
        f'Network Robustness: Space Flight vs. Terrestrial Soil\n'
        f'Top {TOP_N} shared genera | Spearman |r| > {R_THRESHOLD}, p < {P_THRESHOLD}',
        fontsize=13
    )

    for ext in ['png', 'pdf']:
        out = OUTPUT_DIR / f'Main_Fig_NetworkRobustness.{ext}'
        fig.savefig(out, dpi=300 if ext == 'png' else None, bbox_inches='tight')
        print(f'  Saved: {out.name}')
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print('=' * 60)
    print('Network Robustness Simulation')
    print('=' * 60)

    print('\n[1] Reconstructing networks...')
    sp_rel   = aggregate_to_genus(config.FEATURE_TABLE_SPACE, config.TAXONOMY_FILE)
    tr_rel   = aggregate_to_genus(config.FEATURE_TABLE_TERR,  config.TAXONOMY_FILE)
    top_gen  = select_top_shared_genera(sp_rel, tr_rel)
    G_space  = build_network(sp_rel, top_gen)
    G_terr   = build_network(tr_rel, top_gen)
    print(f'  Space       : {G_space.number_of_nodes()} nodes, '
          f'{G_space.number_of_edges()} edges')
    print(f'  Terrestrial : {G_terr.number_of_nodes()} nodes, '
          f'{G_terr.number_of_edges()} edges')

    # Node degree table for reporting
    deg_rows = []
    for g in top_gen:
        deg_rows.append({'Genus': g,
                         'Space_Degree': G_space.degree(g),
                         'Terrestrial_Degree': G_terr.degree(g)})
    deg_df = pd.DataFrame(deg_rows).sort_values('Terrestrial_Degree', ascending=False)
    print('\n  Degree comparison (sorted by terrestrial):')
    print(deg_df.to_string(index=False))

    print('\n[2] Targeted attack...')
    sp_tf, sp_ts = targeted_attack(G_space)
    tr_tf, tr_ts = targeted_attack(G_terr)
    R_sp_t = robustness_score(sp_ts)
    R_tr_t = robustness_score(tr_ts)
    print(f'  Space R_targeted       = {R_sp_t:.4f}')
    print(f'  Terrestrial R_targeted = {R_tr_t:.4f}')

    print(f'\n[3] Random attack ({N_RANDOM_ITER} iterations)...')
    sp_rf, sp_rs, sp_rsd = random_attack(G_space)
    tr_rf, tr_rs, tr_rsd = random_attack(G_terr)
    R_sp_r = robustness_score(sp_rs)
    R_tr_r = robustness_score(tr_rs)
    print(f'  Space R_random         = {R_sp_r:.4f}')
    print(f'  Terrestrial R_random   = {R_tr_r:.4f}')

    # Fragility index = R_random / R_targeted  (>1 means targeted attack hurts more)
    frag_sp = R_sp_r / R_sp_t if R_sp_t > 0 else np.nan
    frag_tr = R_tr_r / R_tr_t if R_tr_t > 0 else np.nan
    print(f'\n  Fragility index (R_random/R_targeted):')
    print(f'  Space       : {frag_sp:.4f}')
    print(f'  Terrestrial : {frag_tr:.4f}')
    print('  (>1 = targeted hub loss is more damaging than random loss)')

    print('\n[4] Saving results...')
    rows = []
    for group, tf, ts, rf, rs, rsd in [
        ('Space_Flight',     sp_tf, sp_ts, sp_rf, sp_rs, sp_rsd),
        ('Terrestrial_Soil', tr_tf, tr_ts, tr_rf, tr_rs, tr_rsd),
    ]:
        for frac, size in zip(tf, ts):
            rows.append({'Group': group, 'Attack': 'Targeted',
                         'Frac_Removed': round(frac, 4),
                         'LCC_Size': round(size, 4),
                         'LCC_Std': None})
        for frac, size, std in zip(rf, rs, rsd):
            rows.append({'Group': group, 'Attack': 'Random',
                         'Frac_Removed': round(frac, 4),
                         'LCC_Size': round(size, 4),
                         'LCC_Std': round(std, 4)})

    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / 'Network_Robustness_Results.csv', index=False)
    print('  Saved: Network_Robustness_Results.csv')

    summary = pd.DataFrame([
        {'Group': 'Space_Flight',     'Attack': 'Targeted',
         'R': round(R_sp_t, 4), 'Fragility_Index': round(frag_sp, 4)},
        {'Group': 'Terrestrial_Soil', 'Attack': 'Targeted',
         'R': round(R_tr_t, 4), 'Fragility_Index': round(frag_tr, 4)},
        {'Group': 'Space_Flight',     'Attack': 'Random',
         'R': round(R_sp_r, 4), 'Fragility_Index': round(frag_sp, 4)},
        {'Group': 'Terrestrial_Soil', 'Attack': 'Random',
         'R': round(R_tr_r, 4), 'Fragility_Index': round(frag_tr, 4)},
    ])
    summary.to_csv(OUTPUT_DIR / 'Network_Robustness_Summary.csv', index=False)
    print('  Saved: Network_Robustness_Summary.csv')

    print('\n--- Summary ---')
    print(summary.to_string(index=False))

    print('\n[5] Plotting...')
    results_dict = {
        'Space_Flight': {
            'targ_f': sp_tf, 'targ_s': sp_ts,
            'rand_f': sp_rf, 'rand_s': sp_rs, 'rand_std': sp_rsd,
            'R_targ': R_sp_t, 'R_rand': R_sp_r, 'fragility': frag_sp,
        },
        'Terrestrial_Soil': {
            'targ_f': tr_tf, 'targ_s': tr_ts,
            'rand_f': tr_rf, 'rand_s': tr_rs, 'rand_std': tr_rsd,
            'R_targ': R_tr_t, 'R_rand': R_tr_r, 'fragility': frag_tr,
        },
    }
    plot_robustness(results_dict)

    print('\nDone.')
