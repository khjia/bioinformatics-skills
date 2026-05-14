#!/usr/bin/env python3
"""
WGCNA module per-pathway TF-Structural gene co-expression network plot.
Inner ring: TF genes (pathway TFs + top-N connected module TFs)
Outer ring: Structural genes (non-TF pathway genes)
Edge width: proportional to WGCNA co-expression weight
TF->Structural edges are bolded.

Usage: Edit CONFIG section, then run: python3 plot_tf_network.py
"""
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import defaultdict

# ---- CONFIG (edit these) ----
EDGE_FILE = '/path/to/module_CytoscapeInput.txt'  # WGCNA edge list (fromNode, toNode, weight, ...)
GO_FILE = '/path/to/module_GO_simplified.txt'      # clusterProfiler GO result
TF_FILE = '/path/to/module_itak_tf.tsv'            # iTAK classification output
OUT_DIR = '/path/to/output/figures'

# NPG color palette
C_TF = '#E64B35'
C_STRUCT = '#4DBBD5'
C_EDGE_TF = '#E64B35'
C_EDGE_SS = '#91D1C2'
C_BG = '#FFFFFF'

# Pathways to plot: {output_suffix: GO_Description_exact_match}
PATHWAYS = {
    'water_deprivation': 'response to water deprivation',
    'ros': 'response to reactive oxygen species',
    'hyperosmotic': 'hyperosmotic salinity response',
}

N_TOP_TF = 5  # Extra module TFs to add per pathway
# ---- END CONFIG ----


def main():
    print("[1] Loading data...")
    edges_df = pd.read_csv(EDGE_FILE, sep='\t')
    go_df = pd.read_csv(GO_FILE, sep='\t')
    tf_df = pd.read_csv(TF_FILE, sep='\t')
    tf_set = set(tf_df['gene_id'])
    tf_family = dict(zip(tf_df['gene_id'], tf_df['TF_family']))
    print(f"  Edges: {len(edges_df)}, TFs: {len(tf_set)}")

    for pathway_key, pathway_desc in PATHWAYS.items():
        print(f"\n[{pathway_key}] Processing...")

        row = go_df[go_df['Description'] == pathway_desc]
        if len(row) == 0:
            print(f"  WARNING: pathway '{pathway_desc}' not found!")
            continue
        pathway_genes = row.iloc[0]['geneID'].split('/')

        pathway_tfs = [g for g in pathway_genes if g in tf_set]
        pathway_struct = [g for g in pathway_genes if g not in tf_set]

        # Find top connected TFs from module
        struct_set = set(pathway_struct)
        tf_weights = defaultdict(float)
        for _, r in edges_df.iterrows():
            n1, n2, w = r['fromNode'], r['toNode'], r['weight']
            if n1 in tf_set and n2 in struct_set:
                tf_weights[n1] += w
            elif n2 in tf_set and n1 in struct_set:
                tf_weights[n2] += w

        for t in pathway_tfs:
            tf_weights.pop(t, None)

        top_tfs = sorted(tf_weights.items(), key=lambda x: -x[1])[:N_TOP_TF]
        extra_tfs = [g for g, w in top_tfs]

        all_tfs = pathway_tfs + extra_tfs
        all_struct = pathway_struct
        all_nodes = set(all_tfs + all_struct)

        print(f"  Pathway TFs: {len(pathway_tfs)}, Extra TFs: {len(extra_tfs)}, Structural: {len(all_struct)}")

        # Build subgraph
        sub_edges = edges_df[
            (edges_df['fromNode'].isin(all_nodes)) & (edges_df['toNode'].isin(all_nodes))
        ]
        G = nx.Graph()
        for node in all_nodes:
            G.add_node(node)
        for _, r in sub_edges.iterrows():
            G.add_edge(r['fromNode'], r['toNode'], weight=r['weight'])

        print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

        # Ring layout
        pos = {}
        n_tf, n_struct = len(all_tfs), len(all_struct)
        r_inner, r_outer = 1.0, 2.2

        for i, node in enumerate(all_tfs):
            angle = 2 * np.pi * i / max(n_tf, 1) - np.pi / 2
            pos[node] = (r_inner * np.cos(angle), r_inner * np.sin(angle))
        for i, node in enumerate(all_struct):
            angle = 2 * np.pi * i / max(n_struct, 1) - np.pi / 2
            pos[node] = (r_outer * np.cos(angle), r_outer * np.sin(angle))

        # Draw
        fig, ax = plt.subplots(1, 1, figsize=(12, 12), facecolor=C_BG)
        ax.set_facecolor(C_BG)

        edge_list = list(G.edges(data=True))
        weights = [d['weight'] for _, _, d in edge_list]
        if weights:
            w_min, w_max = min(weights), max(weights)
            w_range = w_max - w_min if w_max > w_min else 1
        else:
            w_min, w_max, w_range = 0, 1, 1

        for u, v, d in edge_list:
            w = d['weight']
            base_width = 0.5 + 3.0 * (w - w_min) / w_range

            if (u in tf_set and v in struct_set) or (v in tf_set and u in struct_set):
                width = base_width * 1.8
                color, alpha = C_EDGE_TF, 0.6
            else:
                width = base_width * 0.8
                color, alpha = C_EDGE_SS, 0.3

            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color=color, linewidth=width, alpha=alpha, zorder=1)

        for node in all_tfs:
            ax.scatter(*pos[node], s=800, c=C_TF, edgecolors='#333', linewidths=1.5, zorder=3)
        for node in all_struct:
            ax.scatter(*pos[node], s=500, c=C_STRUCT, edgecolors='#333', linewidths=1.0, zorder=3)

        # Labels
        for node in all_nodes:
            x, y = pos[node]
            label = node.replace('LOC_Os', '')
            angle = np.arctan2(y, x)

            if node in tf_set:
                fam = tf_family.get(node, '')
                label = f"{label}\n({fam})"
                dist, fontsize, fw = 0.25, 7.5, 'bold'
            else:
                dist, fontsize, fw = 0.22, 7, 'normal'

            lx = x + dist * np.cos(angle)
            ly = y + dist * np.sin(angle)
            ha = 'center' if abs(np.cos(angle)) < 0.3 else ('left' if np.cos(angle) > 0 else 'right')

            ax.annotate(label, (lx, ly), fontsize=fontsize, ha=ha, va='center',
                        color='black', fontweight=fw, zorder=4)

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=C_TF, markersize=15,
                   markeredgecolor='#333', label='Transcription Factor (TF)'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=C_STRUCT, markersize=12,
                   markeredgecolor='#333', label='Structural Gene'),
            Line2D([0], [0], color=C_EDGE_TF, linewidth=3, alpha=0.6,
                   label='TF → Structural (co-expression)'),
            Line2D([0], [0], color=C_EDGE_SS, linewidth=1.5, alpha=0.4,
                   label='Structural ↔ Structural'),
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=10, framealpha=0.9)
        ax.set_title(f'{pathway_desc}\nTF-Structural Gene Co-expression Network',
                     fontsize=14, fontweight='bold', pad=20)
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(-3.2, 3.2)
        ax.set_aspect('equal')
        ax.axis('off')
        plt.tight_layout()

        for fmt in ['png', 'pdf']:
            plt.savefig(f"{OUT_DIR}/pathway_{pathway_key}.{fmt}",
                        dpi=300, bbox_inches='tight', facecolor=C_BG)
        plt.close()
        print(f"  Saved: {OUT_DIR}/pathway_{pathway_key}.png/pdf")

    print("\nDone!")


if __name__ == '__main__':
    main()
