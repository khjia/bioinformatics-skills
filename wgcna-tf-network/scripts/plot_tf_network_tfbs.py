#!/usr/bin/env python3
"""
Per-pathway TF-structural gene co-expression network plot with TFBS evidence.

Inner ring  : TF genes (red, NPG #E64B35)
Outer ring  : Structural genes (cyan, NPG #4DBBD5)
Edges:
  - TF→Structural WITH TFBS evidence : bold red solid (#E64B35)
  - TF→Structural WITHOUT TFBS       : thin gray dashed (#999999)
  - Structural ↔ Structural          : thin cyan solid (#91D1C2)

Inputs:
  --edges       WGCNA Cytoscape edge file (fromNode toNode weight ...)
  --go          GO enrichment table with columns Description, geneID (slash-sep)
  --tf          TF annotation TSV: gene_id, TF_family
  --tfbs        JSON: gene_id -> [TF_family,...]  (from build_tfbs_evidence.py)
  --pathways    JSON dict: {file_suffix: "GO Description"} — pathways to plot
  --out-dir     output directory
  --top-n       extra top-connected module TFs to add per pathway (default 5)
"""
import argparse
import json
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# NPG palette
C_TF = '#E64B35'
C_STRUCT = '#4DBBD5'
C_EDGE_TFBS = '#E64B35'
C_EDGE_NO_TFBS = '#999999'
C_EDGE_SS = '#91D1C2'
C_BG = '#FFFFFF'

def plot_pathway(pathway_key, pathway_desc, edges_df, go_df, tf_set, tf_family,
                 tfbs_evidence, out_dir, top_n):
    row = go_df[go_df['Description'] == pathway_desc]
    if len(row) == 0:
        print(f"[{pathway_key}] WARNING: pathway not found")
        return
    pathway_genes = row.iloc[0]['geneID'].split('/')
    pathway_tfs = [g for g in pathway_genes if g in tf_set]
    pathway_struct = [g for g in pathway_genes if g not in tf_set]
    struct_set = set(pathway_struct)

    # Top connected module TFs (not already in pathway)
    tf_w = defaultdict(float)
    for _, r in edges_df.iterrows():
        n1, n2, w = r['fromNode'], r['toNode'], r['weight']
        if n1 in tf_set and n2 in struct_set:
            tf_w[n1] += w
        elif n2 in tf_set and n1 in struct_set:
            tf_w[n2] += w
    for t in pathway_tfs:
        tf_w.pop(t, None)
    extra_tfs = [g for g, _ in sorted(tf_w.items(), key=lambda x: -x[1])[:top_n]]

    all_tfs = pathway_tfs + extra_tfs
    all_struct = pathway_struct
    all_nodes = set(all_tfs + all_struct)

    sub = edges_df[edges_df['fromNode'].isin(all_nodes) & edges_df['toNode'].isin(all_nodes)]
    G = nx.Graph()
    G.add_nodes_from(all_nodes)
    for _, r in sub.iterrows():
        G.add_edge(r['fromNode'], r['toNode'], weight=r['weight'])

    # Classify edges
    tfbs_e, no_tfbs_e = set(), set()
    for u, v in G.edges():
        if u in tf_set and v in struct_set:
            tf_g, st_g = u, v
        elif v in tf_set and u in struct_set:
            tf_g, st_g = v, u
        else:
            continue
        fam = tf_family.get(tf_g, '')
        if fam in tfbs_evidence.get(st_g, []):
            tfbs_e.add((u, v))
        else:
            no_tfbs_e.add((u, v))

    print(f"[{pathway_key}] nodes={G.number_of_nodes()} edges={G.number_of_edges()} "
          f"TFBS={len(tfbs_e)} no-TFBS={len(no_tfbs_e)}")

    # Layout
    pos = {}
    r_in, r_out = 1.0, 2.2
    for i, n in enumerate(all_tfs):
        a = 2 * np.pi * i / max(len(all_tfs), 1) - np.pi / 2
        pos[n] = (r_in * np.cos(a), r_in * np.sin(a))
    for i, n in enumerate(all_struct):
        a = 2 * np.pi * i / max(len(all_struct), 1) - np.pi / 2
        pos[n] = (r_out * np.cos(a), r_out * np.sin(a))

    fig, ax = plt.subplots(figsize=(12, 12), facecolor=C_BG)
    ax.set_facecolor(C_BG)

    weights = [d['weight'] for _, _, d in G.edges(data=True)]
    if weights:
        w_min, w_max = min(weights), max(weights)
        w_rng = w_max - w_min if w_max > w_min else 1
    else:
        w_min, w_rng = 0, 1

    for u, v, d in G.edges(data=True):
        bw = 0.5 + 3.0 * (d['weight'] - w_min) / w_rng if w_rng > 0 else 1.0
        if (u, v) in tfbs_e or (v, u) in tfbs_e:
            width, color, alpha, ls = bw * 2.0, C_EDGE_TFBS, 0.7, '-'
        elif (u, v) in no_tfbs_e or (v, u) in no_tfbs_e:
            width, color, alpha, ls = bw * 0.6, C_EDGE_NO_TFBS, 0.4, '--'
        else:
            width, color, alpha, ls = bw * 0.7, C_EDGE_SS, 0.3, '-'
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=color, linewidth=width, alpha=alpha, linestyle=ls, zorder=1)

    for n in all_tfs:
        ax.scatter(*pos[n], s=800, c=C_TF, edgecolors='#333', linewidths=1.5, zorder=3)
    for n in all_struct:
        ax.scatter(*pos[n], s=500, c=C_STRUCT, edgecolors='#333', linewidths=1.0, zorder=3)

    for n in all_nodes:
        x, y = pos[n]
        a = np.arctan2(y, x)
        if n in tf_set:
            label = f"{n}\n({tf_family.get(n, '')})"
            dist, fs, fw = 0.25, 7.5, 'bold'
        else:
            label, dist, fs, fw = n, 0.22, 7, 'normal'
        lx, ly = x + dist * np.cos(a), y + dist * np.sin(a)
        ha = 'center' if abs(np.cos(a)) < 0.3 else ('left' if np.cos(a) > 0 else 'right')
        ax.annotate(label, (lx, ly), fontsize=fs, ha=ha, va='center',
                    color='black', fontweight=fw, zorder=4)

    leg = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_TF, markersize=15,
               markeredgecolor='#333', label='Transcription Factor (TF)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_STRUCT, markersize=12,
               markeredgecolor='#333', label='Structural Gene'),
        Line2D([0], [0], color=C_EDGE_TFBS, linewidth=3.5, alpha=0.7,
               label='TF → Target (TFBS evidence)'),
        Line2D([0], [0], color=C_EDGE_NO_TFBS, linewidth=1.5, alpha=0.5, linestyle='--',
               label='TF → Target (no TFBS)'),
        Line2D([0], [0], color=C_EDGE_SS, linewidth=1.5, alpha=0.4,
               label='Structural ↔ Structural'),
    ]
    ax.legend(handles=leg, loc='upper left', fontsize=10, framealpha=0.9)
    ax.set_title(f'{pathway_desc}\nTF-Structural Gene Network (TFBS validated)',
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xlim(-3.2, 3.2); ax.set_ylim(-3.2, 3.2)
    ax.set_aspect('equal'); ax.axis('off')
    plt.tight_layout()
    for fmt in ('png', 'pdf'):
        plt.savefig(f"{out_dir}/pathway_{pathway_key}_tfbs.{fmt}",
                    dpi=300, bbox_inches='tight', facecolor=C_BG)
    plt.close()
    print(f"[{pathway_key}] saved pathway_{pathway_key}_tfbs.png/pdf")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--edges', required=True)
    ap.add_argument('--go', required=True)
    ap.add_argument('--tf', required=True)
    ap.add_argument('--tfbs', required=True)
    ap.add_argument('--pathways', required=True, help='JSON dict: {key: GO description}')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--top-n', type=int, default=5)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    edges_df = pd.read_csv(args.edges, sep='\t')
    go_df = pd.read_csv(args.go, sep='\t')
    tf_df = pd.read_csv(args.tf, sep='\t')
    tf_set = set(tf_df['gene_id'])
    tf_family = dict(zip(tf_df['gene_id'], tf_df['TF_family']))
    with open(args.tfbs) as f:
        tfbs_evidence = json.load(f)
    with open(args.pathways) as f:
        pathways = json.load(f)

    print(f"edges={len(edges_df)} TFs={len(tf_set)} TFBS-genes={len(tfbs_evidence)}")
    for k, desc in pathways.items():
        plot_pathway(k, desc, edges_df, go_df, tf_set, tf_family,
                     tfbs_evidence, args.out_dir, args.top_n)

if __name__ == '__main__':
    main()
