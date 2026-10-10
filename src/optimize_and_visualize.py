"""
optimize_and_visualize.py
Pipeline to co-optimize mRNA sequence and visualize its secondary structure.
"""

import math
import itertools
import matplotlib.pyplot as plt
import networkx as nx
from typing import List, Tuple, Dict
from constraints import mRNAOptimizationConstraints, STANDARD_CODON_TABLE

# Human Codon Usage Frequencies
HUMAN_CODON_FREQUENCIES = {
    'AUG': 1.00,  # Met
    'GCU': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11,  # Ala
    'CUU': 0.13, 'CUC': 0.20, 'CUA': 0.07, 'CUG': 0.40, 'UUA': 0.07, 'UUG': 0.13,  # Leu
    'UUU': 0.45, 'UUC': 0.55,  # Phe
}


def compute_mfe_and_structure(sequence: str, constraints: mRNAOptimizationConstraints) -> Tuple[float, List[Tuple[int, int]]]:
    """
    Computes MFE energy and recovers paired base indices using Dynamic Programming.
    """
    n = len(sequence)
    dp = [[0.0 for _ in range(n)] for _ in range(n)]
    trace = {}  # Matrix to track backpair trace

    pair_energies = {
        ('G', 'C'): -3.0, ('C', 'G'): -3.0,
        ('A', 'U'): -2.0, ('U', 'A'): -2.0,
        ('G', 'U'): -1.0, ('U', 'G'): -1.0
    }

    min_loop = constraints.min_hairpin_loop

    for length in range(min_loop + 2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1

            # Case 1: Unpaired base at i or j
            dp[i][j] = dp[i + 1][j]
            trace[(i, j)] = ('unpaired_left', i + 1, j)

            if dp[i][j - 1] < dp[i][j]:
                dp[i][j] = dp[i][j - 1]
                trace[(i, j)] = ('unpaired_right', i, j - 1)

            # Case 2: Hairpin Loop Pairing (i, j)
            if constraints.is_valid_base_pair(sequence[i], sequence[j]) and constraints.validate_one_loop(i, j):
                energy = pair_energies.get((sequence[i], sequence[j]), -1.0)
                val = dp[i + 1][j - 1] + energy
                if val < dp[i][j]:
                    dp[i][j] = val
                    trace[(i, j)] = ('pair', i, j)

            # Case 3: Two-Loop (Stacking / Bulges / Internal Loops)
            for k in range(i + 1, j - 1):
                for l in range(k + 1, j):
                    if constraints.validate_two_loop(i, j, k, l):
                        if constraints.is_valid_base_pair(sequence[i], sequence[j]) and constraints.is_valid_base_pair(sequence[k], sequence[l]):
                            e1 = pair_energies.get((sequence[i], sequence[j]), -1.0)
                            e2 = pair_energies.get((sequence[k], sequence[l]), -1.0)
                            val = dp[k][l] + e1 + e2
                            if val < dp[i][j]:
                                dp[i][j] = val
                                trace[(i, j)] = ('two_loop', (i, j), (k, l))

            # Case 4: Bifurcation
            for k in range(i + 1, j):
                val = dp[i][k] + dp[k + 1][j]
                if val < dp[i][j]:
                    dp[i][j] = val
                    trace[(i, j)] = ('bifurcate', k)

    # Traceback paired base indices
    pairs = []
    stack = [(0, n - 1)]
    visited = set()

    while stack:
        i, j = stack.pop()
        if i >= j or (i, j) in visited:
            continue
        visited.add((i, j))

        action = trace.get((i, j))
        if not action:
            continue

        if action[0] == 'unpaired_left':
            stack.append((action[1], action[2]))
        elif action[0] == 'unpaired_right':
            stack.append((action[1], action[2]))
        elif action[0] == 'pair':
            pairs.append((action[1], action[2]))
            stack.append((action[1] + 1, action[2] - 1))
        elif action[0] == 'two_loop':
            pairs.append(action[1])
            pairs.append(action[2])
            stack.append((action[2][0] + 1, action[2][1] - 1))
        elif action[0] == 'bifurcate':
            k = action[1]
            stack.append((i, k))
            stack.append((k + 1, j))

    return dp[0][n - 1], list(set(pairs))


def visualize_secondary_structure(sequence: str, pairs: List[Tuple[int, int]], title: str):
    """
    Plots both Arc Diagram and Planar Graph visualizations.
    """
    n = len(sequence)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7))

    # --- 1. Arc Diagram ---
    for i, base in enumerate(sequence):
        ax1.plot(i, 0, 'o', color='skyblue', markersize=12)
        ax1.text(i, -0.15, f"{base}\n{i+1}", ha='center', va='top', fontsize=9, fontweight='bold')

    ax1.plot(range(n), [0] * n, color='gray', linestyle='--', zorder=1)

    for i, j in pairs:
        center = (i + j) / 2.0
        radius = (j - i) / 2.0
        arc = plt.Circle((center, 0), radius, color='crimson', fill=False, linewidth=2)
        ax1.add_patch(arc)

    ax1.set_xlim(-1, n)
    ax1.set_ylim(-1, n / 2 + 1)
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax1.set_title(f"{title}\nArc Diagram View", fontsize=11, fontweight='bold')

    # --- 2. Planar Network Graph View ---
    G = nx.Graph()
    for i, base in enumerate(sequence):
        G.add_node(i, label=f"{base}_{i+1}")

    # Backbone Edges
    for i in range(n - 1):
        G.add_edge(i, i + 1, edge_type='backbone')

    # Base Pair Edges
    for i, j in pairs:
        G.add_edge(i, j, edge_type='basepair')

    pos = nx.spring_layout(G, seed=42)
    labels = {i: f"{base}\n{i+1}" for i, base in enumerate(sequence)}

    nx.draw_networkx_nodes(G, pos, ax=ax2, node_color='lightblue', node_size=500)
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax2, font_size=8, font_weight='bold')

    backbone_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'backbone']
    basepair_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'basepair']

    nx.draw_networkx_edges(G, pos, edgelist=backbone_edges, ax=ax2, edge_color='gray', width=1.5, style='dashed')
    nx.draw_networkx_edges(G, pos, edgelist=basepair_edges, ax=ax2, edge_color='crimson', width=2.5)

    ax2.set_title("Planar Network Topology View", fontsize=11, fontweight='bold')
    ax2.axis('off')

    plt.tight_layout()
    plt.savefig("mrna_structure_optimized.png", dpi=300)
    print("Structure plot saved to: mrna_structure_optimized.png")
    plt.show()


if __name__ == "__main__":
    target_protein = "MALFDKSWCY"
    lambda_param = 1.0

    constraints = mRNAOptimizationConstraints(
        protein_sequence=target_protein,
        codon_frequencies=HUMAN_CODON_FREQUENCIES
    )

    # Search synonymous codons
    synonymous_choices = [STANDARD_CODON_TABLE[aa] for aa in target_protein]
    all_codon_combinations = list(itertools.product(*synonymous_choices))

    best_score = float('inf')
    best_rna = None
    best_pairs = []
    best_mfe = 0.0

    print("Running Sequence Co-Optimization Search...")
    for codons in all_codon_combinations:
        codons_list = list(codons)
        rna_seq = "".join(codons_list)

        mfe, pairs = compute_mfe_and_structure(rna_seq, constraints)
        cai_penalty = constraints.calculate_cai_penalty(codons_list, lam=lambda_param)
        score = constraints.evaluate_caimfe(codons_list, mfe, lam=lambda_param)

        if score < best_score:
            best_score = score
            best_rna = rna_seq
            best_pairs = pairs
            best_mfe = mfe

    print("\n=== Optimized Sequence & Secondary Structure ===")
    print(f"Optimal mRNA Sequence : {best_rna}")
    print(f"Predicted MFE Energy   : {best_mfe:.2f} kcal/mol")
    print(f"Predicted Base Pairs   : {best_pairs}")

    # Render Visualizations
    visualize_secondary_structure(
        sequence=best_rna,
        pairs=best_pairs,
        title=f"Optimized mRNA Secondary Structure (Target: {target_protein})"
    )