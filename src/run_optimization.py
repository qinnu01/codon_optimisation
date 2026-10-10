"""
run_optimization.py
mRNA Folding Optimization Execution Script
Finds optimal mRNA sequence & secondary structure considering CAI and MFE parameters.
"""

import itertools
import math
from constraints import mRNAOptimizationConstraints, STANDARD_CODON_TABLE

HUMAN_CODON_FREQUENCIES = {
    'AUG': 1.00,  # Met
    'GCU': 0.26, 'GCC': 0.40, 'GCA': 0.23, 'GCG': 0.11,  # Ala
    'CUU': 0.13, 'CUC': 0.20, 'CUA': 0.07, 'CUG': 0.40, 'UUA': 0.07, 'UUG': 0.13,  # Leu
    'UUU': 0.45, 'UUC': 0.55,  # Phe
}


def compute_mfe_score(sequence: str, constraints: mRNAOptimizationConstraints) -> float:
    """
    Computes MFE using dynamic programming over One-Loop, Two-Loop, and Multiloop constraints.
    """
    n = len(sequence)
    dp = [[0.0 for _ in range(n)] for _ in range(n)]

    # Pair stacking energies (kcal/mol)
    pair_energies = {
        ('G', 'C'): -3.0, ('C', 'G'): -3.0,
        ('A', 'U'): -2.0, ('U', 'A'): -2.0,
        ('G', 'U'): -1.0, ('U', 'G'): -1.0
    }

    min_loop = constraints.min_hairpin_loop

    for length in range(min_loop + 2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1

            # Option 1: Unpaired bases
            dp[i][j] = min(dp[i + 1][j], dp[i][j - 1])

            # Option 2: Pair base i and base j (One-Loop / Hairpin validation)
            if constraints.is_valid_base_pair(sequence[i], sequence[j]):
                if constraints.validate_one_loop(i, j):
                    energy = pair_energies.get((sequence[i], sequence[j]), -1.0)
                    dp[i][j] = min(dp[i][j], dp[i + 1][j - 1] + energy)

            # Option 3: Two-Loop (Stack / Bulge / Internal Loop checking)
            for k in range(i + 1, j - 1):
                for l in range(k + 1, j):
                    if constraints.validate_two_loop(i, j, k, l):
                        if constraints.is_valid_base_pair(sequence[i], sequence[j]) and constraints.is_valid_base_pair(sequence[k], sequence[l]):
                            e1 = pair_energies.get((sequence[i], sequence[j]), -1.0)
                            e2 = pair_energies.get((sequence[k], sequence[l]), -1.0)
                            dp[i][j] = min(dp[i][j], dp[k][l] + e1 + e2)

            # Option 4: Multiloop Bifurcation
            for k in range(i + 1, j):
                dp[i][j] = min(dp[i][j], dp[i][k] + dp[k + 1][j])

    return dp[0][n - 1]


def optimize_mrna(
    protein_seq: str,
    lambda_param: float = 1.0,
    forbidden_motifs: list = None
):
    constraints = mRNAOptimizationConstraints(
        protein_sequence=protein_seq,
        codon_frequencies=HUMAN_CODON_FREQUENCIES,
        forbidden_motifs=forbidden_motifs
    )

    # Synonymous codon search space
    synonymous_choices = [STANDARD_CODON_TABLE[aa] for aa in protein_seq]
    all_codon_combinations = list(itertools.product(*synonymous_choices))

    best_score = float('inf')
    best_rna = None
    best_codons = None
    best_mfe = 0.0
    best_cai_penalty = 0.0

    print("--- mRNA Structural & Codon Optimization ---")
    print(f"Target Protein : {protein_seq}")
    print(f"Lambda Factor  : {lambda_param}")
    print(f"Combinations   : {len(all_codon_combinations)}\n")

    for codons in all_codon_combinations:
        codons_list = list(codons)
        rna_seq = "".join(codons_list)

        if not constraints.validate_codon_selection(codons_list):
            continue
        if not constraints.validate_forbidden_motifs(rna_seq):
            continue

        mfe = compute_mfe_score(rna_seq, constraints)
        cai_penalty = constraints.calculate_cai_penalty(codons_list, lam=lambda_param)
        caimfe_score = constraints.evaluate_caimfe(codons_list, mfe, lam=lambda_param)

        if caimfe_score < best_score:
            best_score = caimfe_score
            best_rna = rna_seq
            best_codons = codons_list
            best_mfe = mfe
            best_cai_penalty = cai_penalty

    print("=== Best Sequence Found ===")
    print(f"Optimal Codons : {best_codons}")
    print(f"Optimal mRNA   : {best_rna}")
    print(f"Predicted MFE  : {best_mfe:.2f} kcal/mol")
    print(f"CAI Penalty    : {best_cai_penalty:.4f}")
    print(f"CAIMFE Score   : {best_score:.4f}")

    return best_rna, best_score


if __name__ == "__main__":
    target_protein = "MALFDKSWCY"
    forbidden_list = ["GGUACC"]

    optimize_mrna(
        protein_seq=target_protein,
        lambda_param=1.5,
        forbidden_motifs=forbidden_list
    )