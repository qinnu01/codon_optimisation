"""
run_pipeline.py
Master execution file that runs both sequence optimization (run_optimization.py)
and secondary structure visualization (optimize_and_visualize.py) for any input protein sequence.
"""

import sys
import argparse
from constraints import mRNAOptimizationConstraints, STANDARD_CODON_TABLE
from run_optimization import optimize_mrna
from optimize_and_visualize import compute_mfe_and_structure, visualize_secondary_structure, HUMAN_CODON_FREQUENCIES


def execute_full_pipeline(protein_seq: str, lambda_param: float = 1.0, forbidden_motifs: list = None):
    """
    Executes the end-to-end mRNA design pipeline:
      1. Validates protein input sequence.
      2. Performs CAIMFE co-optimization to find the best mRNA coding sequence.
      3. Computes the MFE secondary structure traceback.
      4. Generates and saves Arc Diagram & Network Topology visualization plots.
    """
    protein_seq = protein_seq.upper().strip()
    forbidden_motifs = forbidden_motifs or []

    # Check for valid standard amino acids
    invalid_chars = [aa for aa in protein_seq if aa not in STANDARD_CODON_TABLE]
    if invalid_chars:
        raise ValueError(f"Invalid amino acid characters in input protein sequence: {set(invalid_chars)}")

    print("=" * 60)
    print("      mRNA STRUCTURE & CODON CO-OPTIMIZATION PIPELINE")
    print("=" * 60)
    print(f"Input Target Protein : {protein_seq} (Length: {len(protein_seq)} AAs)")
    print(f"Target mRNA Length   : {len(protein_seq) * 3} nucleotides")
    print(f"Lambda Factor        : {lambda_param}")
    print(f"Forbidden Motifs     : {forbidden_motifs if forbidden_motifs else 'None'}\n")

    # Step 1: Execute Optimization (run_optimization.py logic)
    print("[1/2] Running mRNA Sequence Co-Optimization...")
    best_rna, best_score = optimize_mrna(
        protein_seq=protein_seq,
        lambda_param=lambda_param,
        forbidden_motifs=forbidden_motifs
    )

    # Step 2: Extract Structure & Generate Visualizations (optimize_and_visualize.py logic)
    print("\n[2/2] Extracting Secondary Structure & Generating Visualizations...")
    constraints = mRNAOptimizationConstraints(
        protein_sequence=protein_seq,
        codon_frequencies=HUMAN_CODON_FREQUENCIES,
        forbidden_motifs=forbidden_motifs
    )

    mfe_energy, base_pairs = compute_mfe_and_structure(best_rna, constraints)

    print("\n" + "=" * 60)
    print("                     PIPELINE RESULTS")
    print("=" * 60)
    print(f"Target Protein  : {protein_seq}")
    print(f"Optimized mRNA  : {best_rna}")
    print(f"MFE Free Energy : {mfe_energy:.2f} kcal/mol")
    print(f"CAIMFE Score    : {best_score:.4f}")
    print(f"Base Pair List  : {base_pairs}")
    print("=" * 60)

    # Render Visual Plot
    plot_filename = f"mrna_structure_{protein_seq}.png"
    visualize_secondary_structure(
        sequence=best_rna,
        pairs=base_pairs,
        title=f"Optimized mRNA Secondary Structure (Target: {protein_seq})"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run complete mRNA sequence co-optimization & visualization pipeline.")
    parser.add_argument("--protein", type=str, default="MALFDKSWCY", help="Target protein amino acid sequence (e.g., MALFDKSWCY)")
    parser.add_argument("--lambda_param", type=float, default=1.0, help="Lambda trade-off weight between MFE and CAI")
    parser.add_argument("--forbidden", nargs="*", default=["GGUACC"], help="List of forbidden sequence motifs to exclude")

    args = parser.parse_args()

    execute_full_pipeline(
        protein_seq=args.protein,
        lambda_param=args.lambda_param,
        forbidden_motifs=args.forbidden
    )