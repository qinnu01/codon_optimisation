"""
main.py
Execution entrypoint for Quantum mRNA Structure Optimization.
"""

from quantum_mrna_solver import QuantummRNASolver

HUMAN_CODON_FREQUENCIES = {
    'AUG': 1.00,
    'UUU': 0.45, 'UUC': 0.55,
    'CUU': 0.13, 'CUC': 0.20, 'CUA': 0.07, 'CUG': 0.40, 'UUA': 0.07, 'UUG': 0.13
}

def main():
    target_protein = "MFL"
    print(f"=== Quantum Optimization for mRNA Secondary Structure ===")
    print(f"Target Protein: {target_protein}")

    solver = QuantummRNASolver(
        protein_seq=target_protein,
        codon_freqs=HUMAN_CODON_FREQUENCIES,
        lambda_cai=1.5
    )

    solver.build_qubo_model()
    print("Running QAOA Optimization on Quantum Simulator...")
    res = solver.solve_qaoa()

    print("\n=== Quantum Solver Output ===")
    print(f"Optimal QUBO Objective Energy: {res['optimal_value']:.4f}")


if __name__ == "__main__":
    main()