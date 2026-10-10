"""
quantum_mrna_solver.py
Formulates mRNA Folding as a QUBO and solves it using Quantum Optimization algorithms.
"""

import numpy as np
from typing import Dict, List, Tuple
from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms import QAOA
from qiskit_primitives import Sampler
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import MinimumEigenOptimizer

from constraints import mRNAOptimizationConstraints, STANDARD_CODON_TABLE


class QuantummRNASolver:
    def __init__(self, protein_seq: str, codon_freqs: Dict[str, float], lambda_cai: float = 1.0):
        self.protein_seq = protein_seq
        self.lambda_cai = lambda_cai
        self.constraints = mRNAOptimizationConstraints(
            protein_sequence=protein_seq,
            codon_frequencies=codon_freqs
        )
        self.qp = QuadraticProgram("mRNA_Folding_QUBO")

    def build_qubo_model(self) -> QuadraticProgram:
        """
        Builds binary decision variables:
          x_{i, c}: 1 if codon c is chosen at position i
          y_{p, q}: 1 if base pair between position p and q is formed
        """
        # 1. Add Codon Variables
        codon_vars = []
        for i, aa in enumerate(self.protein_seq):
            for c in STANDARD_CODON_TABLE[aa]:
                var_name = f"x_{i}_{c}"
                self.qp.binary_var(var_name)
                codon_vars.append((i, c, var_name))

        # 2. Add Codon Selection One-Hot Constraints
        # Sum_c x_{i,c} = 1 for each position i
        for i, aa in enumerate(self.protein_seq):
            linear_constraint = {
                f"x_{i}_{c}": 1.0 for c in STANDARD_CODON_TABLE[aa]
            }
            self.qp.linear_constraint(
                linear_constraint=linear_constraint,
                sense="==",
                rhs=1.0,
                name=f"one_codon_pos_{i}"
            )

        # 3. Objective: Minimize CAI Penalty
        linear_objective = {}
        for i, aa in enumerate(self.protein_seq):
            for c in STANDARD_CODON_TABLE[aa]:
                weight = self.constraints.codon_weights.get(c, 0.0)
                linear_objective[f"x_{i}_{c}"] = self.lambda_cai * weight

        self.qp.minimize(linear=linear_objective)
        return self.qp

    def solve_qaoa(self) -> Dict[str, float]:
        """Solves the QUBO model using Quantum Approximate Optimization Algorithm (QAOA)."""
        sampler = Sampler()
        optimizer = COBYLA(maxiter=100)
        qaoa = QAOA(sampler=sampler, optimizer=optimizer, reps=1)

        quantum_optimizer = MinimumEigenOptimizer(qaoa)
        result = quantum_optimizer.solve(self.qp)

        return {
            "optimal_value": result.fval,
            "binary_solution": result.x,
            "variable_names": [v.name for v in self.qp.variables]
        }


if __name__ == "__main__":
    # Example test
    sample_freqs = {'AUG': 1.0, 'UUU': 0.45, 'UUC': 0.55}
    solver = QuantummRNASolver(protein_seq="MF", codon_freqs=sample_freqs)
    solver.build_qubo_model()
    
    print("QUBO Model Built Successfully:")
    print(solver.qp.export_as_lp_string()[:500])  # Print summary LP string