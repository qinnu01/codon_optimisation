"""
constraints.py
mRNA Secondary Structure & Codon Optimization Constraints
Implements One-Loop, Two-Loop, Multiloop, and Sequence rules.
"""

import math
from typing import Dict, List, Tuple, Set, Optional

STANDARD_CODON_TABLE: Dict[str, List[str]] = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAU', 'AAC'],
    'D': ['GAU', 'GAC'],
    'C': ['UGU', 'UGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'K': ['AAA', 'AAG'],
    'M': ['AUG'],
    'F': ['UUU', 'UUC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    '*': ['UAA', 'UGA', 'UAG']
}

VALID_BASE_PAIRS: Set[Tuple[str, str]] = {
    ('A', 'U'), ('U', 'A'),
    ('G', 'C'), ('C', 'G'),
    ('G', 'U'), ('U', 'G')
}


class mRNAOptimizationConstraints:
    def __init__(
        self,
        protein_sequence: str,
        codon_frequencies: Dict[str, float],
        forbidden_motifs: Optional[List[str]] = None,
        min_hairpin_loop: int = 3,
        max_two_loop_size: int = 30,
        ml_init: float = 3.4,
        ml_u: float = 0.4,
        ml_p: float = 0.1
    ):
        self.protein_sequence = protein_sequence.upper()
        self.codon_freqs = codon_frequencies
        self.forbidden_motifs = forbidden_motifs or []
        self.min_hairpin_loop = min_hairpin_loop
        self.max_two_loop_size = max_two_loop_size

        # Multiloop thermodynamic energy parameters (kcal/mol)
        self.ml_init = ml_init
        self.ml_u = ml_u
        self.ml_p = ml_p

        self.num_amino_acids = len(self.protein_sequence)
        self.mrna_length = self.num_amino_acids * 3
        self.codon_weights = self._precompute_codon_weights()

    def _precompute_codon_weights(self) -> Dict[str, float]:
        weights = {}
        for aa, codons in STANDARD_CODON_TABLE.items():
            max_f = max(self.codon_freqs.get(c, 1e-6) for c in codons)
            for c in codons:
                f_c = self.codon_freqs.get(c, 1e-6)
                weights[c] = -math.log(f_c / max_f)
        return weights

    def validate_codon_selection(self, selected_codons: List[str]) -> bool:
        if len(selected_codons) != self.num_amino_acids:
            return False
        for aa, codon in zip(self.protein_sequence, selected_codons):
            if codon not in STANDARD_CODON_TABLE.get(aa, []):
                return False
        return True

    def validate_forbidden_motifs(self, rna_sequence: str) -> bool:
        for motif in self.forbidden_motifs:
            if motif in rna_sequence:
                return False
        return True

    def is_valid_base_pair(self, b1: str, b2: str) -> bool:
        return (b1.upper(), b2.upper()) in VALID_BASE_PAIRS

    def validate_one_loop(self, i: int, j: int) -> bool:
        """One-loop (Hairpin): Requires at least min_hairpin_loop unpaired bases."""
        return (j - i - 1) >= self.min_hairpin_loop

    def validate_two_loop(self, i: int, j: int, k: int, l: int) -> bool:
        """Two-loop (Stack / Bulge / Internal): Bounds total unpaired nucleotides."""
        if not (i < k < l < j):
            return False
        unpaired = (k - i - 1) + (j - l - 1)
        return unpaired <= self.max_two_loop_size

    def calculate_multiloop_penalty(self, num_unpaired: int, num_branches: int) -> float:
        """Multiloop contribution: ML_init + ML_u * u + ML_p * p."""
        return self.ml_init + (self.ml_u * num_unpaired) + (self.ml_p * num_branches)

    def calculate_cai_penalty(self, selected_codons: List[str], lam: float = 1.0) -> float:
        total_weight = sum(self.codon_weights.get(c, 0.0) for c in selected_codons)
        return lam * total_weight

    def evaluate_caimfe(self, selected_codons: List[str], mfe_energy: float, lam: float = 1.0) -> float:
        return mfe_energy + self.calculate_cai_penalty(selected_codons, lam=lam)