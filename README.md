# mRNA Secondary Structure & Codon Co-Optimization Pipeline

An end-to-end computational biology workflow for co-optimizing mRNA coding sequences to balance **thermodynamic secondary structure stability** (Minimum Free Energy / MFE) and **host translation efficiency** (Codon Adaptation Index / CAI), based on the mRNA folding review literature (*Ward et al., 2025*).

This repository contains both standard dynamic programming algorithms (similar to *LinearDesign* and *CDSfold*) and a quantum optimization extension (QUBO / QAOA via Qiskit) for structural folding over codon candidate graphs.

---

## 📌 Problem Formulation & Objective

The primary objective of mRNA therapeutic design is to identify an optimal nucleotide sequence $\pi$ encoding a target amino acid sequence $\alpha$ that minimizes a joint thermodynamic and translation penalty score:

$$\pi^* = \arg\min_{\pi \in \text{CDS}(\alpha)} \left[ \text{MFE}(\pi) - \lambda \cdot \log(\text{CAI}(\pi)) \right]$$

Where:
* **$\text{CDS}(\alpha)$**: The set of all valid mRNA sequences that translate into protein $\alpha$.
* **$\text{MFE}(\pi)$**: Minimum Free Energy of the predicted secondary structure (kcal/mol).
* **$\text{CAI}(\pi)$**: Codon Adaptation Index relative to host organism usage frequency tables.
* **$\lambda$**: User-defined weighting parameter balancing structural stability against codon preference.

---

## 🧬 Incorporated Structural Constraints

The optimization engine models nearest-neighbor thermodynamic secondary structures across specific loop classifications and constraints:

1. **One-Loops (Hairpin Loops):** Closed by a single base pair $(i, j)$ with a mandatory minimum hairpin size constraint of $u \ge 3$ unpaired nucleotides.
2. **Two-Loops (Stacks, Bulges, Internal Loops):** Closed by two base pairs $(i, j)$ and $(k, l)$ where $i < k < l < j$. Includes a upper bound constraint on internal loop size ($u \le 30$ unpaired bases).
3. **Multiloops:** Evaluated using initiation and branching penalties: $e_{\text{ML}} = \text{ML}_{\text{init}} + \text{ML}_u \cdot u + \text{ML}_p \cdot p$.
4. **Base-Pairing Rules:** Enforces canonical base pairs ($\text{A-U}$, $\text{G-C}$, $\text{G-U}$) and orientation symmetry.
5. **Nested Structure Rule:** Ensures secondary structures contain no pseudoknots (crossing arcs are strictly forbidden).
6. **Forbidden Motif Exclusion:** Prevents unwanted sequence patterns (e.g., restriction sites like `GGUACC` or poly-U tracks) across chosen candidate codons.

---

## 📁 Repository Structure

```text
.
├── constraints.py             # Defines structural constraints, loop bounds, and CAI weights
├── run_optimization.py        # Core Dynamic Programming co-optimization algorithm
├── optimize_and_visualize.py  # Secondary structure traceback and visualization plot generator
├── run_pipeline.py            # Master execution script to run the full workflow
├── quantum_mrna_solver.py     # QUBO model builder & QAOA quantum solver (Qiskit extension)
└── README.md                  # Project documentation

```

---

## ⚙️ Installation & Prerequisites

### Prerequisites

* **Python**: Version 3.9 or higher
* **Git**: For cloning the repository

### Installation Steps

1. **Clone the Repository:**
```bash
git clone [https://github.com/your-username/mrna-folding-optimization.git](https://github.com/your-username/mrna-folding-optimization.git)
cd mrna-folding-optimization

```


2. **Create and Activate a Virtual Environment:**
* **Linux/macOS:**
```bash
python3 -m venv venv
source venv/bin/activate

```


* **Windows:**
```cmd
python -m venv venv
venv\Scripts\activate

```




3. **Install Dependencies:**
```bash
pip install matplotlib networkx numpy qiskit qiskit-algorithms qiskit-optimization

```



*(Optional)* For production-grade thermodynamic calculations matching ViennaRNA parity, you can also install the ViennaRNA Python bindings:

```bash
pip install ViennaRNA

```

---

## 🚀 Usage Instructions

### 1. Run the Complete End-to-End Pipeline

The master pipeline (`run_pipeline.py`) takes an amino acid sequence input, computes the optimal synonymous codon sequence, recovers the secondary structure base-pairing traceback, and exports a structural plot.

```bash
python run_pipeline.py --protein MALFDKSWCY --lambda_param 1.0 --forbidden GGUACC

```

#### Command-Line Arguments:

* `--protein`: Target protein sequence in standard single-letter amino acid code (Default: `MALFDKSWCY`).
* `--lambda_param`: Weight parameter $\lambda$ balancing MFE vs. CAI penalty (Default: `1.0`).
* `--forbidden`: Space-separated list of motif strings to exclude from candidate sequences (Default: `GGUACC`).

---

### 2. Run Individual Modules

#### A. Sequence & CAIMFE Optimization Only

Executes the dynamic programming search over the codon graph space and prints the optimal codon selection:

```bash
python run_optimization.py

```

#### B. Secondary Structure Traceback & Plotting

Computes paired indices for the optimal sequence and displays both an **Arc Diagram** and a **Planar Network Graph**:

```bash
python optimize_and_visualize.py

```

#### C. Quantum QUBO / QAOA Optimization

Converts the sequence design choices into a Quadratic Unconstrained Binary Optimization (QUBO) problem and solves it using Qiskit's QAOA solver:

```bash
python quantum_mrna_solver.py

```

---

## 📊 Generated Output Files

When running the pipeline, the following artifacts are generated:

1. **Console Summary:** Displays optimal codons, full nucleotide mRNA sequence, predicted MFE energy (kcal/mol), CAI penalty score, and base pair list.
2. **Visual Plot (`mrna_structure_<PROTEIN>.png`):**
* **Arc Diagram:** Displays backbone nodes and arc connections for base-paired regions.
* **Planar Topology Graph:** Visualizes network graph representation of backbone and pairing interactions.



---

## 📄 References & Literature

* Ward, M., Richardson, M., & Metkar, M. (2025). *"mRNA Folding Algorithms for Structure and Codon Optimization"*. arXiv:2503.19273v1 [q-bio.BM].
* Zhang, H., et al. (2023). *"Algorithm for optimized mRNA design improves stability and immunogenicity"*. Nature, 621, 396–403 (LinearDesign).
* Terai, G., Kamegai, S., & Asai, K. (2016). *"CDSfold: an algorithm for designing a protein-coding sequence with the most stable secondary structure"*. Bioinformatics, 32(6), 828–834.

```

```