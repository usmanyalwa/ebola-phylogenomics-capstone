# Ebola Virus 2014 Phylogenomics & Mutation Analysis

## Overview
This repository contains a Python-based computational framework for analyzing genomic variation and substitution metrics across 2014 West African *Zaire ebolavirus* strains, specifically evaluating accession **KJ660348**. 

The pipeline replaces manual alignment inspections with automated sequence processing using Biopython.

## Project Structure
```text
ebola-phylogenomics-capstone/
├── data/
│   └── ebola_sequences.fasta     # Aligned genomic FASTA sequences
├── results/
│   └── analysis_summary.txt      # Generated outputs and logs
├── scripts/
│   └── analyze_ebola.py          # Core Python execution script
├── .gitignore
├── requirements.txt              # Environment dependencies
└── README.md                     # Project documentation

```

## Methodology & Findings

1. **Sequence Alignment:** 10 full-length viral genomes (19,220 bp alignment length).
2. **Terminal Node Analysis vs. Pairwise Global Sums:**
* **Ancestral Terminal Mutation Count:** Direct ancestral reconstruction yields **0** mutations along the specific lineage branch for target accession `KJ660348`.
* **Global Pairwise Substitution Metric:** Pairwise distance sums across all sequence terminal nodes accumulate to **10,644** total site mismatches across the global matrix.



## Getting Started

### Prerequisites

* Python 3.8+
* Biopython & Pandas

### Installation & Execution

1. Clone the repository:
```bash
git clone [https://github.com/usmanyalwa/ebola-phylogenomics-capstone.git](https://github.com/usmanyalwa
/ebola-phylogenomics-capstone.git)
cd ebola-phylogenomics-capstone

```


2. Install dependencies:
```bash
pip install -r requirements.txt

```


3. Run the mutation analysis pipeline:
```bash
python scripts/analyze_ebola.py

```



## Author

* **Usman Yalwa Ibrahim** — M.Sc. Bioinformatics