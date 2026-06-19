# Biolab 🧬

A Python-based bioinformatics analysis toolkit built from scratch.
Analyzes DNA sequences and protein properties of cancer-associated genes.

## Project Overview

This project was built progressively from zero Python knowledge to a 
working multi-gene analysis pipeline. It uses BioPython to fetch, parse, 
and analyze real human gene sequences from NCBI.

**Genes analyzed:** TP53, BRCA1, KRAS (cancer-associated genes)

## Tools

### 1. sequence_analyzer.py
Command-line FASTA sequence analysis tool.

**Usage:**
python3 scripts/sequence_analyzer.py <input.fasta> <output.txt>

**Features:**
- Parses single and multi-line FASTA files
- Calculates GC and AT content
- Identifies highest and lowest GC sequences
- Generates formatted analysis report

---

### 2. multi_gene_analysis.py
Full pipeline — analyzes multiple genes in one run.

**Features:**
- Fetches protein sequences directly from NCBI
- Calculates nucleotide properties (GC%, AT%, length)
- Calculates protein properties (molecular weight, pI, stability)
- Generates comprehensive comparison report

---

### 3. compare_genes.py
Side-by-side GC content comparison across multiple genes.

---

### 4. translate_genes.py
Translates mRNA sequences to protein using GenBank CDS annotations.

### 5. gene_pipeline.py
Complete automated pipeline — fetch, analyze, and export.

**Usage:**
```bash
python3 scripts/gene_pipeline.py TP53 BRCA1 KRAS EGFR
```

**Features:**
- Fetches any genes automatically from NCBI
- Calculates nucleotide and protein properties
- Loads results into pandas DataFrame
- Prints comparison table and statistics
- Exports timestamped CSV and summary report

## Sample Results (multi_gene_analysis.py)

| Gene  | mRNA  | GC%   | Protein | MW (kDa) | pI   | Stability |
|-------|-------|-------|---------|----------|------|-----------|
| TP53  | 2512  | 53.38 | 393 aa  | 43.7     | 6.33 | Unstable  |
| BRCA1 | 7088  | 41.77 | 1863 aa | 207.7    | 5.29 | Unstable  |
| KRAS  | 5306  | 36.28 | 188 aa  | 21.4     | 8.24 | Unstable  |
| EGFR  | 9905  | 47.78 | 1210 aa | 134.3    | 6.26 | Unstable  |

## Requirements
Python 3.x

BioPython

Install BioPython:
pip3 install biopython --break-system-packages

## Skills Demonstrated

- FASTA file parsing
- GenBank format handling
- NCBI Entrez API usage
- Protein sequence analysis
- Multi-gene comparative genomics
- Command-line tool design
- Git version control

## Background

Built as part of a self-directed bioinformatics learning path
targeting computational drug discovery (CADD) and bioinformatics
roles. Genes selected for their clinical relevance in oncology.
