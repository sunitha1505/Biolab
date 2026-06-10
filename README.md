# Biolab

Python scripts for biological sequence analysis.

## Tools

### sequence_analyzer.py
A command-line FASTA sequence analysis tool.

**Usage:**
python3 sequence_analyzer.py <input.fasta> <output.txt>

**Features:**
- Parses single and multi-line FASTA files
- Calculates GC and AT content for each sequence
- Identifies highest and lowest GC sequences
- Generates a formatted analysis report

**Tested on:**
- Human TP53 tumor suppressor gene (NM_000546.6)
- Length: 2512 bp | GC Content: 53.38%

## Requirements
Python 3.x
