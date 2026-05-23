# plot_chimeras
Generate combined chimeric protein sequence maps from an MSA.

The segments are inferred by Dynamic Programming (Viterbi algorithm) to find the sequence of parental origins that minimizes the total number of segment switches (crossovers).

## Inputs
* A multiple sequence alignment (MSA) in FASTA format (no '*' for stop codon) containing the target chimeric protein and the parental sequences (required).
* An annotation file (in CSV) format containing general, or target specific annotations (optional).

## Output
A PNG file containing the annotated sequence maps of the targets.

## Installation
As a python script you can just run the plot_chimeras.py file, or put a symbolic link in any directory of your PATH.
The script depends on the following Python libraries:
* argparse
* os
* matplotlib

## Usage
```
plot_chimeras [-h] -i FASTA -t TARGETS [TARGETS ...] -p PARENTALS [PARENTALS ...] [-f FEATURES] [-o OUTPUT] [-v]

options:
  -h, --help            show this help message and exit
  -i, --fasta FASTA     Input MSA FASTA file
  -t, --targets TARGETS [TARGETS ...]
                        List of target sequence IDs (space-separated)
  -p, --parentals PARENTALS [PARENTALS ...]
                        List of parental sequence IDs (space-separated)
  -f, --features FEATURES
                        Optional CSV file with features
  -o, --output OUTPUT   Output PNG file path
  -v, --version         Show program's version number and exit
```

## Example
Check the example folder.
