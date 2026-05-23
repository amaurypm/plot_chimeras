# plot_chimeras
Generate combined chimeric protein sequence maps from an MSA.

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

