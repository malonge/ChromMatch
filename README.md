# ChromMatch

A tool for assigning chromosome labels based on a reference genome. This method is intended to be more sensitive than whole-genome alignment.

## Dependencies
- Python 3
- [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
- [minimap2](https://github.com/lh3/minimap2)
- [RagTag](https://github.com/malonge/RagTag)
- numpy
- pysam
- networkx
- gffutils

## Installation
No installation is needed. Just install the above dependencies and run `python3 chrom_match.py`

## Usage
Suppose there is a "target" genome assembly with 12 chromosome-scale sequences and the goal is to assign chromosome names based on a related reference genome, also with 12 chromsomes.

```
python3 chrom_match.py target.fa reference.fa reference.genes.gff3
```

Suppose that the target and reference assemblies have additional sequences, such as unplaced contigs/scaffolds. Then, supply the sequences to be matched with `-t` and `-r`.

In all cases, the number of target and reference sequencs must match in order to be matched

## Pipeline
1. Write reference transcripts to a FASTA file
2. Align these transcripts to the target assembly with minimap2
3. Process these alignments to build a bipartite graph
    - For each gene, only consider its longest (representative) transcript
    - Only consider representative transcripts if they align with mapq > 10 and coverage >= 85%
    - Target sequences make one set of nodes
    - Reference sequences make the other set of nodes
    - Edges connecting these nodes are decremented by one if they share a transcript
4. Compute a minimum weight full matching for the graph
5. Output the matching solution in AGP format
