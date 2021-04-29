#!/usr/bin/env python

import os
import sys
import time
import argparse

import numpy as np
import networkx as nx
import pysam
import gffutils

from ragtag_utilities.AlignmentReader import PAFReader
from ragtag_utilities.AGPFile import AGPFile
from ragtag_utilities.utilities import run_oe


def log(level, message):
    """ Log messages to standard error. """
    level = level.upper()
    if level not in {"VERSION", "CMD", "INFO", "WARNING", "DEBUG"}:
        raise ValueError("Invalid logging level: {}".format(level))

    sys.stderr.write(time.ctime() + " --- " + level + ": " + message + "\n")
    sys.stderr.flush()


def main():
    parser = argparse.ArgumentParser(description='Make joins and fill gaps in one assembly (target.fa) using sequences from another assembly (query.fa)', usage="ragtag.py patch <target.fa> <query.fa>")

    parser.add_argument("target_asm_fn", metavar="<target.fa>", nargs='?', default="", type=str, help="target FASTA file (uncompressed or bgzipped)")
    parser.add_argument("reference_asm_fn", metavar="<reference.fa>", nargs='?', default="", type=str, help="reference FASTA file (uncompressed or bgzipped)")
    parser.add_argument("reference_genes_fn", metavar="<reference.gff3>", nargs='?', default="", type=str, help="reference genes GFF3 file (uncompressed)")
    parser.add_argument("-t", metavar="PATH", default="", type=str, help="list of target sequences to match [null]")
    parser.add_argument("-r", metavar="PATH", default="", type=str, help="list of reference sequences to match [null]")
    parser.add_argument("-g", metavar="PATH", default="gffread", type=str, help="path to gffread [gffread]")
    parser.add_argument("-m", metavar="PATH", default="minimap2", type=str, help="path to minimap2 [minimap2]")

    args = parser.parse_args()

    # Check for the required arguments
    if not args.target_asm_fn:
        parser.print_help()
        print("\n** The target FASTA file is required **")
        sys.exit(1)

    if not args.reference_asm_fn:
        parser.print_help()
        print("\n** The reference FASTA file is required **")
        sys.exit(1)

    if not args.reference_genes_fn:
        parser.print_help()
        print("\n** The reference GFF3 file is required **")
        sys.exit(1)

    # Parse command line arguments
    target_asm_fn = os.path.abspath(args.target_asm_fn)
    reference_asm_fn = os.path.abspath(args.reference_asm_fn)
    reference_genes_fn = os.path.abspath(args.reference_genes_fn)
    target_seqs_fn = args.t
    reference_seqs_fn = args.r
    gffread_path = args.g
    mm2_path = args.m

    # Check provided PATHs to executables
    if gffread_path:
        if not gffread_path.endswith("gffread"):
            raise ValueError("Invalid PATH provided to -g: {}".format(gffread_path))

    if mm2_path:
        if not mm2_path.endswith("minimap2"):
            raise ValueError("Invalid PATH provided to -m: {}".format(mm2_path))

    # Get the list of sequences to match
    target_fai = pysam.FastaFile(target_asm_fn)
    reference_fai = pysam.FastaFile(reference_asm_fn)

    target_seqs = set(target_fai.references)
    reference_seqs = set(reference_fai.references)

    target_fai.close()
    reference_fai.close()

    # Check if only specific sequences are to be matched
    if target_seqs_fn:
        with open(target_seqs_fn) as f:
            specified_target_seqs = set(f.read().rstrip().split("\n"))
            if not specified_target_seqs.issubset(target_seqs):
                raise ValueError("Sequences provided in '-t' are not a subset of the target assembly sequences.")
            target_seqs = specified_target_seqs

    if reference_seqs_fn:
        with open(reference_seqs_fn) as f:
            specified_reference_seqs = set(f.read().rstrip().split("\n"))
            if not specified_reference_seqs.issubset(reference_seqs):
                raise ValueError("Sequences provided in '-r' are not a subset of the reference assembly sequences.")
            reference_seqs = specified_reference_seqs

    # Check that there is the same number of target and reference sequences
    if len(target_seqs) != len(reference_seqs):
        raise RuntimeError("The number of target ({}) and reference ({}) sequences do not match.".format(len(target_seqs), len(reference_seqs)))

    # Get the current working directory
    output_path = os.getcwd() + "/"

    # Run gffread
    cmd = [
        gffread_path,
        "-w",
        output_path + "transcripts.fa",
        "-g",
        reference_asm_fn,
        reference_genes_fn
    ]
    run_oe(cmd, output_path + "gffread.out", output_path + "gffread.err")

    # Run minimap2
    cmd = [
        mm2_path,
        "-t",
        "8",
        "-cx",
        "splice",
        target_asm_fn,
        output_path + "transcripts.fa"
    ]
    run_oe(cmd, output_path + "alns.paf", output_path + "mm2.err")

    # Process the gff file
    log("INFO", "Processing the GFF file")
    tx2ref_seq = {}  # transcript -> reference sequence
    tx2strand = {}  # transcript -> strand

    db = gffutils.create_db(reference_genes_fn, reference_genes_fn + "_db", merge_strategy="create_unique", force=True, verbose=True)
    for feature in db.all_features():
        # Go through each gene and get the list of its transcripts
        if feature.featuretype == "gene":
            gene_id, chrom, strand = feature.id, feature.chrom, feature.strand
            if chrom in reference_seqs:
                transcripts = [i.id for i in db.children(gene_id, featuretype='mRNA')]
                if transcripts:
                    # Save the transcript if there is only one
                    if len(transcripts) == 1:
                        tx2ref_seq[transcripts[0]] = chrom
                        tx2strand[transcripts[0]] = strand
                    else:
                        # Pick the longest transcript
                        tx_lens = np.zeros(len(transcripts))
                        for i in range(len(transcripts)):
                            exons = [j for j in db.children(transcripts[i], featuretype='exon')]
                            for e in exons:
                                tx_lens[i] += e.stop - (e.start - 1)

                        tx2ref_seq[transcripts[int(np.argmax(tx_lens))]] = chrom
                        tx2strand[transcripts[int(np.argmax(tx_lens))]] = strand

    log("INFO", "Using {} representative reference transcripts".format(len(tx2ref_seq.keys())))

    # Process the mm2 alignments and create the bipartite graph
    log("INFO", "Building graph from minimap2 alignments")
    aligned_tx = set()
    G = nx.Graph()
    G.add_nodes_from(target_seqs, bipartite=0)
    G.add_nodes_from(reference_seqs, bipartite=1)

    paf_reader = PAFReader(output_path + "alns.paf")
    for aln in paf_reader.parse_alignments():
        if aln.query_header in tx2ref_seq and aln.ref_header in target_seqs:
            # Only consider alignments where > 85% of the transcript maps with mapq >= 10
            if aln.mapq >= 10 and (aln.query_end - aln.query_start)/aln.query_len > 0.85:
                aligned_tx.add(aln.query_header)
                this_ref_seq = tx2ref_seq[aln.query_header]
                this_strand = tx2strand[aln.query_header]
                if G.has_edge(this_ref_seq, aln.ref_header):
                    G[this_ref_seq][aln.ref_header]["weight"] -= 1
                    if this_strand == aln.strand:
                        G[this_ref_seq][aln.ref_header]["plus"] -= 1
                    else:
                        G[this_ref_seq][aln.ref_header]["minus"] -= 1
                else:
                    G.add_edge(this_ref_seq, aln.ref_header, weight=-1, plus=0, minus=0)
                    if this_strand == aln.strand:
                        G[this_ref_seq][aln.ref_header]["plus"] -= 1
                    else:
                        G[this_ref_seq][aln.ref_header]["minus"] -= 1

    log("INFO", "{:.2f}% of representative transcripts aligned with mapq >= 10 and coverage > 85%".format(((len(aligned_tx))/(len(tx2ref_seq.keys())))*100))

    # Write the graph as a gml file
    log("INFO", "Writing GML file")
    nx.readwrite.gml.write_gml(G, output_path + "graph.gml")

    # Check that each sequence has at least one edge
    connected_components = [i for i in nx.connected_components(G=G)]
    if len(connected_components) > 1:
        raise RuntimeError("Not enough information to assign sequences. Check GML file.")

    # Get the solution
    log("INFO", "Computing a matching solution")
    matching = nx.bipartite.minimum_weight_full_matching(G)

    # Sort the target sequences by length
    log("INFO", "Writing AGP output")
    target_fai = pysam.FastaFile(target_asm_fn)
    ref_lens = []
    for i in target_fai.references:
        ref_lens.append((target_fai.get_reference_length(i), i))
    ref_lens = sorted(ref_lens, reverse=True)

    # Write the solution in AGP format
    agp = AGPFile(output_path + "solution.agp", "w")
    agp.add_pragma()
    for i in ref_lens:
        length, header = i[0], i[1]
        if header in matching:
            strand = "+"
            if G[header][matching[header]]["minus"] < G[header][matching[header]]["plus"]:
                strand = "-"
            agp.add_seq_line(matching[header], 1, length, 1, "W", header, 1, length, strand)
        else:
            agp.add_seq_line(header, 1, length, 1, "W", header, 1, length, "+")

    agp.write()

    log("INFO", "use 'ragtag.py agp2fa' to produce a FASTA file")
    log("INFO", "Goodbye")


if __name__ == "__main__":
    main()
