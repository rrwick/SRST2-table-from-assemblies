"""
This script reformats a standard database of MLST alleles as per the format for genotyping in SRST2.
Run this script before using Ryan's script srst2_table_from_assemblies.py (https://github.com/rrwick/SRST2-table-from-assemblies)
for sequence typing.

Usage:
    cat mlst/*.fas | python reformat_mlst_db.py > srst2_mlst.fna
    cat mlst/*.fas | python reformat_mlst_db.py '-' > srst2_mlst.fna  # when the MLST delimiter is a dash

Author: Yu Wan (wanyuac@gmail.com, GitHub: https://github.com/wanyuac)
Edition history: 11/8/2016
Python v3.5.2 (This is my first script in Python 3)
License: GNU GPL 2.1
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_mlst_delimiter():
    if len(sys.argv) > 1:
        delimiter = sys.argv[1]
    else:
        delimiter = "_"  # the default
    return delimiter

def main():
    mlst_delim = get_mlst_delimiter()
    cluster_id = -1
    seq_num = 0
    gene_list = {}
    genes = []

    # read the whole MLST database from the stdin
    for seq in SeqIO.parse(sys.stdin, "fasta"):  # read the input FASTA file from stdin
        gene, allele_id = seq.id.split(mlst_delim)
        seq_num += 1
        if not gene in genes:  # if this gene is not recorded yet
            cluster_id += 1  # assign a new cluster ID
            gene_list[gene] = cluster_id  # push the new item {gene : cluster_id} into the dictionary
            genes = list(gene_list.keys())
        seq.id = "__".join([str(gene_list[gene]), gene, allele_id, str(seq_num)])  # assign a new sequence ID
        print(seq.format("fasta"))  # write the current sequence to the stdout

if __name__ == '__main__':
    main()
    