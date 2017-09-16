"""
This script reformats a standard database of MLST alleles as per the format for genotyping in SRST2.
Run this script before using screen_genes_in_assemblies.py to format your MLST database for sequence typing.

Usage:
    cat mlst/*.fas | python reformat_mlst_db.py > srst2_mlst.fna
    cat mlst/*.fas | python reformat_mlst_db.py '-' > srst2_mlst.fna  # when the MLST delimiter is a dash

Python versions 2.7 and 3 compatible (This is my first script in Python 3)

Copyright (C) 2016-2017 Yu Wan <wanyuac@gmail.com, GitHub: https://github.com/wanyuac>
Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
First edition: 11 Aug 2016, latest edition: 11 Sep 2017
"""

from __future__ import print_function
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

if __name__ == "__main__":
    main()
    