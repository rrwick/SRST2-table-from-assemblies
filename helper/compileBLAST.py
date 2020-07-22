#!/usr/bin/env python

"""
Compile and parse nucleotide BLAST+ outputs (outfmt 6) for a number of genomes. The output must be
generated using command line:
blastn -query [FASTA file] -db [Database name] -out [sample name][delimiter][suffix].tsv -task megablast \
    -evalue [1e-5] -perc_identity [70] -max_target_seqs [10] \
    -outfmt '6 qseqid sseqid slen pident qcovhsp length mismatch gapopen qstart qend sstart send sstrand evalue bitscore sseq'

Input: Files of BLAST+ outputs in output format 6.
Outputs:
    A tab-delimited table compiling hits in all sample genomes
    FASTA files of subject DNA sequences, one file per query sequence
    FASTA files of translated subject DNA sequences based on a given codon table (for bacterial genes: 11), one file per query sequence
Example command:
    python compileBLAST.py --input *__megaBLAST.tsv --delimiter '__' --genes 'gene1,gene2,gene3' --output study1 --translate --codon_table 11 --add_sample_name > compile_blast.log
Dependencies: Python 3, BioPython

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 11 June 2020; the latest modification: 22 July 2020
"""

import os
import re
from collections import namedtuple
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_arguments():
    parser = ArgumentParser(description = 'Compile and parse nucleotide BLAST outputs')
    parser.add_argument('--input', '-i', dest = 'input', nargs = '+', required = True, type = str,\
        help = 'Tab-delimited BLAST outputs in format 6')
    parser.add_argument('--delimiter', '-d', dest = 'delimiter', type = str, required = True, \
        help = 'Delimiter character(s) separating the sample name and the suffix in every input filename')
    parser.add_argument('--genes', '-g', dest = 'genes', type = str, required = True, \
        help = 'Comma delimited vector of gene names (BLAST queries) in BLAST output files')
    parser.add_argument('--output', '-o', dest = 'output', type = str, required = False, default = 'blast', \
        help = 'Common prefix of output files. Default: blast')
    parser.add_argument('--translate', '-t', dest = 'translate', action = "store_true", required = False,\
        help = 'Translate DNA sequences into protein sequences')
    parser.add_argument('--codon_table', '-c', dest = 'codon_table', required = False, type = int, default = 11, \
        help = 'Index of codon tables for translating coding sequences into protein sequences')
    parser.add_argument('--add_sample_name', '-a', action = "store_true", required = False,\
        help = 'Append a sample name to each sequence ID')

    return parser.parse_args()


HIT_ATTRS = ['qseqid', 'sseqid', 'slen', 'pident', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',\
             'sstart', 'send', 'sstrand', 'evalue', 'bitscore']  # 15 attribute fields


def main():
    args = parse_arguments()

    # Create output files, keep them open, and save their handles
    genes = args.genes.split(',')
    out_tsv = open(args.output + '__aln.tsv', 'w')  # A summary of alignments
    out_tsv.write('\t'.join(['sample'] + HIT_ATTRS) + '\n')  # The write method does not automatically append a newline character to the end of the output line.
    out_fna_handles = {}
    if args.translate:
        genes_prot = {}  # Protein names from gene names
        out_faa_handles = {}

    for g in genes:
        p = g[0].upper() + g[1 : ]  # Function capitalize() changes other letters to lower case
        fna_out = open('%s__%s.fna' % (args.output, g), 'w')
        out_fna_handles[g] = fna_out
        
        if args.translate:
            genes_prot[g] = p  # Protein name
            faa_out = open('%s__%s.faa' % (args.output, p), 'w')
            out_faa_handles[p] = faa_out

    # Pool and parse BLAST outputs
    for tsv in args.input:
        sample_name = extract_sample_name(tsv, args.delimiter)
        print('Parsing BLAST output of sample %s.' % sample_name)
        with open(tsv, 'r') as f:
            lines = f.read().splitlines()
        
        # Parse each line and write data into output files
        for line in lines:
            hit = parse_blast_line(line)
            g = hit.qseqid  # Name of the query gene
            if g in genes:
                # Write alignment information
                out_tsv.write('\t'.join([sample_name, hit.qseqid, hit.sseqid, hit.slen, hit.pident, hit.qcovhsp, hit.length,\
                              hit.mismatch, hit.gapopen, hit.qstart, hit.qend, hit.sstart, hit.send, hit.sstrand, hit.evalue,\
                              hit.bitscore]) + '\n')
                
                # Prepare a DNA record
                seq_id = '.'.join([g, sample_name]) if args.add_sample_name else g
                descr = '|'.join([g, sample_name, hit.sseqid, hit.sstart + '-' + hit.send, hit.sstrand])  # It changes later when args.translate = True.
                seq_dna = SeqRecord(Seq(hit.sseq), id = seq_id, name = '', description = '')  # The description will be filled later.

                # Write the protein sequence when it is required
                if args.translate:
                    # Translate DNA till the first stop codon
                    try:
                        seq_prot = seq_dna.translate(table = args.codon_table, id = seq_id[0].upper() + seq_id[1 : ],\
                                   description = '', to_stop = True, cds = False)  # Set cds = True to check error: 'First codon is not a start codon'.
                    except KeyError:
                        print('Warning: sequence of %s in %s cannot be translated.' % (g, sample_name))
                        seq_prot = SeqRecord(Seq(''), id = '', name = '', description = '')
                    
                    if len(seq_prot.seq) > 0:  # Sometimes partial sequences do not have any protein product.
                        descr = descr + '|CDS'  # This sequence can be translated into a protein (or polypeptide) sequence
                        seq_prot.description = descr
                        write_fasta(out_faa_handles[genes_prot[g]], seq_prot)
                    else:
                        descr = descr + '|NA'  # No translation is available and to be written, but change the sequence description

                # Finally, write the DNA sequence
                seq_dna.description = descr
                write_fasta(out_fna_handles[g], seq_dna)

    # Close output files
    out_tsv.close()
    for g in genes:
        out_fna_handles[g].close()
        if args.translate:
            out_faa_handles[genes_prot[g]].close()

    return


def parse_blast_line(line):
    """
    Parse each line of BLAST outputs into a Hit object (namedtuple)
    """
    Hit = namedtuple("Hit", HIT_ATTRS + ['sseq'])
    fields = line.split('\t')

    # Use a regular expression to remove gaps ('-') in the DNA sequence.
    # See stackoverflow.com/questions/3939361/remove-specific-characters-from-a-string-in-python
    hit = Hit(qseqid = fields[0], sseqid = fields[1], slen = fields[2], pident = fields[3], qcovhsp = fields[4],\
              length = fields[5], mismatch = fields[6], gapopen = fields[7], qstart = fields[8], qend = fields[9],\
              sstart = fields[10], send = fields[11], sstrand = '+' if fields[12] == 'plus' else '-', evalue = fields[13],\
              bitscore = fields[14], sseq = re.sub('-', '', fields[15]))

    return hit


def extract_sample_name(filename, delim):
    """
    Extract sample name from a file name.
    """
    b = os.path.basename(filename)
    sample_name = b.split(delim)[0]

    return sample_name


def write_fasta(file_handle, seq_rec):
    """
    Write sequences without wrapping.
    Method write_record of FastaIO.FastaWriter does not work when multiple file handles are kept open.
    """
    print('>%s %s' % (seq_rec.id, seq_rec.description), file = file_handle)
    print(seq_rec.seq, file = file_handle)

    return


if __name__ == '__main__':
    main()
