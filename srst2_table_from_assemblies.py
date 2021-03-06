#!/usr/bin/env python
'''
SRST2 results from assemblies

This is a tool to screen for genes in a collection of assemblies and output
the results in a table which mimics those produced by SRST2.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import sys
import argparse
import subprocess
import os
import gzip
from distutils import spawn


def main():
    args = get_arguments()
    check_for_blast()
    check_file_exists(args.gene_db)
    check_algorithm(args.algorithm)

    unique_allele_symbols = determine_allele_symbol_uniqueness(args.gene_db)
    gene_db_name = os.path.splitext(os.path.basename(args.gene_db))[0]

    if args.report_new_consensus:
        new_consensus_path, new_consensus_name = os.path.split(args.report_new_consensus)
        if new_consensus_path and not os.path.exists(new_consensus_path):
            os.makedirs(new_consensus_path)
        new_consensus_alleles = open(args.report_new_consensus, 'w')
    if args.report_all_consensus:
        all_consensus_path, all_consensus_name = os.path.split(args.report_all_consensus)
        if all_consensus_path and not os.path.exists(all_consensus_path):
            os.makedirs(all_consensus_path)
        all_consensus_alleles = open(args.report_all_consensus, 'w')

    output_path, output_prefix = os.path.split(args.output)
    if output_path and not os.path.exists(output_path):
        os.makedirs(output_path)
    output_table = os.path.join(output_path, output_prefix + '__genes__' + gene_db_name + '__results.txt')

    all_clusters = set()
    all_results = {} # key = assembly_name, value = dict of cluster and allele
    for assembly in args.assemblies:
        assembly_name = remove_extension_from_assembly_file(os.path.basename(assembly))
        all_results[assembly_name] = {}
        blast_results = blast_assembly(assembly, args.gene_db, args.algorithm, unique_allele_symbols)
        filtered_blast_results = filter_blast_results(blast_results, args.min_coverage, args.max_divergence)
        best_hits = get_best_match_for_each_cluster(filtered_blast_results)
        for cluster, best_hit in best_hits.iteritems():
            all_clusters.add(cluster)
            query_name = best_hit[0]
            allele_name = best_hit[1]
            identity = best_hit[2]
            coverage = best_hit[3]
            hit_seq = best_hit[4]
            imperfect_match = identity < 100.0 or coverage < 100.0
            if imperfect_match:
                allele_name += '*'
            all_results[assembly_name][cluster] = (query_name, allele_name, hit_seq, imperfect_match)

    sorted_clusters = sorted(list(all_clusters))
    sorted_assemblies = sorted(list(all_results.keys()))

    out_file = open(output_table, 'w')
    out_file.write('\t'.join(['Sample'] + sorted_clusters))
    out_file.write('\n')
    for assembly in sorted_assemblies:
        assemly_line = assembly + '\t'
        results = []
        for cluster in sorted_clusters:
            if cluster not in all_results[assembly]:
                results.append('-')
            else:
                query_name = all_results[assembly][cluster][0]
                allele_name = all_results[assembly][cluster][1]
                hit_seq = all_results[assembly][cluster][2]
                imperfect_match = all_results[assembly][cluster][3]
                results.append(allele_name)
                if args.report_all_consensus:
                    full_query_name = query_name
                    if imperfect_match:
                        full_query_name += '.variant'
                    else:
                        full_query_name += '.consensus'
                    add_fasta_to_file(full_query_name + ' ' + assembly, hit_seq, all_consensus_alleles)
                if args.report_new_consensus and imperfect_match:
                    add_fasta_to_file(query_name + '.variant ' + assembly, hit_seq, new_consensus_alleles)
        out_file.write('\t'.join([assembly] + results))
        out_file.write('\n')


def check_file_exists(filename):
    if not os.path.isfile(filename):
        print('Error: could not load ' + filename, file=sys.stderr)
        quit()


def check_algorithm(algorithm):
    if algorithm != 'blastn' and algorithm != 'blastn-short' and algorithm != 'megablast' and algorithm != 'dc-megablast':
        print('Error: algorithm must be "blastn", "blastn-short", "megablast" or "dc-megablast"', file=sys.stderr)
        quit()


def check_for_blast():
    makeblastdb_path = spawn.find_executable('makeblastdb')
    blastn_path = spawn.find_executable('blastn')
    blast_installed = (makeblastdb_path != None and blastn_path != None)
    if not blast_installed:
        print('Error: could not find BLAST program', file=sys.stderr)
        quit()


def get_arguments():
    parser = argparse.ArgumentParser(description='SRST2 table from assemblies')
    parser.add_argument('--assemblies', nargs='+', type=str, required=True, help='Fasta file/s for assembled contigs')
    parser.add_argument('--gene_db', type=str, required=True, help='Fasta file for gene databases')
    parser.add_argument('--output', type=str, required=True, help='Output prefix for the table of results')
    parser.add_argument('--min_coverage', type=float, required=False, help='Minimum %%coverage cutoff for gene reporting (default 90)',default=90)
    parser.add_argument('--max_divergence', type=float, required=False, help='Maximum %%divergence cutoff for gene reporting (default 10)',default=10)
    parser.add_argument('--report_new_consensus', type=str, required=False, help='When matching alleles are not found, report the found alleles in this file')
    parser.add_argument('--report_all_consensus', type=str, required=False, help='Report all found alleles in this file')
    parser.add_argument('--algorithm', action="store", help="blast algorithm (blastn)", default="blastn")
    return parser.parse_args()


def blast_assembly(assembly, gene_db, algorithm, unique_allele_symbols):
    check_file_exists(assembly)

    # If the contigs are in a gz file, make a temporary decompressed FASTA file.
    if get_compression_type(assembly) == 'gz':
        new_assembly = assembly + '_temp_decompress.fasta'
        decompress_file(assembly, new_assembly)
        assembly = new_assembly
        temp_decompress = True
    else:
        temp_decompress = False

    # Make the BLAST database if it doesn't already exist.
    if not os.path.isfile(assembly + ".nin"):
        makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', assembly]
        process = subprocess.Popen(makeblastdb_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        if len(err) > 0:
            print('\nmakeblastdb encountered an error:', file=sys.stderr)
            print(err, file=sys.stderr)
            quit()

    blastn_command = ['blastn', '-task', algorithm, '-db', assembly, '-query', gene_db, '-outfmt', '6 qseqid pident qlen length sseq bitscore']
    process = subprocess.Popen(blastn_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if len(err) > 0:
        print('\nblastn encountered an error:', file=sys.stderr)
        print(err, file=sys.stderr)
        quit()

    blast_results = []

    for line in blast_results_iterator(out):
        line_parts = line.split('\t')
        if len(line_parts) < 5:
            continue
        query_name = line_parts[0]
        identity = float(line_parts[1])
        query_length = float(line_parts[2])
        hit_length = float(line_parts[3])
        coverage = 100.0 * hit_length / query_length
        hit_seq = line_parts[4]
        bit_score = float(line_parts[5])
        query_name_parts = query_name.split()[0].split('__')
        if len(query_name_parts) < 4:
            print('Error: gene database names must be in the following format:', file=sys.stderr)
            print('[clusterUniqueIdentifier]__[clusterSymbol]__[alleleSymbol]__[alleleUniqueIdentifier]', file=sys.stderr)
            quit()
        cluster_name = query_name_parts[1]
        allele_name = query_name_parts[2]
        allele_id = query_name_parts[3]

        if not unique_allele_symbols:
            allele_name += '_' + allele_id

        blast_results.append((query_name, cluster_name, allele_name, identity, coverage, hit_seq, bit_score))

    # If we've been working on a temporary decompressed file, delete it now.
    if temp_decompress:
        os.remove(assembly)

    return blast_results


def filter_blast_results(blast_results, coverage_cutoff, max_divergence):
    identity_cutoff = 100.0 - max_divergence
    filtered_blast_results = []
    for result in blast_results:
        identity = result[3]
        coverage = result[4]
        if identity >= identity_cutoff and coverage >= coverage_cutoff:
            filtered_blast_results.append(result)
    return filtered_blast_results


def get_best_match_for_each_cluster(filtered_blast_results):
    best_hits = {}
    for result in filtered_blast_results:
        query_name = result[0]
        cluster_name = result[1]
        allele_name = result[2]
        identity = result[3]
        coverage = result[4]
        hit_seq = result[5].replace('-', '')
        bit_score = result[6]
        this_hit = (query_name, allele_name, identity, coverage, hit_seq, bit_score)

        # If a hit is the first one for its cluster, then it is by definition
        # the best hit.
        if cluster_name not in best_hits:
            best_hits[cluster_name] = this_hit

        else:
            current_best_identity = best_hits[cluster_name][2]
            current_best_coverage = best_hits[cluster_name][3]
            current_best_length = len(best_hits[cluster_name][4])
            current_best_bit_score = best_hits[cluster_name][5]

            hit_is_perfect = (identity == 100.0 and coverage == 100.0)
            current_best_is_perfect = (current_best_identity == 100.0 and current_best_coverage == 100.0)

            # If this hit is perfect and the current best isn't perfect, then
            # this is the new best hit.
            if hit_is_perfect and not current_best_is_perfect:
                best_hits[cluster_name] = this_hit

            # If this hit is perfect and the current best is also perfect, then
            # this is the new best hit only if it is longer.
            elif hit_is_perfect and current_best_is_perfect and len(hit_seq) > current_best_length:
                best_hits[cluster_name] = this_hit

            # If neither this hit nor the current best are perfect, this is the
            # new best only if it has the higher bit score.
            elif not hit_is_perfect and not current_best_is_perfect and bit_score > current_best_bit_score:
                best_hits[cluster_name] = this_hit

    return best_hits


# http://stackoverflow.com/questions/3054604/iterate-over-the-lines-of-a-string
def blast_results_iterator(results):
    prevnl = -1
    while True:
      nextnl = results.find('\n', prevnl + 1)
      if nextnl < 0: break
      yield results[prevnl + 1:nextnl]
      prevnl = nextnl


def remove_extension_from_assembly_file(assembly_filename):
    assembly_filename = rchop(assembly_filename, '.fasta')
    assembly_filename = rchop(assembly_filename, '.fa')
    return assembly_filename


# http://stackoverflow.com/questions/3663450/python-remove-substring-only-at-the-end-of-string
def rchop(thestring, ending):
    if thestring.endswith(ending):
        return thestring[:-len(ending)]
    return thestring


def add_fasta_to_file(name, seq, file):
    file.write('>')
    file.write(name)
    file.write('\n')
    file.write(seq)
    file.write('\n')


def determine_allele_symbol_uniqueness(gene_db_filename):
    '''
    This function determines whether any two alleles in the gene database have
    the same allele identifier.  If this is the case, then every allele in the
    results will have its allele identifier included.
    It returns True if all allele names are unique and false if there is at
    least one duplicate.
    This mimics the behaviour of SRST2, which does the same.
    '''
    allele_names = set()
    gene_db = open(gene_db_filename, 'r')
    for line in gene_db:
        if not line.startswith('>'):
            continue
        name_parts = line.split()[0].split('__')
        allele_name = name_parts[2]
        if allele_name in allele_names:
            return False
        allele_names.add(allele_name)
    return True


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('cannot use zip format - use gzip instead')
    return compression_type


def decompress_file(in_file, out_file):
    with gzip.GzipFile(in_file, 'rb') as i, open(out_file, 'wb') as o:
        s = i.read()
        o.write(s)


if __name__ == '__main__':
    main()
