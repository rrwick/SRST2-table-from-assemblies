#!/usr/bin/env python
"""
SRST2 results from assemblies (contigs)

This is a tool to screen for genes in a collection of contigs and output the results in a table which mimics those produced by SRST2.
This script is able to process multiple input FASTA files.

Subject sequences (the BLAST database): the assembly
Query sequences: an SRST2-formatted reference gene database

Python version: 3

Copyright (C) 2015-2017 Ryan Wick <rrwick@gmail.com>, Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
Latest edition: 8-9 Sep 2017
"""

import os
import sys
sys.dont_write_bytecode = True  # Do not write .pyc files (in the __pycache__ directory) on the import of source modules. Use it before "from parseBLAST import ...".
import gzip
import argparse
import subprocess
from distutils import spawn
from collections import namedtuple
from parseBLAST import Hit, Hap, Cluster, Assembly  # for classes that process BLAST outputs


def get_arguments():
    parser = argparse.ArgumentParser(description = "SRST2 table from assemblies")
    parser.add_argument("--assemblies", nargs = "+", type = str, required = True, help = "Fasta file/s for assembled contigs")
    parser.add_argument("--gene_db", type = str, required = True, help = "Fasta file for gene databases")
    parser.add_argument("--prefix", type = str, required = False, default = "BLAST", help = "Output prefix for the table of results")
    parser.add_argument("--suffix", type = str, required = False, default = ".fasta", help = "Characters to be chopped off from the end of every assembly name in order to get a sample name")
    parser.add_argument("--outdir", type =str, required = False, default = ".", help = "Output directory for the table of results")
    parser.add_argument("--min_coverage", type = float, required = False, default = 90.0, help = "Minimum %%coverage cutoff for gene reporting (default 90)")
    parser.add_argument("--max_divergence", type = float, required = False, default = 10.0, help = "Maximum %%divergence cutoff for gene reporting (default 10)")
    parser.add_argument("--report_new_consensus", action = "store_true", required = False, help = "Configure it to save consensus sequences of variants")
    parser.add_argument("--report_all_consensus", action = "store_true", required = False, help = "Configure it to save all consensus sequences")
    parser.add_argument("--algorithm", action = "store", required = False, help = "blast algorithm (blastn)", default = "megablast")
    parser.add_argument("--mlst", action = "store_true", required = False, help = "Turn it on to find MLST genes")
    parser.add_argument("--incl_alt", action = "store_true", required = False, help = "Flag it to include all putative alternative calls for each allele")
    parser.add_argument("--max_overlapping_nt", type = int, required = False, default = 0, help = "Maximal number of overlapping nucleotides allowed to treat two hits as physically separate")
    parser.add_argument("--del_blast", action = "store_true", required = False, help = "Flag it to delete the text file of BLAST outputs")
    return parser.parse_args()


def main():
    args = get_arguments()
    
    # Check and configure the runtime environment ###############
    check_for_blast()
    check_file_exists(args.gene_db)
    check_algorithm(args.algorithm)
    if args.outdir and not os.path.exists(args.outdir):  # This logical statement becomes false when output_path = "".
        os.makedirs(args.outdir)
    report_best_hits = not args.incl_alt
    unique_allele_symbols = determine_allele_symbol_uniqueness(args.gene_db)  # return True of False
    gene_db_name = os.path.splitext(os.path.basename(args.gene_db))[0]

    # Run BLAST to determine the presence of alleles ###############
    assemblies = {} # key = assembly_name, value = an Assembly object
    assembly_names = []
    cluster_names = set()
    for assembly in args.assemblies:  # for each FASTA file
        assembly_name = rchop(os.path.basename(assembly), args.suffix)  # eg. ERR08090_spades.fasta => ERR08090 given args.suffix = "_spades.fasta".
        assembly_names.append(assembly_name)
        
        # run BLAST and parse outputs line by line, returning a dictionary of clusters of the current assembly
        clusters = blast_assembly(assembly = assembly, gene_db = args.gene_db, algorithm = args.algorithm, \
                                  unique_allele_symbols = unique_allele_symbols, mlst_run = args.mlst, \
                                  blast_out = os.path.join(args.outdir, "__".join([assembly_name, args.prefix, \
                                                                                   "blast", gene_db_name, "results.txt"])), \
                                  del_blast = args.del_blast)
        
        # filter hits of each cluster based on thresholds of query coverage and nucleotide divergence
        for _, cluster in clusters.items():
            cluster.filter_hits(args.min_coverage, args.max_divergence)
        
        # organise remaining clusters into an Assembly object
        if not assembly_name in list(assemblies.keys()):
            assemblies[assembly_name] = Assembly(assembly_name, clusters, True)  # create an Assembly object using non-empty clusters of filtered hits
        
        # add cluster names of the current assembly to the overall set of names with duplicates removed
        cluster_names.update(assemblies[assembly_name].cluster_names)  # use the update method to append iterable values into a set
        
    sorted_cluster_names = sorted(list(cluster_names))
    sorted_assembly_names = sorted(assembly_names)
    
    # Set up output files ###############
    if args.report_new_consensus:
        new_consensus_alleles = open(os.path.join(args.outdir, args.prefix + ".new_consensus_alleles.fasta"), "w")  # file handle for consensus sequences of called variants
    if args.report_all_consensus:
        all_consensus_alleles = open(os.path.join(args.outdir, args.prefix + ".all_consensus_alleles.fasta"), "w")  # file handle for consensus sequences of all called alleles    
    genotype_file = open(os.path.join(args.outdir, args.prefix + "__genes__" + gene_db_name + "__results.txt"), "w")  # mimic the per-sample result from SRST2
    
    # Write results into files ################
    genotype_file.write("\t".join(["Sample"] + sorted_cluster_names) + "\n")  # write the header line of the (merged) genotype file (for all samples)
    
    for assembly in sorted_assembly_names:
        line_fields = [assembly]  # value for the first column in the current row
        if report_best_hits:
            """ Only report the best hit for each cluster, which is exactly the behaviour of SRST2. """
            best_hits = assemblies[assembly].find_best_hits()  # return a dictionary of Hit objects, each for a cluster; may return an empty dictionary but it does affect the algorithm.
            clusters_present = list(best_hits.keys())  # clusters present in the current assembly; return [] when best_hits = {}
            for cluster in sorted_cluster_names:
                if cluster in clusters_present:  # always return False when clusters_present = []
                    hit = best_hits[cluster]  # a Hit object, which is the best hit within the current cluster
                    hit_seq = hit.hit_seq.replace("-", "")  # remove signs of deletions from the consensus sequence
                    if len(hit) > 0:  # len(hit) = 0 when there is no qualified (above the thresholds for query coverage and nucleotide identity) allele calls at all.
                        line_fields.append(hit.allele)  # add an allele name for the current column
                        if args.report_all_consensus:
                            full_query_name = hit.query  # such as 80__TetG_Tet__TetG__632
                            if hit.perfect_match:
                                full_query_name += ".consensus"
                            else:
                                full_query_name += ".variant"
                            add_fasta_to_file(full_query_name + " " + assembly, hit_seq, all_consensus_alleles)
                        if args.report_new_consensus and not hit.perfect_match:  # only write consensus sequences of variants
                            add_fasta_to_file(hit.query + ".variant " + assembly, hit_seq, new_consensus_alleles)
                    else:
                        line_fields.append("-")
                else:  # when the current assembly does not have this cluster called
                    line_fields.append("-")
        else:
            """ Alternatively, report all physically separated hits per cluster, which yields multiple allele calls for each cluster. """
            all_hits = assemblies[assembly].find_all_copies(args.max_overlapping_nt)  # return a dictionary {cluster : [hit1, ..., hitn]}
            clusters_present = list(all_hits.keys())  # equals [] when all_hits is an empty dictionary
            for cluster in sorted_cluster_names:
                if cluster in clusters_present:  # False when clusters_present is an empty list
                    haps = merge_hits_by_haplotypes(all_hits[cluster])  # Command "all_hits[cluster]" returns a list of non-overlapping valid hits within this cluster.
                    if len(haps) > 0:
                        # determine the content for the current item in the genotype matrix
                        allele_names = []  # for the current column
                        
                        # define a counter for variants of each allele
                        all_queries = set()
                        for hap in haps:  # hap: a Hap namedtuple
                            all_queries.add(hap.rep_hit.query)
                        var_counts = {k : 0 for k in list(all_queries)}  # define a dictionary with default values given a list of keys; {query : variant count}
                        
                        for hap in haps:  # go through every haplotype of the current cluster
                            if hap.copy_num > 1:
                                allele_name_suffix = "[" + str(hap.copy_num) + "]"
                            else:
                                allele_name_suffix = ""
                            allele_names.append(hap.rep_hit.allele + allele_name_suffix)  # such as "AadA24_1609*[2]"
                            
                            # count the number of variants
                            if not hap.rep_hit.perfect_match:
                                var_counts[hap.rep_hit.query] += 1
                            
                            # write consensus sequences into a file
                            hap_seq = hap.rep_hit.hit_seq.replace("-", "")
                            if args.report_all_consensus:
                                full_query_name = hap.rep_hit.query  # such as 80__TetG_Tet__TetG__632
                                if hap.rep_hit.perfect_match:
                                    full_query_name += ".consensus"
                                else:
                                    """
                                    Some alleles of the same cluster may have different names. Therefore, variant count only
                                    appends to variants that share the same allele name.
                                    """
                                    n = var_counts[hap.rep_hit.query]
                                    if n > 1:
                                        full_query_name += ".variant" + str(n)
                                    else:
                                        full_query_name += ".variant"
                                add_fasta_to_file(full_query_name + " " + assembly, hap_seq, all_consensus_alleles)
                            if args.report_new_consensus and not hap.rep_hit.perfect_match:  # only write consensus sequences of variants
                                n = var_counts[hap.rep_hit.query]
                                if n > 1:
                                    full_query_name = hap.rep_hit.query + ".variant" + str(n)
                                else:
                                    full_query_name = hap.rep_hit.query + ".variant"
                                add_fasta_to_file(full_query_name + " " + assembly, hap_seq, new_consensus_alleles)
                        line_fields.append(",".join(allele_names))  # eg. "AadA24_1609*[2],AadA3*,AadA23[3]"
                    else:
                        line_fields.append("-")
                else:
                    line_fields.append("-")
                
        genotype_file.write("\t".join(line_fields) + "\n")  # make one line for each assembly (sample) in the genotype profile
    
    genotype_file.close()
    if args.report_all_consensus:
        all_consensus_alleles.close()
    if args.report_new_consensus:
        new_consensus_alleles.close()
    return
    

# Run BLAST and parse its output for a given assembly
def blast_assembly(assembly, gene_db, algorithm, unique_allele_symbols, mlst_run, blast_out, del_blast):
    check_file_exists(assembly)
        
    # If the contigs are in a gz file, make a temporary decompressed FASTA file.
    if get_compression_type(assembly) == "gz":
        new_assembly = assembly + "_temp_decompress.fasta"
        decompress_file(assembly, new_assembly)
        assembly = new_assembly
        temp_decompress = True
    else:
        temp_decompress = False

    # Make the BLAST database if it doesn"t already exist.
    if not os.path.isfile(assembly + ".nin"):
        makeblastdb_command = ["makeblastdb", "-dbtype", "nucl", "-in", assembly]
        process = subprocess.Popen(makeblastdb_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        _, err = process.communicate()  # return bytes not strings
        err = err.decode()  # convert bytes into a string
        if len(err) > 0:
            print("\nmakeblastdb encountered an error:", file = sys.stderr)
            print(err, file = sys.stderr)
            quit()
            
    """
    Run BLAST and fetch results from the subprocess handel for an assembly. Here, we align every allele in the database as a query against every subject ie. contigs.
    sseq: aligned part of the subject sequence. In this application, it is the aligned part in a contig against a query sequence (a reference allele).
    length: alignment length. The output may contain multiple lines for each allele.
    """
    header = "qseqid sseqid sstart send pident qlen length bitscore sseq"
    n_columns = header.count(" ") + 1
    blastn_command = ["blastn", "-task", algorithm, "-db", assembly, "-query", gene_db, "-outfmt", "6 " + header]
    process = subprocess.Popen(blastn_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = process.communicate()  # obtain outcomes (out and err) as two long strings (in bytes class though)
    out = out.decode()
    err = err.decode()
    if len(err) > 0:
        print("\nblastn encountered an error:", file = sys.stderr)
        print(err, file = sys.stderr)
        quit()
    
    if not del_blast:
        with open(blast_out, "w") as f:
            f.write(header.replace(" ", "\t") + "\n")
            f.write(out)

    # Parse BLAST results line by line
    clusters = {}  # a dictionary of Cluster objects
    for line in blast_results_iterator(out):  # Each line is a hit.
        n = line.count("\t") + 1
        if n < n_columns:  # ignore "no-hit" lines
            continue
        else:
            query_name, subject_name, subject_start, subject_end, identity, query_length, hit_length, bit_score, hit_seq = line.split("\t")
            if query_name.count("__") < 3:
                print("Error: gene database names must be in the following format:", file = sys.stderr)
                print("[clusterUniqueIdentifier]__[clusterSymbol]__[alleleSymbol]__[alleleUniqueIdentifier]", file = sys.stderr)
                sys.exit(1)
                
            """
            Extract fields from the query name. Eg. 74__TetD_Tet__TetD__1047 => [74, TetD_Tet, TetD, 1047].
            cluster_name: gene name; seqid: an integer distinguishing alleles of the same allele_name.            
            """
            _, cluster_name, allele_name, seq_id = query_name.split("__")
            if not unique_allele_symbols and not mlst_run:
                allele_name += "_" + seq_id  # eg. TetD => TetD_1047
                
            # determine whether the current hit represents a perfect match
            identity = float(identity)
            coverage = float(hit_length) / float(query_length) * 100.0
            is_imperfect_match = identity < 100.0 or coverage < 100.0
            if is_imperfect_match:
                allele_name += "*"  # eg. TetD_1047 => TetD_1047*
                
            # initialise an empty cluster object if the current cluster is new
            if not cluster_name in list(clusters.keys()):
                clusters[cluster_name] = Cluster(cluster_name, [])
            
            """
            Append each line of the BLAST output as a named tuple to a list. The named tuple class Hit is defined
            in the module parseBLAST. The query coverage may be greater than 100% when there are any insertions.
            """
            hit = Hit(query = query_name, cluster = cluster_name, allele = allele_name, \
                      query_length = int(query_length), coverage = coverage, contig = subject_name, \
                      hit_length = int(hit_length), start = int(subject_start), end = int(subject_end), \
                      identity = identity, bit_score = float(bit_score), perfect_match = not is_imperfect_match, \
                      hit_seq = hit_seq)
            clusters[cluster_name].add_hit(hit)  # add this hit into a corresponding cluster
        
    # If we"ve been working on a temporary decompressed file, delete it now.
    if temp_decompress:
        os.remove(assembly)

    return clusters


# Chop a long string into elements of a generator based on newline characters. Ref: http://stackoverflow.com/questions/3054604/iterate-over-the-lines-of-a-string
def blast_results_iterator(results):
    prevnl = -1  # not found
    while True:  # keep searchng for newline characters until reaching the last one
        """
        The find method "returns the lowest index in the string where substring sub is found within the slice s[start : end]".
        Here, the search starts from the position 0 in the first iteration.
        """
        nextnl = results.find("\n", prevnl + 1)
        if nextnl < 0:  # after the last "\n" at the end of the string, or results do not have multiple lines at all
            break  # terminate the while loop
        
        """
        The yield command makes the function to return a "generator" (data type), which can only be iterated once.
        Every element in this generator corresponds to a line in the BLAST output.
        """
        yield results[prevnl + 1 : nextnl]  # extract all characters within a certain region defined by (prevnl + 1) and nextnl
        prevnl = nextnl
    return


def merge_hits_by_haplotypes(hits):
    """ Identify haplotypes (hap) that a list of Hit objects (hits) represents """
    haps = []  # a list of Hap objects
    n = len(hits)
    while n > 1:
        hap = hits[0]  # a Hit object
        hits = [h for h in hits[1 : ] if h.hit_seq != hap.hit_seq]  # filter out hits sharing the same sequence as the current hap
        n_filtered = len(hits)  # may equal zero when new hits = []
        haps.append(Hap(rep_hit = hap, copy_num = n - n_filtered))  # Hap.rep_hit is a namedtuple as a value of another namedtuple.
        n = n_filtered
    # Python 3 does not allow the while...elif... syntax.
    if n == 1:  # for the last hit after filtering or hits consists of a single hit at the beginning
        haps.append(Hap(rep_hit = hits[0], copy_num = 1))
    return haps        


# http://stackoverflow.com/questions/3663450/python-remove-substring-only-at-the-end-of-string
def rchop(thestring, ending):
    if thestring.endswith(ending):
        return thestring[ : -len(ending)]
    return thestring


def add_fasta_to_file(name, seq, file):  # Here, file is a handle.
    file.write(">" + name + "\n")
    file.write(seq + "\n")
    return


def determine_allele_symbol_uniqueness(gene_db_filename):
    """
    This function determines whether any two alleles in the gene database have
    the same allele identifier.  If this is the case, then every allele in the
    results will have its allele identifier included.
    It returns True if all allele names are unique and false if there is at
    least one duplicate.
    This mimics the behaviour of SRST2, which does the same.
    """
    allele_names = set()
    gene_db = open(gene_db_filename, "r")
    for line in gene_db:
        if not line.startswith(">"):
            continue
        name_parts = line.split()[0].split("__")
        allele_name = name_parts[2]
        if allele_name in allele_names:
            return False
        allele_names.add(allele_name)
    return True


def check_for_blast():
    makeblastdb_path = spawn.find_executable("makeblastdb")
    blastn_path = spawn.find_executable("blastn")
    blast_installed = (makeblastdb_path != None and blastn_path != None)
    if not blast_installed:
        sys.exit("Error: could not find BLAST program")
    return


def check_file_exists(filename):
    if not os.path.isfile(filename):
        sys.exit("Error: could not load " + filename)
    return


def check_algorithm(algorithm):
    if not algorithm in ["blastn", "blastn-short", "megablast", "dc-megablast"]:
        sys.exit("Error: algorithm must be blastn, blastn-short, megablast or dc-megablast")
    return


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {"gz": (b"\x1f", b"\x8b", b"\x08"),
                  "bz2": (b"\x42", b"\x5a", b"\x68"),
                  "zip": (b"\x50", b"\x4b", b"\x03", b"\x04")}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in list(magic_dict.items()):
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("cannot use bzip2 format - use gzip instead")
    if compression_type == "zip":
        sys.exit("cannot use zip format - use gzip instead")
    return compression_type


def decompress_file(in_file, out_file):  # output_file is another file handle.
    with gzip.GzipFile(in_file, "rb") as i, open(out_file, "wb") as o:
        s = i.read()
        o.write(s)
    return


if __name__ == "__main__":  # Execute the main function when this script is not imported as a module.
    main()
