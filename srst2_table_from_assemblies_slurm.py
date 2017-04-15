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
import os

sys.dont_write_bytecode = True
from srst2_table_from_assemblies import check_file_exists
from srst2_table_from_assemblies import remove_extension_from_assembly_file
from srst2_table_from_assemblies import check_algorithm


def main():
    args = get_arguments()

    check_file_exists(args.gene_db)
    args.gene_db = os.path.abspath(args.gene_db)

    if args.algorithm:
        check_algorithm(args.algorithm)

    if not args.rundir:
        args.rundir = os.getcwd()

    if args.script:
        if args.script.endswith('srst2_table_from_assemblies.py'):
            script_path = args.script
        else:
            script_path = os.path.join(args.script, 'srst2_table_from_assemblies.py')
        check_file_exists(script_path)
    else:
        script_path_cwd = os.path.join(os.getcwd(), 'srst2_table_from_assemblies.py')
        if os.path.isfile(script_path_cwd):
            script_path = script_path_cwd
        else:
            script_path_file_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'srst2_table_from_assemblies.py')
            if os.path.isfile(script_path_file_directory):
                script_path = script_path_file_directory
            else:
                print('Error: could not find srst2_table_from_assemblies.py')
                quit()

    for assembly_filename in args.assemblies:
        assembly_name = remove_extension_from_assembly_file(os.path.basename(assembly_filename))

        output_path, output_name = os.path.split(args.output)
        if not output_path:
            output_path = args.rundir
        output_name_and_path = os.path.join(output_path, assembly_name + '_' + output_name)
    
        cmd = '#!/bin/bash'
        cmd += '\n#SBATCH -p sysgen'
        cmd += '\n#SBATCH --job-name=srst2_table_' + assembly_name
        cmd += '\n#SBATCH --ntasks=1'
        cmd += '\n#SBATCH --mem-per-cpu=' + args.memory
        cmd += '\n#SBATCH --time=' + args.walltime
        cmd += '\ncd ' + args.rundir
        cmd += '\nmodule load BLAST+/2.2.30-vlsci_intel-2015.08.25-Python-2.7.10'
        cmd += '\nmodule load Python/2.7.10-vlsci_intel-2015.08.25-SG'
        cmd += '\n' + script_path
        cmd += ' --assemblies ' + assembly_filename
        cmd += ' --gene_db ' + args.gene_db
        cmd += ' --output ' + output_name_and_path
        if args.min_coverage:
            cmd += ' --min_coverage ' + str(args.min_coverage)
        if args.max_divergence:
            cmd += ' --max_divergence ' + str(args.max_divergence)
        if args.algorithm:
            cmd += ' --algorithm ' + args.algorithm
        if args.mlst:
            cmd += ' --mlst'

        if args.report_new_consensus:
            new_consensus_path, new_consensus_name = os.path.split(args.report_new_consensus)
            if not new_consensus_path:
                new_consensus_path = args.rundir
            new_consensus_name_and_path = os.path.join(new_consensus_path, assembly_name + '_' + new_consensus_name)
            cmd += ' --report_new_consensus ' + new_consensus_name_and_path
        if args.report_all_consensus:
            all_consensus_path, all_consensus_name = os.path.split(args.report_all_consensus)
            if not all_consensus_path:
                all_consensus_path = args.rundir
            all_consensus_name_and_path = os.path.join(all_consensus_path, assembly_name + '_' + all_consensus_name)
            cmd += ' --report_all_consensus ' + all_consensus_name_and_path

        print(cmd)
        print()
        os.system('echo "' + cmd + '" | sbatch')
        print()



def get_arguments():
    parser = argparse.ArgumentParser(description='SRST2 table from assemblies - SLURM job generator')
    parser.add_argument('--walltime', type=str, required=False, help='wall time (default 0-0:30 = 30 min)', default='0-0:30')
    parser.add_argument('--memory', type=str, required=False, help='mem (default 4096 = 4gb)', default='4096')
    parser.add_argument('--rundir', type=str, required=False, help='directory to run in (default current dir)')
    parser.add_argument('--script', type=str, required=False, help="path to srst2_table_from_assemblies.py, if not in current directory or this script's directory")
    parser.add_argument('--assemblies', nargs='+', type=str, required=True, help='Fasta file/s for assembled contigs')
    parser.add_argument('--gene_db', type=str, required=True, help='Fasta file for gene databases')
    parser.add_argument('--output', type=str, required=True, help='Identifier for outputs (will be combined with assembly identifiers)')
    parser.add_argument('--min_coverage', type=float, required=False, help='Minimum %%coverage cutoff for gene reporting (default 90)')
    parser.add_argument('--max_divergence', type=float, required=False, help='Maximum %%divergence cutoff for gene reporting (default 10)')
    parser.add_argument('--report_new_consensus', type=str, required=False, help='When matching alleles are not found, report the found alleles in this file (will be combined with assembly identifiers)')
    parser.add_argument('--report_all_consensus', type=str, required=False, help='Report all found alleles in this file (will be combined with assembly identifiers)')
    parser.add_argument('--algorithm', action="store", help="blast algorithm (blastn)")
    parser.add_argument('--mlst', action='store_true', required=False, help="Turn it on to find MLST genes")
    return parser.parse_args()


if __name__ == '__main__':
    main()
