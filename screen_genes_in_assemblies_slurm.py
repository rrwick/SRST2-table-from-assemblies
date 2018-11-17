#!/usr/bin/env python
"""
SRST2 results from assemblies

This is a tool to screen for genes in a collection of assemblies and output
the results in a table which mimics those produced by SRST2.

Python versions 2.7 and 3 compatible.
Previous name: srst2_table_from_assemblies_slurm.py

Copyright (C) 2015-2017 Ryan Wick <rrwick@gmail.com>, Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
First edition: 9-10 Sep 2017; the latest edition: 17 Nov 2018
"""

from __future__ import print_function
import sys
sys.dont_write_bytecode = True  # Do not write .pyc files on the import of source modules.
import argparse
import os
import time
from screen_genes_in_assemblies import check_file_exists, check_algorithm, rchop


def get_arguments():
    parser = argparse.ArgumentParser(description = "Screen genes in assemblies - a SLURM job generator")
    
    # SLURM arguments
    parser.add_argument("--script", type = str, required = False, default = "./screen_genes_in_assemblies.py", help = "path to screen_genes_in_assemblies.py, if not in the current directory or this script's directory")
    parser.add_argument("--walltime", type = str, required = False, default = "0-0:30:0", help = "wall time for each job bundle (default: 0-0:30:0 = 30 min)")
    parser.add_argument("--memory", type = str, required = False, default = "512", help = "memory assigned to every job (default: 512 MB)")
    parser.add_argument("--account", type = str, required = False, default = "", help = "SLURM account name (default: null)")
    parser.add_argument("--partition", type = str, required = False, default = "main", help = "name of the queue partition (default: main)")
    parser.add_argument("--bundle_name_prefix", type = str, required = False, default = "genotyping", help = "prefix for the name of every bundle of SLURM jobs (default: genotyping)")
    parser.add_argument("--bundle_size", type = int, required = False, default = 16, help = "number of individual jobs per bundle (default: 16)")
    parser.add_argument("--outdir", type = str, required = False, default = os.getcwd(), help = "output directory for results (default: current dir)")
    parser.add_argument("--dont_run", action = "store_true", required = False, help = "Flag it to only print SLURM scripts without submission")
    
    # script arguments
    parser.add_argument("--assemblies", nargs = "+", type = str, required = True, help = "Fasta file/s for assembled contigs")
    parser.add_argument("--gene_db", type = str, required = True, help = "Fasta file for gene databases")
    parser.add_argument("--algorithm", type = str, required = False, default = "megablast", help = "blast algorithm (megablast)")
    parser.add_argument("--prefix", type = str, required = False, default = "BLAST", help = "Output prefix for the table of results")
    parser.add_argument("--suffix", type = str, required = False, default = ".fasta", help = "Characters to be chopped off from the end of every assembly name in order to get a sample name")
    parser.add_argument("--other_args", type = str, required = False, help = "A single string consisting of other arguments to be passed directly to screen_genes_in_assemblies.py")
    parser.add_argument("--serial", action = "store_true", required = False, help = "Run jobs within the same bunch in a serial manner.")
    
    # environment settings
    parser.add_argument("--blast", type = str, required = False, default = "", help = "Module name of BLAST")
    parser.add_argument("--python", type = str, required = False, default = "", help = "Module name of Python")  # load it after BLAST to avoid changing of Python versions in some systems
    parser.add_argument("--other_modules", type = str, required = False, default = "", help = "Comma-delimited list of modules")
    
    return parser.parse_args()


def main():
    args = get_arguments()
    run = not args.dont_run

    # Environmental settings ###############
    check_file_exists(args.gene_db)
    args.gene_db = os.path.abspath(args.gene_db)
    if args.algorithm:
        check_algorithm(args.algorithm)
    if args.script:
        if args.script.endswith("screen_genes_in_assemblies.py"):
            script_path = args.script
        else:
            script_path = os.path.join(args.script, "screen_genes_in_assemblies.py")
        check_file_exists(script_path)
    else:
        script_path_cwd = os.path.join(os.getcwd(), "screen_genes_in_assemblies.py")
        if os.path.isfile(script_path_cwd):
            script_path = script_path_cwd
        else:
            script_path_file_directory = os.path.join(os.path.dirname(os.path.abspath(__file__)), "screen_genes_in_assemblies.py")
            if os.path.isfile(script_path_file_directory):
                script_path = script_path_file_directory
            else:
                sys.exit("Error: could not find the script screen_genes_in_assemblies.py")

    # Generate and submit a SLURM script for each assembly ###############
    """
    We use job bundles to solve the problem of losing outputs when a large number of individual short jobs
    are submitted to the SLURM system. The Melbourne Bioinformatics (https://www.melbournebioinformatics.org
    .au/documentation/running_jobs/slurm_x86/) suggests to use job bundles to address this problem. In this
    method, the total amount of memory required for each bundle equals args.memory * args.bundle_size with a
    unit of MB.
    """
    bundle_count = 0
    job_count = 0
    job_left = len(args.assemblies)
    tasks = ""
    
    if args.serial:
        job_line_terminal = "\n"
    else:
        job_line_terminal = " &\n"  # run jobs of the same bundle in parallel
    
    for assembly_filename in args.assemblies:
        assembly_name = rchop(os.path.basename(assembly_filename), args.suffix)
        
        # make a the command line for a single task
        tasks += "srun --nodes=1 --ntasks=1 --cpus-per-task=1 python " + script_path
        tasks += " --assemblies " + assembly_filename
        tasks += " --gene_db " + args.gene_db
        tasks += " --outdir " + args.outdir
        tasks += " --algorithm " + args.algorithm
        """
        The following concatenation of assembly_name and args.prefix is essential for the "srst2 --prev_output" command
        to run properly, otherwise, it reports an error that "Couldn't decide what to do with file results".
        """
        tasks += " --prefix " + assembly_name + "_" + args.prefix  # to comply with the SRST2 convention: sample_prefix__[database name]...
        if args.suffix != "":
            tasks += " --suffix " + args.suffix
        if args.other_args != "":
            tasks += " " + args.other_args + job_line_terminal  # send this task to the background
        else:
            tasks += job_line_terminal
        job_count += 1
        job_left -= 1
        
        # submit a bundle when there are args.bundle_size jobs
        if job_count == args.bundle_size or job_left == 0:
            bundle_count += 1
            bundle_cmd = "#!/bin/bash"
            if args.account != "":
                bundle_cmd += "\n#SBATCH --account=" + args.account
            bundle_cmd += "\n#SBATCH --partition=" + args.partition
            bundle_cmd += "\n#SBATCH --job-name=" + args.bundle_name_prefix + "_" + str(bundle_count)
            if args.serial:
                bundle_cmd += "\n#SBATCH --ntasks=1"
            else:
                bundle_cmd += "\n#SBATCH --ntasks=" + str(job_count)
            bundle_cmd += "\n#SBATCH --cpus-per-task=1"
            bundle_cmd += "\n#SBATCH --mem-per-cpu=" + args.memory  # memory per task (CPU) when there has to be one processor (CPU) per task
            bundle_cmd += "\n#SBATCH --time=" + args.walltime
            if args.blast != "":
                bundle_cmd += "\nmodule load " + args.blast  # change to your own module names
            if args.python != "":
                bundle_cmd += "\nmodule load " + args.python
            if args.other_modules != "":
                add_modules = args.other_modules.split(",")
                for m in add_modules:
                    bundle_cmd += "\nmodule load " + m + "\n"
            bundle_cmd += tasks + "wait\n"  # terminate this bundle only when all tasks end
            if run:
                os.system("echo '" + bundle_cmd + "' | sbatch")  # submit this bundle of jobs
            print(bundle_cmd)  # put this and the following two commands after the submission to increase the interval between two submissions
            job_count = 0
            tasks = ""
            time.sleep(1)  # futher delay the submission process in order to give the scheduler enough time to process all commands


if __name__ == "__main__":
    main()
