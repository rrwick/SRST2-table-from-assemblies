#!/bin/bash
# Pool FASTA files of each strain into a single multi-FASTA file, which simplifies the gene screen using
# the script screen_genes_in_assemblies.py.
#
# Usage:
#   pool_seqs.sh --in_dir="./fasta" --out_dir="./fasta_pooled/" --suffix "__*.fasta" strain1 strain2 ...
#
# Copyright (C) 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
# First edition: 7 Apr 2018; the latest edition: 8 Apr 2018

# Defaults
strains=()
in_dir="./fasta"
out_dir="./fasta_pooled"
suffix="__*.fasta"

# Read arguments
for i in "$@"; do  # loop through every argument (e.g. "--gbk=*.gbk" is treated as a single word by the $@ operator)
    case $i in
		--in_dir=*)  # where individual FASTA files are stored
		in_dir="${i#*=}"  # script that extracts features
		;;
        --out_dir=*)  # directory for multi-FASTA files: [strain name].fasta
        out_dir="${i#*=}"  # output directory without the forward slash
        ;;
        --suffix=*)  # input FASTA files should be named in the format: [strain name][suffix]
        suffix="${i#*=}"
        ;;
        *)  # a vector of space-delimited strain names
        strains+=(${i})  # Other arguments are treated as GenBank filenames.
        ;;
    esac
done

# Pool sequences
for s in "${strains[@]}"; do
    echo "Pooling sequences of the strain ${s}."
    cat ${in_dir}/${s}${suffix} > ${out_dir}/${s}.fasta
done
