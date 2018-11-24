# Gene screen for genome assemblies

[![DOI](https://zenodo.org/badge/65195290.svg)](https://zenodo.org/badge/latestdoi/65195290)

This is a tool for conducting a gene screen on assemblies. It produces a table in the format of [SRST2](https://github.com/katholt/srst2)'s compiled results. It is particularly useful when you want to screen for genes in samples where some of them have reads available while the others only have assemblies. For reads, you can use SRST2 to perform the gene screen; and for assemblies, you may use this script to produce results that can be compiled with those obtained from reads. Nonetheless, note that if you have both reads and assemblies for the same sample, it is preferable to use SRST2 on the reads instead of using this script on assemblies, as SRST2 has better sensitivity (variants may get lost for various reasons during the assembly process).

Although this tool is designed for screening genes in haploid organisms, it is in principle applicable to other organisms as well.  

Dependency: this tool is Python 2 and 3 compatible. It requires a local [BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK279690/) installation to conduct nucleotide-level gene searches.

Citation: Ryan R. Wick, Yu Wan, _gene\_screen\_in\_assemblies_, https://github.com/wanyuac/screen_genes_in_assemblies.git, doi: 10.5281/zenodo.893164.

## Arguments and options

```
python screener.py -h
       --assemblies ASSEMBLIES
       [ASSEMBLIES ...] --gene_db GENE_DB
       [--prefix PREFIX] [--suffix SUFFIX]
       [--outdir OUTDIR]
       [--min_coverage MIN_COVERAGE]
       [--max_divergence MAX_DIVERGENCE]
       [--report_new_consensus]
       [--report_all_consensus]
       [--algorithm ALGORITHM] [--mlst]
       [--incl_alt]
       [--max_overlapping_nt MAX_OVERLAPPING_NT]
       [--del_blast]
```

* `--assemblies`: FASTA files of all assemblies to screen. The sample name will be taken from the assembly filename (without the `.fasta` or `.fa` extension). If a BLAST database does not exist for each assembly (`.nhr`, `.nin` and `.nsq` files), then it will be made using `makeblastdb`. Since doing so creates new files, you will need write permission to the directory of the assembly.
* `--gene_db`: a gene database to search for in the [SRST2-compatible format](https://github.com/katholt/srst2#generating-srst2-compatible-clustered-database-from-raw-sequences).
* `--prefix`: Output prefix for the table of results. It works in the same way as SRST2.
* `--suffix`: Characters to be chopped off from the end of every assembly name in order to get a sample name. For example, *strain* is extracted from the file name *strain_spades.fasta* given `--suffix '_spades.fasta'`.  
* `--outdir`: Output directory for the table of results.
* `--min_coverage`: the minimum allowed BLAST hit percent coverage for an allele (default = 90).
* `--max_divergence`: the maximum allowed percent divergence between a BLAST hit and an allele sequence (default = 10).
* `--report_new_consensus`: When matching alleles are not found, report the found alleles in this file (default: do not save allele sequences to file).
* `--report_all_consensus`: Report all found alleles in this file (default: do not save allele sequences to file).
* `--algorithm`: which BLAST+ algorithm to use for DNA sequence alignment (`blastn`, `blastn-short`, `megablast` or `dc-megablast`, default = `megablast`).
* `--mlst`: Turn it on to find MLST (multi-locus sequence typing) genes.
* `--incl_alt`: Flag it to include all putative alternative calls for each allele.
* `--max_overlapping_nt`: Maximal number of overlapping nucleotides allowed to treat two hits as physically separate.
* `--del_blast`: Flag it to avoid saving raw BLAST outputs.

## Parallel gene screen through the SLURM queueing system

Another script, screener\_slurm.py, is included to generate a series of SLURM jobs to efficiently run many gene screens in parallel, although the script screener.py is also able to handle multiple FASTA files (but in a series manner). Since this script produces a separate output table for each assembly, so a user may want to [compile them together using SRST2](https://github.com/katholt/srst2#running-lots-of-jobs-and-compiling-results) afterwards.

## Outputs
Assuming assembly files are named in the format of \[sample name\]\[suffix\].fasta, and this tool is run with the correct `--suffix` specification.
 
### screener.py
There are four files generated for all samples when options `--report_new_consensus`, `--report_all_consensus` and `--del_blast` are turned on.

* [prefix]\_\_genes\_\_[gene\_db]\_\_results.txt: gene profiles for all samples. It is equivalent to the compiled gene profiles produced by the SRST2 command `python srst2.py --prev_output`.
* [prefix].all\_consensus\_alleles.fasta: consensus sequences, namely, the alignments on sample contigs, of all samples.
* [prefix].new\_consensus\_alleles.fasta: consensus sequences of all valid imperfect matches (namely, variants) from all samples.
* [prefix]\_\_blast\_\_[gene\_db]\_\_results.txt: raw BLAST outputs with a header line attached for each sample.

### screener\_slurm.py
This script makes output files for each sample, which is different from the behaviour of screener.py.

* [sample]\_[prefix]\_\_genes\_\_[gene\_db]\_\_results.txt: a gene profile is generated for each sample.
* [sample]\_[prefix].all\_consensus\_alleles.fasta: consensus sequences of the current sample.
* [sample]\_[prefix].new\_consensus\_alleles.fasta: consensus sequences of variants from the current sample.
* [sample name]\_\_[sample name]\_[prefix]\_\_blast\_\_[gene\_db]\_\_results.txt: raw BLAST outputs with a header line attached for the current sample.

## Examples

### 1. Screening for the best allele call per gene
This example strictly follows the output format of SRST2 to make a single (but the best) allele call for each sequence cluster (representing a gene) in the reference gene database. In this example, the option `--incl_alt` is left off.  

`screener.py --assemblies *.fasta --gene_db gene_db.fasta --algorithm blastn --output test --report_new_consensus`

This command (1) screens every one of the assemblies (all `*.fasta` files) for each of the genes in `resistance_genes.mfasta` using `blastn`; (2) saves a table of results to test\_\_genes\_\_gene\_db\_\_results.txt; and (3) saves a FASTA file of any new alleles into new\_alleles.fasta.

**An example of reference gene databases**  
Sequences of this database are aligned against every assembly as queries.  

```
>0__abcA__abcA_1__0
TCGCAGGGCGAGCGGCGCGTCTCACGGAATGACCATGTCCTGCATCATAAATTAACGTAA
>0__abcA__abcA_2__1
TCGCAGGGCGAGCGCCGCGTCTCACGGAATGACCATGTCCTGCATGATAAATTAACGTAA
>1__abcB__abcB_1__2
TAATAGTGATGGGTATTGAGGGCTCCCCTTGAAGCCTCGCAGAAAGCAGATCAATTTCAA
>1__abcB__abcB_2__3
TAATAGTGATGGGTTTTGAGGGCTCGCCTTGAAGCCTCGCAGATAGCAGATCAATTTCAA
```

For this tool, the important parts are the gene cluster name and the allele name, the second and third pieces of the sequence header delimited with double underlines "\_\_". Ideally, every allele in your database has a unique name.  If this is the case, this script's output table will use only the allele names.  If it is not the case (i.e. there is at least one duplicate allele name in the database), then the output table will use a combination of the allele name and the allele unique ID (the last piece of the sequence header).  This is the same behaviour as SRST2.

**Example output table**

Sample | abcA | abcB
--- | --- | --- |
sample1 | abcA_1 | abcB_1 |
sample2 | abcA_1 | abcB_2* |
sample3 | abcA_2 | - |

As is the case for SRST2:

* Imperfect matches (containing at least one mismatch or indel) are indicated with '*' after the allele name.
* Absent genes are indicated with '-'.

When no perfect allele match is present but there are multiple possible alleles (which satisfy the coverage and divergence thresholds), this script will choose the one with the highest BLAST bit score. Using the above table as an example, if sample2 matched abcB\_1 with a bit score of 123.4 and abcB\_2 with a bit score of 134.5, this script will choose abcB\_2 and include a '*' to indicate the match was not exact.

### 2. MLST
This is a special case of calling the best alleles for a panel of seven genes given a reference DNA database. A script format\_mlst\_db.py was developed to convert headers (also known as the sequence definition lines) of reference sequences into an SRST2-compatible format.

An example command is:
```
python screener_slurm.py --walltime '0-1:0:0' --outdir mlst --script screener.py --assemblies assemblies/*.fasta --gene_db mlst_db.fasta --prefix test --suffix ".fasta" --report_all_consensus --mlst > mlst.log
```

### 3. Screening for all valid alleles per gene
This is an extension of the SRST2 output format, where multiple allele calls are enabled for each cluster. This is particularly useful to identify all copy-number variants of each gene in assemblies (more powerful in finished-grade genomes). The option `--incl_alt` and the argument `--max_overlapping_nt` work together to control this utility. To be more specific, allele calls of the same cluster in the same sample are not overlapping by more than `--max_overlapping_nt` bp (accordingly they are considered as physically separated).

**Example command lines**  
A single-job version for two assemblies

```
python screener.py --gene_db ARGannot_r2.fasta --prefix demo1 --suffix '_spades.fasta' --outdir genes --report_new_consensus --report_all_consensus --algorithm megablast --incl_alt --max_overlapping_nt 0 --assemblies strain1_spades.fasta strain2_spades.fasta
```

A parallel version for a large number of assemblies

```
python screener_slurm.py --script screener.py --algorithm megablast --walltime "0-0:30:0" --memory 1024 --partition project1 --bundle_name_prefix test --bundle_size 16 --outdir genes --assemblies data/assemblies/*.fasta --prefix demo2 --suffix '_spades.fasta' --gene_db ARGannot_r2.fasta --other_args "--incl_alt --max_overlapping_nt 0 --report_all_consensus" > gene_screen.log
```

**Example output table**

Sample | abcA | abcB
--- | --- | --- |
sample1 | abcA\_1,abcA\_1\*,abcA\_1\* | abcB\_1 |
sample2 | abcA\_1 | abcB\_2\*,abc\_2[2] |
sample3 | abcA\_2\*[2] | - |

This table illustrates that there are three alleles of the cluster abcA identified in the first sample, including an exact match to the abcA\_1 sequence and two relevant variants of different sequences. The number within a pair of square brackets denotes the copy number of the exact sequence found in different places (either on different contigs or distinct positions of the same contig. The latter scenario is controlled by the `--max_overlapping_nt` argument) in the same sample. For instance, there are two copies of the same abcA\_2 variants identified in the sample 3. The order of allele calls under each cluster depends on the Python iterator, and consensus sequences of the same sample are printed into a corresponding FASTA file in the same order.

The output \*\_\_gene\_\_\*.txt files can be compiled into a single table in the usual way with SRST2.  

## Licence
[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.en.html)