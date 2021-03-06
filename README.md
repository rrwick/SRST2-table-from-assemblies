# SRST2 table from assemblies

This is a tool for conducting a gene screen on assemblies.  It produces a table which mimics the format of [SRST2](https://github.com/katholt/srst2)'s compiled results.

It can be useful when there are some samples for which you have reads and other samples for which you only have assemblies.  You can use SRST2 to conduct a gene screen on the reads and this script to conduct a gene screen on the assemblies.  The results can then be combined using SRST2.  Note that if you have both reads and assemblies, it is preferrable to use SRST2 on the reads instead of this script, as SRST2 has better sensitivity.

It uses [BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK279690/) to conduct the gene searches and requires BLAST+ to be installed.

## Usage:
```
srst2_table_from_assemblies.py [-h]
                               --assemblies ASSEMBLIES [ASSEMBLIES ...]
                               --gene_db GENE_DB
                               --output OUTPUT
                               [--min_coverage MIN_COVERAGE]
                               [--max_divergence MAX_DIVERGENCE]
                               [--report_new_consensus REPORT_NEW_CONSENSUS]
                               [--report_all_consensus REPORT_ALL_CONSENSUS]
                               [--algorithm ALGORITHM]
```

* `--assemblies`: FASTA files of all assemblies to screen.  The sample name will be taken from the assembly filename (without the `.fasta` or `.fa` extension).  If a BLAST database does not exist for each assembly (`.nhr`, `.nin` and `.nsq` files), then it will be made using `makeblastdb`.  Since doing so creates new files, you will need write permission to the directory of the assembly.
* `--gene_db`: a gene database to search for in [SRST2 format](https://github.com/katholt/srst2#generating-srst2-compatible-clustered-database-from-raw-sequences).
* `--output`: The output prefix for the table of results.  The table's full name will be [prefix]__genes__[gene_db]__results.txt
* `--min_coverage`: the minimum allowed BLAST hit percent coverage for an allele (default = 90).
* `--max_divergence`: the maximum allowed percent divergence between a BLAST hit and an allele sequence (default = 10).
* `--report_new_consensus`: When matching alleles are not found, report the found alleles in this file (default: do not save allele sequences to file).
* `--report_all_consensus`: Report all found alleles in this file (default: do not save allele sequences to file).
* `--algorithm`: which BLAST+ algorithm to use (`blastn`, `blastn-short`, `megablast` or `dc-megablast`, default = `blastn`).

### Example command

`srst2_table_from_assemblies.py --assemblies *.fasta --gene_db gene_db.mfasta --output test --report_new_consensus new_alleles.fasta`

This command will:
* Screen every one of the assemblies (all `*.fasta` files) for each of the genes in `resistance_genes.mfasta` using `blastn`.
* Save a table of results to `test__genes__gene_db__results.txt`.
* Save a FASTA file of any new alleles to `new_alleles.fasta`.

### Example query gene database

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

For this script, the important parts are the gene cluster name and the allele name, the second and third pieces of the sequence header delimited with `__`.

Ideally, every allele in your database has a unique name.  If this is the case, this script's output table will use only the allele names.  If it is not the case (i.e. there is at least one duplicate allele name in the database), then the output table will use a combination of the allele name and the allele unique ID (the last piece of the sequence header).  This is the same behaviour as SRST2.

### Example output table

Sample | abcA | abcB
--- | --- | --- |
sample1 | abcA_1 | abcB_1 |
sample2 | abcA_1 | abcB_2* |
sample3 | abcA_2 | - |

As is the case for SRST2:
* Imperfect matches (containing at least one mismatch or indel) are indicated with '*' after the allele name.
* Absent genes are indicated with '-'.

When no perfect allele match is present but there are multiple possible alleles (which satisfy the coverage and divergence thresholds), this script will choose the one with the highest BLAST bit score.  Using the above table as an example, if sample2 matched abcB_1 with a bit score of 123.4 and abcB_2 with a bit score of 134.5, this script will choose abcB_2 and include a '*' to indicate the match was not exact.

## SLURM queueing system

Another script, `srst2_table_from_assemblies_slurm.py`, is included to generate SLURM jobs to run many gene screens in parallel.  This will generate a separate output table for each assembly, so you may want to then [compile them together using SRST2](https://github.com/katholt/srst2#running-lots-of-jobs-and-compiling-results).

## License

GNU General Public License, version 3
