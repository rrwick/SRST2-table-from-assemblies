# SRST2 table from assemblies

This is a tool for conducting a gene screen on assemblies.  It produces a table which mimics the format of [SRST2](https://github.com/katholt/srst2)'s compiled results.

It can be useful when there are some samples for which you have reads and other samples for which you only have assemblies.  You can then use SRST2 to conduct a gene screen on the reads and this script to conduct a gene screen on the assemblies.  The results can then be combined using SRST2.  Note that if you have both reads and assemblies, it is preferrable to use SRST2 on the reads instead of this script (better sensitivity).

It uses [BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK279690/) to conduct the gene searches and requires BLAST+ to be installed.

## Usage:
```
srst2_table_from_assemblies.py [-h] --assemblies ASSEMBLIES [ASSEMBLIES ...]
                               --gene_db GENE_DB
                               --output OUTPUT
                               [--min_coverage MIN_COVERAGE]
                               [--max_divergence MAX_DIVERGENCE]
                               [--report_new_consensus]
                               [--report_all_consensus]
                               [--algorithm ALGORITHM]
```

* `--assemblies`: FASTA files of all assemblies to screen.  The sample name will be taken from the assembly filename (without the `.fasta` or `.fa` extension).  If a BLAST database does not exist for each assembly (`.nhr`, `.nin` and `.nsq` files), then it will be made using `makeblastdb`.
* `--gene_db`: a gene database to search for in [SRST2 format](https://github.com/katholt/srst2#generating-srst2-compatible-clustered-database-from-raw-sequences).
* `--output`: the filename for the resulting tab-delimited table.
* `--min_coverage`: the minimum allowed BLAST hit percent coverage for an allele (default = 90).
* `--max_divergence`: the maximum allowed percent divergence between a BLAST hit and an allele sequence (default = 10).
* `--report_new_consensus`: save all novel alleles (those which do not exactly match any in the gene database) to a FASTA file named `new_consensus_alleles.fasta`.
* `--report_all_consensus`: save all found alleles (both those in the gene database and novel alleles) to a FASTA file named `all_consensus_alleles.fasta`.
* `--algorithm`: which BLAST+ algorithm to use (`blastn`, `blastn-short`, `megablast` or `dc-megablast`, default = `blastn`).

### Example command

`srst2_table_from_assemblies.py --assemblies *.fasta --gene_db resistance_genes.mfasta --output table.txt --report_new_consensus`

This command will:
* Screen every one of the assemblies (all `*.fasta` files) for each of the genes in `resistance_genes.mfasta` using `blastn`.
* Save a table of results to `table.txt`.
* Save a FASTA file of any new alleles to `new_consensus_alleles.fasta`.

### Example query gene database

```
>0__abcA__abcA_1__0
ATGGACTTTTCCCGCTTTTATATCGACAGGCCGATCTTC
>0__abcA__abcA_2__1
ATGGACTTTTCCCGCTTTTTTATCGACAGGCCGATTTTC
>1__abcB__abcB_1__2
ATGAGAAATAAAGGAATCGATCAATTTTGTGTGATTGCA
>1__abcB__abcB_2__3
ATGAAAAATAAAGGAATCGATCAATTTCGTGTGATTGCA
```

For this script, the important parts are the gene cluster name and the allele name, the second and third pieces of the sequence header delimited with `__`.

### Example output table

Sample | abcA | abcB
--- | --- | --- |
sample1 | abcA_1 | abcB_1 |
sample2 | abcA_1 | abcB_2* |
sample3 | abcA_2 | - |

As is the case for SRST2:
* Imperfect matches (containing at least one mismatch or indel) are indicated with '*' after the allele name.
* Absent genes are indicated with '-'.

## License

GNU General Public License, version 3
