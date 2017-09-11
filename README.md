# Gene screen for haploid genome assemblies

This is a tool for conducting a gene screen on assemblies. It produces a table in the format of [SRST2](https://github.com/katholt/srst2)'s compiled results.

It is particularly useful when you want to screen for genes in samples where some of them have reads available while the others only have assemblies. For reads, you can use SRST2 to perform the gene screen; and for assemblies, you may use this script to produce results that can be compiled with those obtained from reads. Nonetheless, note that if you have both reads and assemblies for the same sample, it is preferable to use SRST2 on the reads instead of using this script on assemblies, as SRST2 has better sensitivity (variants may get lost for various reasons during the assembly process).

This tool is Python 2 and 3 compatible. It uses [BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK279690/) to conduct the nucleotide-level gene searches and hence requires BLAST+ to be installed.

## Arguments and options  

```
python srst2_table_from_assemblies.py -h
usage: srst2_table_from_assemblies.py [-h] --assemblies ASSEMBLIES
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

  --min_coverage MIN_COVERAGE
                        Minimum %coverage cutoff for gene reporting (default
                        90)
  --max_divergence MAX_DIVERGENCE
                        Maximum %divergence cutoff for gene reporting (default
                        10)
  --report_new_consensus
                        Configure it to save consensus sequences of variants
  --report_all_consensus
                        Configure it to save all consensus sequences
  --algorithm ALGORITHM
                        blast algorithm (blastn)
  --mlst                Turn it on to find MLST genes
  --incl_alt            Flag it to include all putative alternative calls for
                        each allele
  --max_overlapping_nt MAX_OVERLAPPING_NT
                        Maximal number of overlapping nucleotides allowed to
                        treat two hits as physically separate
  --del_blast           Flag it to delete the text file of BLAST outputs
```

* `--assemblies`: FASTA files of all assemblies to screen.  The sample name will be taken from the assembly filename (without the `.fasta` or `.fa` extension).  If a BLAST database does not exist for each assembly (`.nhr`, `.nin` and `.nsq` files), then it will be made using `makeblastdb`.  Since doing so creates new files, you will need write permission to the directory of the assembly.
* `--gene_db`: a gene database to search for in [SRST2 format](https://github.com/katholt/srst2#generating-srst2-compatible-clustered-database-from-raw-sequences).
* `--prefix`: Output prefix for the table of results. It works in the same way as SRST2.
* `--suffix`: Characters to be chopped off from the end of every assembly name in order to get a sample name. For example, ''strain'' is extracted from the file name ''strain_spades.fasta'' given `--suffix '_spades.fasta'`.  
* `--outdir`: Output directory for the table of results.
* `--min_coverage`: the minimum allowed BLAST hit percent coverage for an allele (default = 90).
* `--max_divergence`: the maximum allowed percent divergence between a BLAST hit and an allele sequence (default = 10).
* `--report_new_consensus`: When matching alleles are not found, report the found alleles in this file (default: do not save allele sequences to file).
* `--report_all_consensus`: Report all found alleles in this file (default: do not save allele sequences to file).
* `--algorithm`: which BLAST+ algorithm to use (`blastn`, `blastn-short`, `megablast` or `dc-megablast`, default = `blastn`).

## Outputs
The table's full name will be [prefix]__genes__[gene_db]__results.txt

## Examples

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

[Apache license, version 2](http://www.apache.org/licenses/LICENSE-2.0)
