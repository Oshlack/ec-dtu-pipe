EC DTU pipeline
===============

This repository contains a pipeline to run DTU testing in four difference ways:

* using salmon equivalence classes
* using salmon transcript quantifications
* using dexseq's exon counts
* using featureCount's exon counts

This pipeline implements the commands used to generate the quantification data for our paper: [Fast and accurate differential transcript usage by testing equivalence class counts](https://doi.org/10.1101/501106]).

Setup
-----

Install the following software:

* [bpipe](https://github.com/ssadedin/bpipe/releases)
* [salmon](https://github.com/COMBINE-lab/salmon) (you will need to compile the develop branch in order to use the `--skipQuant` flag)
* [subread](https://sourceforge.net/projects/subread/files/subread-1.6.3/)
* [DEXSeq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html) (you will need to install DEXSeq in R via bioconductor, *and* download the source package somewhere to your machine)
* [STAR](https://github.com/alexdobin/STAR)

Additionally, clone the following [respository](https://github.com/markrobinsonuzh/diff_splice_paper.git) to your machine:

```
git clone https://github.com/markrobinsonuzh/diff_splice_paper.git
```

Now download the following references for your organism of choice:

* genome fasta
* transcriptome fasta
* transcriptome GTF (containing exon locations)
* reference lookup containing gene, transcript, exon number, gene name and exon ID (space-separated); for example:

```
ENSG00000198888 ENST00000361390 1 MT-ND1 ENSE00001435714
ENSG00000198763 ENST00000361453 1 MT-ND2 ENSE00001435686
ENSG00000198804 ENST00000361624 1 MT-CO1 ENSE00001435647
ENSG00000198712 ENST00000361739 1 MT-CO2 ENSE00001435613
```

Lastly, you'll need a transcript ID lookup file that looks like this:

```
41361 ENST00000417324
41362 ENST00000461467
41363 ENST00000594647
41364 ENST00000335137
```

This is because salmon will take the transcript ID to be the first `name` after the `>`. For example, the `ENST00000417324` transcript will be referred to as `41361` if your fasta format is as follows:

```
>41361 ENST00000417324 1- 34554-35174,35277-35481,35721-36081
AGCAGCAGGAGTGTTTTAATTAAAAGAAGGCAGTTGCTGTAACCAACTATAAACAAATAA
```

This file shouldn't be necessary if your transcript IDs come right after the `>`, but for the meantime the pipeline requires this file. In the case that your transcript IDs are the same as salmon's IDs, create a file with identical columns, separated by a space:

```
ENST00000417324 ENST00000417324
ENST00000461467 ENST00000461467
ENST00000594647 ENST00000594647
ENST00000335137 ENST00000335137
```

The correct paths and arguments will then need to be set up in the `params.txt` file (see `sample_params.txt` for an example for setup using the drosophila data Soneson[1].

Fastq format
------------

It is important that your fastq files follow your fastq files follow your fastq mask. For example, if our files are in the fornmat:

```
Dm_sample_1_1.fq.gz
Dm_sample_1_2.fq.gz
Dm_sample_2_1.fq.gz
Dm_sample_2_2.fq.gz
```

Our mask would be: `%_sample_%_*.fq.gz`. See the [bpipe documentation](http://docs.bpipe.org/Overview/Introduction/) for more details.

Your samples must all share a common pattern. For example, all the files above contain `Dm`. This will be defined as `sample_regex` in the params file. The pattern does not need to be at the start. Additionally, to perform DTU, you will have to list all samples belonging to your second experimental group under `-p group=` in the `params.txt` file (this will correspond to the names of the samples in the `%` part of the file mask.

Running the pipeline
--------------------

The code in this repository allows you to run four different pipelines.

To perform DTU analysis using equivalence classes, run the pipeline as follows:

```
bpipe run @params.txt ec-dtu-pipe.groovy <fastq files>
```

To perform DTU analysis using transcript quantifications, run the pipeline as follows:

```
bpipe run @params.txt tx-dtu-pipe.groovy <fastq files>
```

To perform DTU analysis using exon counts derived from DEXSeq's count, run the pipeline as follows:

```
bpipe run @params.txt ex-dtu-pipe.groovy <fastq files>
```

To perform DTU analysis using exon counts derived from featureCounts, run the pipeline as follows:

```
bpipe run @params.txt ex-fc-dtu-pipe.groovy <fastq files>
```

Results
-------

The pipeline will save the results from DTU in an RData file containing `feat_data` and `results` objects. `feat_data` contains the input data to DEXseq, and `results` contains the DEXSeq output. Access the results as follows:

* `results[['dexseq_object']]` is the DEXSeq object
* `results[['dexseq_results']]` is the DEXSeq results object
* `results[['gene_FDR']]` is the results from DEXSeq's `perGeneQValue` method

References
----------
[1] Soneson, C., Matthes, K. L., Nowicka, M., Law, C. W., & Robinson, M. D. (2016). Isoform prefiltering improves performance of count-based methods for analysis of differential transcript usage. Genome Biology, 17(1), 1–15. https://doi.org/10.1186/s13059-015-0862-3
