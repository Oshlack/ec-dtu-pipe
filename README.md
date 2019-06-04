EC DTU pipeline
===============

This repository contains a pipeline to run DTU testing in four different ways:

* using salmon equivalence classes
* using salmon transcript quantifications
* using dexseq's exon counts
* using featureCount's exon counts

This pipeline implements the commands used to generate the quantification data for our paper: [Fast and accurate differential transcript usage by testing equivalence class counts](https://doi.org/10.1101/501106).

Setup
-----

Install the following software:

* [bpipe](https://github.com/ssadedin/bpipe/releases)
* [Anaconda](https://www.anaconda.com/distribution/#download-section) (if running exon-based analysis, you will require python 2.7+ to run DEXSeq's counting, otherwise either version is compatible)
  * (equivalently, [Python](https://www.python.org/downloads/) with [Numpy](http://www.numpy.org/) and [Pandas](https://pandas.pydata.org/) can be installed)
* [salmon](https://github.com/COMBINE-lab/salmon) (required for EC and transcript-based analysis; you will need to compile the develop branch in order to use the `--skipQuant` flag)
* [subread](https://sourceforge.net/projects/subread/files/subread-1.6.3/) (required for exon-based analysis using featureCounts)
* [DEXSeq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html) (required for exon-based analysis using DEXSeq's exon counting; download the source package somewhere to your machine)
* [STAR](https://github.com/alexdobin/STAR) (required for exon-based analysis)

Install the following R packages:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DEXSeq")
BiocManager::install("DRIMSeq")
BiocManager::install("tximport")
BiocManager::install("biomaRt")
install.packages('data.table')
install.packages('dplyr')
```

Clone the pipeline repository to your machine:

```
git clone git@github.com:Oshlack/ec-dtu-pipe.git
```

If you reqiure exon-based analysses, clone the following [respository](https://github.com/markrobinsonuzh/diff_splice_paper.git) to your machine:

```
git clone https://github.com/markrobinsonuzh/diff_splice_paper.git
```

Ensure that the `params.txt` file is set up with the correct paths and options. Refer to `sample_params.txt` for an example setup using the drosophila data from Soneson et al[1].

References
----------

Download the following references for your organism of choice:

* genome fasta
* transcriptome fasta
* transcriptome GTF (containing exon locations)

The pipeline needs to match transcript IDs to gene IDs in order to run _DEXSeq_. There are two ways to do this. Either make sure [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) is installed and then add the following to your `params.txt` file:

```
-p bmart_dset=<dataset>
```

For example, your dataset may be `hsapiens_gene_ensembl` if you are using human data.

Alternatively, you can create a reference lookup file containing gene, transcript, exon number, gene name and exon ID (space-separated); for example:

```
ENSG00000198888 ENST00000361390 1 MT-ND1 ENSE00001435714
ENSG00000198763 ENST00000361453 1 MT-ND2 ENSE00001435686
ENSG00000198804 ENST00000361624 1 MT-CO1 ENSE00001435647
ENSG00000198712 ENST00000361739 1 MT-CO2 ENSE00001435613
```

Ensure that this file referenced in your `params.txt` file, e.g.:

```
-p tx_lookup=<path to lookup file>
```

In the case where your transcript IDs in your fasta reference file do not contain your desired transcript reference IDs, you will need to create a lookup file. This may be required as salmon will take the transcript ID to be the first `name` after the `>`. For example, the `ENST00000417324` transcript will be referred to as `41361` if your fasta format is as follows:

```
>41361 ENST00000417324 1- 34554-35174,35277-35481,35721-36081
AGCAGCAGGAGTGTTTTAATTAAAAGAAGGCAGTTGCTGTAACCAACTATAAACAAATAA
```

This file isn't necessary if your transcript IDs come right after the `>`.

If a lookup file is required, it must have the IDs from the fasta in the first column, and the transcript IDs in the second column. The transcript IDs must match whatever reference you are using in biomaRt or in your transcript reference file (ensembl, for example). For example:

```
41361 ENST00000417324
41362 ENST00000461467
41363 ENST00000594647
41364 ENST00000335137
```

Fastq format
------------

It is important that your fastq files follow your fastq mask. For example, if our files are in the format:

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
