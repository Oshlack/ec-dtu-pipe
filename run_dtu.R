# plotting
library(ggplot2)
library(gridExtra)
# library(UpSetR)
library(VennDiagram)
library(RColorBrewer)

# DE/DTU
library(edgeR)
library(DEXSeq)
library(DRIMSeq)
library(tximport)

# utils
library(data.table)
library(dplyr)
library(readr)

#sources
source('R/load_data.R')
source('R/dtu.R')
source('R/util.R')
source('R/plotting.R')

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4 | args[1] == '-h') {
    print('Usage: Rscript run_dtu.R <feature> <data> <group> <sample_regex> <tx_ref> <tx_lookup>')
}

feature <- args[1]
dat <- args[2]
group <- args[3]
sample_regex <- args[4]
if (length(args) > 4) {
    tx_ref_file <- args[5]
    tx_lookup_file <- args[6]
}

group <- as.numeric(strsplit(group, ',')[[1]])

if (feature == 'ec') {
    feat_data <- load_ec_data(dat, tx_lookup_file, tx_ref_file)
} else if (feature == 'tx') {
    feat_data <- load_tx_data(dat, tx_lookup_file, tx_ref_file)
} else if (feature == 'ex') {
    feat_data <- load_ex_data(dat, sample_regex)
}

results <- run_diffsplice(feat_data, group, sample_regex, feature = feature)
save(list = c('results', 'feat_data'), file='results.RData')
