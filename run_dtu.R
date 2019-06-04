# DE/DTU
library(DEXSeq)
library(DRIMSeq)
library(tximport)

# utils
library(data.table)
library(dplyr)
library(readr)

# load helper methods
args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("--file=", args, value = TRUE)
data_helper <- gsub("--file=(.*)run_dtu.R","\\1R/load_data.R", file.arg)
dtu_helper <- gsub("--file=(.*)run_dtu.R","\\1R/dtu.R", file.arg)
source(data_helper)
source(dtu_helper)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5 | args[1] == '-h') {
    print('Usage: Rscript run_dtu.R <feature> <data> <group> <sample_regex> <outfile> [--ref <tx_ref>] or [--biomart <dataset>] and optionally: [--lookup <tx_lookup_file>], [--cores <N>]')
}

feature <- args[1]
dat <- args[2]
group <- args[3]
sample_regex <- args[4]
outfile <- args[5]

# initialise and set optional vars
tx_ref_file <- ''; tx_lookup_file <- ''; dataset <- ''
if (length(args) > 5) {
    if ('--biomart' %in% args) {
        dataset <- args[which(args == '--biomart') + 1]
    } else {
        tx_ref_file <- args[which(args == '--ref') + 1]
    }
    if ('--lookup' %in% args) {
        tx_lookup_file <- args[which(args == '--lookup') + 1]
    }
}
cores <- 1
if ('--cores' %in% args) {
    cores <- as.numeric(args[which(args == '--cores') + 1])
}

if (feature == 'ec') {
    feat_data <- load_ec_data(dat, tx_lookup=tx_lookup_file, reference=tx_ref_file, bmart_dset=dataset)
} else if (feature == 'tx') {
    feat_data <- load_tx_data(dat, tx_lookup=tx_lookup_file, reference=tx_ref_file, bmart_dset=dataset)
} else if (feature %in% c('ex', 'ex_fc')) {
    feat_data <- load_ex_data(dat, sample_regex)
    feature = 'ex'
}

group <- as.character(strsplit(group, ',')[[1]])
samples <- colnames(feat_data)[grep(sample_regex, colnames(feat_data))]
samples <- as.character(sapply(samples, function(x){strsplit(x, 'Aligned')[[1]][1]}))
group <- as.numeric(samples %in% group)

results <- run_diffsplice(feat_data, group, sample_regex, feature=feature, cores=cores)
save(list = c('results', 'feat_data'), file = outfile)
