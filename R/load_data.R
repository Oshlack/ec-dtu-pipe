get_ref <- function(reference, bmart_dset) {
    if (reference == '') {
        library(biomaRt)
        mart <- useMart('ENSEMBL_MART_ENSEMBL',
                        dataset=bmart_dset)
        ensg_attr <- c('external_gene_name',
                       'ensembl_gene_id',
                       'ensembl_transcript_id')
        gtx <- getBM(mart=mart, attributes=ensg_attr)
        colnames(gtx) <- c('symbol', 'gene_id', 'ensembl_id')
    } else {
        gtx <- read.delim(gzfile(reference), header=F, sep=' ')
        colnames(gtx) <- c('gene_id', 'ensembl_id', 'exon', 'symbol', 'exon_id')
    }
    return(gtx)
}


load_ec_data <- function(ec_matrix_file, tx_lookup='', reference='', bmart_dset='', filter_multi_ecs=T, salmon=T) {
    is_gzipped <- tail(strsplit(ec_matrix_file, '\\.')[[1]], 1) == 'gz'
    ec_matrix_file <- ifelse(is_gzipped, paste('zcat <', ec_matrix_file), ec_matrix_file)
    df <- fread(ec_matrix_file)
    gtx <- get_ref(reference, bmart_dset)

    if(salmon & tx_lookup != ''){
        lookup <- read.delim(gzfile(tx_lookup), sep=' ', header=F)
        colnames(lookup) <- c('transcript', 'ensembl_id')
        df <- inner_join(df, lookup, by='transcript')
    } else {
        colnames(df)[colnames(df) == 'transcript'] <- 'ensembl_id'
    }

    df <- inner_join(df, gtx, by='ensembl_id')

    if(filter_multi_ecs) {
        multi_ecs <- data.table(df)[, length(unique(gene_id)), keyby=ec_names]
        multi_ecs <- multi_ecs[multi_ecs$V1 > 1]
        multi_ecs <- multi_ecs$ec_names
        df <- df[!df$ec_names %in% multi_ecs,]
    }
    df <- df[,!colnames(df) %in% c('tx_id', 'transcript', 'ensembl_id', 'exon', 'exon_id', 'symbol')]
    df <- distinct(data.table(df))

    # remove suffixes from sample names
    sample_names <- !colnames(df) %in% c('ec_names', 'ensembl_id', 'gene_id', 'symbol')
    colnames(df)[sample_names] <- sapply(colnames(df)[sample_names], function(x){strsplit(x, '_')[[1]][1]})

    return(df)
}

load_tx_data <- function(dir, tx_lookup='', reference='', bmart_dset='', scaling="scaledTPM") {
    gtx <- get_ref(reference, bmart_dset)

    files <- list.files(dir)
    files <- as.character(sapply(files, function(x){paste0(dir, x)}))
    txi <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance=scaling)
    tx_counts <- data.frame(txi$counts)
    samples <- sapply(files, function(x){strsplit(basename(x), '_quant')[[1]][1]})
    colnames(tx_counts) <- as.character(samples)
    tx_counts$Name <- as.character(rownames(tx_counts))

    if(tx_lookup != '') {
        lookup <- read.delim(tx_lookup, sep=' ', header=F)
        colnames(lookup) <- c('Name', 'ensembl_id')
        lookup$Name <- as.character(lookup$Name)
        df <- inner_join(tx_counts, lookup, by='Name')
    } else {
        df <- tx_counts
        colnames(df)[colnames(df)=='Name'] <- 'ensembl_id'
    }
    df <- inner_join(df, gtx, by='ensembl_id')

    df[is.na(df)] <- 0
    df <- distinct(df)

    return(data.table(df))
}

get_tx_info <- function(df, samps, sample_regex) {
    counts <- data.frame(gene_id=df$gene_id, feature_id=df$ensembl_id,
                         df[,.SD,.SDcols = names(df) %like% sample_regex])
    counts <- distinct(counts)
    d <- dmDSdata(counts=counts, samples=samps)
    return(d)
}

load_ex_data <- function(ex_dir, sample_regex) {
    exs <- list.files(ex_dir)
    exs <- exs[grep(sample_regex, exs)]
    exs <- exs[grep("txt$", exs)]
    exons <- NULL
    featurecounts <- F
    tmp <- read.delim(paste(ex_dir, exs[1], sep='/'), sep='\t', header=T, comment='#')
    if(ncol(tmp) > 2) {featurecounts <- T}

    for(ex in exs) {
        sample <- strsplit(ex, '_')[[1]][1]
        if(featurecounts) {
            e <- read.delim(paste(ex_dir, ex, sep='/'), sep='\t', header=T, comment='#')
            colnames(e)[1] <- 'exon_id'
            e <- e[,c('exon_id', tail(colnames(e), 1))]
        } else {
            e <- read.delim(paste(ex_dir, ex, sep='/'), sep='\t', header=F,
                            col.names = c('exon_id', sample))
        }
        if(is.null(exons)){exons <- e}else{exons <- merge(exons, e, by='exon_id')}
    }

    exons$exon_id <- as.character(exons$exon_id)
    exons <- exons[grep('_', as.character(exons$exon_id), invert=T),]
    exons$gene_id <- sapply(exons$exon_id, function(x){strsplit(x, ':')[[1]][1]})

    return(data.table(exons))
}
