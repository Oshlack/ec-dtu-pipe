run_dexseq <- function(counts, genes, group, cpm_cutoff, n_sample_cutoff, cores=8) {
    BPPARAM=MulticoreParam(workers=cores)
    start.time <- Sys.time(); print(start.time)

    counts <- as.matrix(counts[,!colnames(counts)%in%c('ec_names','exon_id','gene_id')])

    sampleTable <- data.frame(condition=as.factor(group))
    rownames(sampleTable) <- colnames(counts)

    dx <- DEXSeqDataSet(counts, sampleData = sampleTable,
                        design=~sample + exon + condition:exon,
                        groupID = genes$gene_id, featureID = genes$feature_id)
    dx <- estimateSizeFactors(dx)
    dx <- estimateDispersions(dx, BPPARAM=BPPARAM)
    dx <- testForDEU(dx, BPPARAM=BPPARAM)
    dx <- estimateExonFoldChanges(dx, BPPARAM=BPPARAM)

    dxr <- DEXSeqResults(dx)
    pgq <- perGeneQValue(dxr); gc()

    end.time <- Sys.time(); print(start.time); gc()
    print(end.time - start.time)

    return(list(dexseq_object=dx,
                dexseq_results=dxr,
                gene_FDR=data.frame(gene=names(pgq), FDR=pgq)))
}

run_diffsplice <- function(df, group, sample_regex, feature=c('tx', 'ec', 'ex'), cpm_cutoff=0, n_sample_cutoff=0) {
    if (feature == 'tx') {
        samps <- data.frame(sample_id = colnames(df)[grep(sample_regex, colnames(df))])
        samps$condition <- as.numeric(as.factor(group)) - 1

        d <- get_tx_info(df, samps, sample_regex)
        samples <- DRIMSeq::samples(d)
        counts <- round(as.matrix(counts(d)[,-c(1:2)]))
        genes <- counts(d)[,c(1:2)]
    } else {
        fname <- ifelse(feature == 'ec', 'ec_names', 'exon_id')
        counts <- distinct(df[,.SD,.SDcols = names(df) %like% paste0(sample_regex, '|gene|', fname)])
        genes <- data.frame(counts)[,c('gene_id', fname)]
        colnames(genes)[2] <- 'feature_id'
        counts <- data.frame(counts)[,grep(sample_regex, colnames(counts))]
    }
    results <- run_dexseq(counts, genes, group, cpm_cutoff, n_sample_cutoff)
    return(results)
}

get_random_comp <- function(samples, n, group1='D2', group2='B6') {
    g1s <- samples[samples$type==group1,]$sample
    g2s <- samples[samples$type==group2,]$sample

    g1 <- sample(g1s, n)
    g2 <- sample(g2s, n)

    return(samples[samples$sample%in%c(g1, g2),])
}
