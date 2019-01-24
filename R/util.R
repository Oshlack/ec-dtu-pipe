count_features <- function(ecs, txs, exs, species, sample_regex) {
    # make sure there is at least one count (don't need to do ECs as salmon already filters)
    txs <- txs[rowSums(txs[,.SD,.SDcols = names(txs) %like% sample_regex])>0,]
    exs <- exs[rowSums(exs[,.SD,.SDcols = names(exs) %like% sample_regex])>0,]

    ecs_per_gene <- ecs[, length(unique(ec_names)), keyby=list(gene_id)]
    txs_per_gene <- txs[, length(unique(ensembl_id)), keyby=list(gene_id)]
    exs_per_gene <- exs[, length(unique(exon_id)), keyby=list(gene_id)]

    output <- rbind(data.frame(ecs_per_gene, feature='Equivalence classes'),
                    data.frame(txs_per_gene, feature='Transcripts'),
                    data.frame(exs_per_gene, feature='Exons'))
    output$species <- species
    return(output)
}

calculate_feature_variance <- function(df, sample_regex, species='drosophila',
                                       feature='Equivalence classes',
                                       group1=c(1:3), group2=c(4:6)) {
    fname <- 'ec_names'; if(feature!='Equivalence classes'){fname <- 'ensembl'}

    fcols <- paste(sample_regex, fname, 'gene', sep='|')
    counts <- distinct(df[,.SD,.SDcols = names(df) %like% fcols])
    counts <- data.frame(counts[,.SD,.SDcols = names(counts) %like% sample_regex])

    c1 <- cpm(counts)[,group1]
    c2 <- cpm(counts)[,group2]

    c1 <- c1[apply(c1, 1, sum) > 0 & rowSums(c1 > 1) >= 1,]
    c2 <- c2[apply(c2, 1, sum) > 0 & rowSums(c2 > 1) >= 1,]

    conds <- c(rep('c1', nrow(c1)), rep('c2', nrow(c2)))
    vars <- c(apply(c1, 1, var), apply(c2, 1, var))
    means <- c(apply(c1, 1, mean), apply(c2, 1, mean))
    ec <- data.frame(data=feature, condition=conds, variance=vars, mean=means, species=species)
    return(ec)
}

calculate_stats <- function(thresholds, truth, top, gene_id='gene_id', fdr_id='FDR') {
    ist <- intersect(truth$gene, top[,gene_id])
    tr <- subset(truth, gene %in% ist)

    tmp <- top[top[,gene_id] %in% ist, ]
    tmp <- data.table(distinct(tmp[,c(gene_id, fdr_id)]))
    tmp <- suppressWarnings(tmp[,min(get(fdr_id),na.rm=T),by=eval(gene_id)])
    colnames(tmp)[2] <- fdr_id; tmp <- data.frame(tmp)
    if(any(is.infinite(tmp[,fdr_id]))){tmp[is.infinite(tmp[,fdr_id]),fdr_id] <- NA}

    not_ds_genes <- tr$gene[which(tr$ds_status == 0)]
    ds_genes <- tr$gene[which(tr$ds_status == 1)]

    stats <- sapply(thresholds, function(i) {
        sig_genes <- tmp[which(tmp[,fdr_id] <= i), gene_id]
        fp <- length(intersect(sig_genes, not_ds_genes))
        tp <- length(intersect(sig_genes, ds_genes))
        all_tp <- length(ds_genes)
        return(c(fp/length(sig_genes),tp/all_tp))})
    return(stats)
}

get_fdr_tpr_stats <- function(test, truth, results, thresholds, species) {
    res <- NULL
    x <- inner_join(test, truth, species, by='gene')
    methods <- colnames(x)[grep('adjP',colnames(x))]
    for(method in methods) {
        stats <- calculate_stats(thresholds, truth, test, 'gene', method)
        fdrs <- stats[1,]; tprs <- stats[2,]
        res <- rbind(res, data.frame(method=method, FDR=fdrs, TPR=tprs,
                                     thresholds=thresholds, species=species))
    }
    for(i in names(results)) {
        stats <- calculate_stats(thresholds, truth, results[[i]], gene_id = 'gene', fdr_id = 'FDR')
        fdrs <- stats[1,]; tprs <- stats[2,]
        res <- rbind(res, data.frame(method=i, FDR=fdrs, TPR=tprs,
                                     thresholds=thresholds, species=species))
    }
    return(res)
}

get_rank_orders <- function(res, false_dtu, feature, pick_iter=1, glimit=500) {
    sig <- res[[pick_iter]]
    sig <- sig[order(sig$FDR),]
    sig$gene_rank <- rank(sig$FDR)
    sig$false_positives <- NA
    ranks <- unique(sig$gene_rank)
    for(rnk in ranks[ranks <= glimit]) {
        sig[sig$gene_rank == rnk, 'false_positives'] <- sum(sig[sig$gene_rank <= rnk, ]$gene %in% false_dtu)
    }
    sig <- sig[!is.na(sig$false_positives),]
    sig$feature <- feature
    return(sig)
}

get_subset_tests_results <- function(res_ec, res_tx, res_ex,
                                     true_ec, true_tx, true_ex,
                                     cutoff=0.05,
                                     method=c(FDR = 'union', TPR = 'intersect')) {
    results <- NULL
    union_genes <- union(union(true_ec, true_tx), true_ex)
    itx_genes <- intersect(intersect(true_ec, true_tx), true_ex)
    if(method['FDR'] == 'union') {
        true_ec_fdr <- union_genes
        true_tx_fdr <- union_genes
        true_ex_fdr <- union_genes
    } else if(method['FDR'] == 'intersect') {
        true_ec_fdr <- itx_genes
        true_tx_fdr <- itx_genes
        true_ex_fdr <- itx_genes
    } else { # assume individual
        true_ec_fdr <- true_ec
        true_tx_fdr <- true_tx
        true_ex_fdr <- true_ex
    }
    if(method['TPR'] == 'union') {
        true_ec_tpr <- union_genes
        true_tx_tpr <- union_genes
        true_ex_tpr <- union_genes
    } else if(method['TPR'] == 'intersect') {
        true_ec_tpr <- itx_genes
        true_tx_tpr <- itx_genes
        true_ex_tpr <- itx_genes
    } else { # assume individual
        true_ec_tpr <- true_ec
        true_tx_tpr <- true_tx
        true_ex_tpr <- true_ex
    }

    for(i in 1:length(res_ec)) {
        sig_ec <- res_ec[[i]]
        sig_ec <- unique(sig_ec[sig_ec$FDR < cutoff,]$gene)
        sig_ec_true <- sum(sig_ec%in%true_ec_tpr)
        sig_ec_false <- sum(!sig_ec%in%true_ec_fdr)
        results <- rbind(results, data.frame(feature='Equivalence classes',
                                             iter=i,
                                             true_sig=sig_ec_true,
                                             false_sig=sig_ec_false,
                                             FDR=sig_ec_false / length(sig_ec),
                                             TPR=sig_ec_true / length(true_ec_tpr)))

        sig_tx <- res_tx[[i]]
        sig_tx <- unique(sig_tx[sig_tx$FDR < cutoff,]$gene)
        sig_tx_true <- sum(sig_tx%in%true_tx_tpr)
        sig_tx_false <- sum(!sig_tx%in%true_tx_fdr)
        results <- rbind(results, data.frame(feature='Transcripts',
                                             iter=i,
                                             true_sig=sig_tx_true,
                                             false_sig=sig_tx_false,
                                             FDR=sig_tx_false / length(sig_tx),
                                             TPR=sig_tx_true / length(true_tx_tpr)))

        sig_ex <- res_ex[[i]]
        sig_ex <- unique(sig_ex[sig_ex$FDR < cutoff,]$gene)
        sig_ex_true <- sum(sig_ex%in%true_ex_tpr)
        sig_ex_false <- sum(!sig_ex%in%true_ex_fdr)
        results <- rbind(results, data.frame(feature='Exons',
                                             iter=i,
                                             true_sig=sig_ex_true,
                                             false_sig=sig_ex_false,
                                             FDR=sig_ex_false / length(sig_ex),
                                             TPR=sig_ex_true / length(true_ex_tpr)))
    }
    return(results)
}
