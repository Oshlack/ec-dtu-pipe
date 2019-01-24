plot_bottomly_boxplot <- function(results, cols, title='', toplot=c('FDR', 'TPR'), hline=NA, lines=F) {
    results$feature <- factor(results$feature, levels=names(cols))
    p <- ggplot(results, aes_string('feature', toplot, group='iter', colour='feature')) +
            geom_boxplot(group = 'feature', outlier.shape = NA) +
            geom_jitter(width=0.2, size=0.4) +
            theme_bw() +
            theme(legend.position = 'none',
                  plot.title = element_text(hjust=0.5),
                  text = element_text(size = 18)) +
            ylim(0,1) +
            ylab('') +
            xlab('') +
            ggtitle(title) +
            geom_hline(yintercept = 0.05, colour='grey',  linetype='dotted') +
            scale_color_manual(values = cols)
    if(lines) {
        p <- p + geom_line(alpha=0.3, colour='grey')
    }
    if(!is.na(hline)) {
        p <- p + geom_hline(yintercept = hline, colour='grey',  linetype='dotted')
    }
    return(p)
}
