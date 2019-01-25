code_base='/group/bioi1/marekc/20180202_ec_dtu/ec_dtu_pipeline'

salmon='salmon'
star='STAR'
time='/usr/bin/time -v'
python2='/group/bioi1/marekc/apps/conda2/bin/python'
featurecounts='/group/bioi1/marekc/apps/subread-1.6.3-Linux-x86_64/bin/featureCounts'

threads=8
genome_mem=5050000000

make_salmon_index = {
    output.dir = salmon_index

    produce('hash.bin') {
        exec """
        $time $salmon index -t $txome -i $salmon_index
        """
    }
}

run_salmon = {
    def workingDir = System.getProperty('user.dir');
    def base_outdir = branch.name + '/salmon_out'
    output.dir = base_outdir + '/aux_info'

    produce('eq_classes.txt') {
        exec """
        $time $salmon quant --dumpEq --seqBias -i $salmon_index -l A -r $inputs -p $threads -o $base_outdir ;
        mkdir -p quant; ln -s $workingDir/$base_outdir/quant.sf quant/${branch.name}_quant.sf
        """, 'run_salmon'
    }
}

flatten_gtf = {
    idx = tx_gtf.lastIndexOf('.')
    basename = idx != -1 ? tx_gtf[0..<idx] : tx_gtf
    basename = basename.split('/')[-1]

    produce(basename + '.gff') {
        exec """
        $time $python2 $dexseq/dexseq_prepare_annotation.py --aggregate='no' $tx_gtf $output
        """
    }
}

make_star_index = {
    output.dir = star_index

    produce('Genome') {
        exec """
        $time $star --runMode genomeGenerate \
              --genomeDir $star_index \
              --sjdbGTFfile $tx_gtf \
              --genomeFastaFiles $genome \
              --runThreadN $threads \
              --genomeSAindexNbases 5
        """, 'make_star_index'
    }
}

star_align = {
    produce(branch.name + 'Aligned.sortedByCoord.out.bam'){
        exec """
        rm -rf ${out_prefix}_STARtmp ;
        $time $star --genomeDir $star_index
           --readFilesCommand zcat
           --readFilesIn $inputs
           --outSAMtype BAM SortedByCoordinate
           --outFileNamePrefix $branch.name
           --runThreadN $threads
           --limitBAMsortRAM $genome_mem ;
        $time samtools index $output
        """, 'star_align'
   }
}

dexseq_count = {
    output.dir = 'exon_counts'

    produce(branch.name + '_count.txt') {
        exec """
        mkdir -p $output.dir ;
        $time $python2 $dexseq/dexseq_count.py \
            -p yes -r pos -s no -f bam $input.gff \
            $input.bam $output
        """, 'dexseq_count'
    }
}

featurecounts_count = {
    output.dir = 'fc_counts'

    produce(branch.name + '_count.txt') {
        exec """
        mkdir -p $output.dir ;
        $time $featurecounts -T $threads -t exon -g exon_id -a $tx_gtf -o $output $input.bam
        """, 'featurecounts_count'
    }
}

make_ec_matrix = {
    def sample_names = inputs.txt.split().collect { it.split('/')[-4] }
    sample_names = sample_names.join(',')

    produce('ec_count_matrix.txt') {
        exec """
        $time python $code_base/create_salmon_ec_count_matrix.py $inputs $sample_names $output
        """, 'create_ec_count_matrix'
    }
}

run_dtu = {
    //TODO: fix group -- should be sample names in group 1 rather than binary (can't guarantee order for txs and exs)
    def dat = input
    def outname = feature
    if (feature == 'tx') {
        dat = 'quant/'
    } else if (feature == 'ex') {
        dat = 'exon_counts/'
    } else if (feature == 'ex_fc') {
        feature = 'ex'
        dat = 'fc_counts/'
    }

    produce(outname + '_results.RData') {
        exec """
        $time Rscript $code_base/run_dtu.R $feature $dat $group $sample_regex $output $tx_ref $tx_lookup ;
        """, 'run_dtu'
    }
}
