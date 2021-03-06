make_salmon_index = {
    output.dir = salmon_index

    produce('*.bin') {
        exec """
        $time $salmon index -t $txome -i $salmon_index
        """
    }
}

run_salmon = {
    def base_outdir = branch.name + '_' + feature + '/salmon_out'
    def skipquantarg = feature == 'ec' && skipQuant.toBoolean() ? '--skipQuant' : ''

    def reads = '-r ' + inputs
    if (reads.size() > 1) {
        reads = '-1 ' + inputs[0] + ' -2 ' + inputs[1]
    }

    output.dir = base_outdir + '/aux_info'
    produce('eq_classes.txt*') {
        exec """
        $time $salmon quant --seqBias --gcBias --dumpEq \
            --index $salmon_index -l A $reads -p $cores -o $base_outdir $skipquantarg
        """, 'run_salmon'
    }
}

link_quant_file = {
    def workingDir = System.getProperty('user.dir');
    def base_outdir = branch.name + '_tx/salmon_out'
    output.dir = 'quant'

    produce(branch.name + '_quant.sf') {
        exec """
        ln -s $workingDir/$base_outdir/quant.sf $output
        """, 'link_quant_file'
    }
}

flatten_gtf = {
    idx = tx_gtf.lastIndexOf('.')
    basename = idx != -1 ? tx_gtf[0..<idx] : tx_gtf
    basename = basename.split('/')[-1]

    produce(basename + '.gff') {
        exec """
        $time $python $dexseq/dexseq_prepare_annotation.py --aggregate='no' $tx_gtf $output
        """
    }
}

flatten_gtf_featurecounts = {
    def workingDir = System.getProperty('user.dir');
    idx = tx_gtf.lastIndexOf('.')
    basename = idx != -1 ? tx_gtf[0..<idx] : tx_gtf
    basename = basename.split('/')[-1]

    produce(basename + '.featurecounts.gtf') {
        exec """
        R CMD BATCH --no-restore --no-save "--args input_gtf=\'$tx_gtf\' output_gtf=\'$output\' ignore_strand=TRUE" $RCODEGEN/generate_flattened_gtf.R $workingDir/generate_flattened_gtf.Rout ;
        """
    }
}

make_star_index = {
    output.dir = star_index

    produce('Genome') {
        exec """
        $time STAR --runMode genomeGenerate \
              --genomeDir $star_index \
              --sjdbGTFfile $tx_gtf \
              --genomeFastaFiles $genome \
              --runThreadN $cores \
              --genomeSAindexNbases 5
        """, 'make_star_index'
    }
}

star_align = {
    produce(branch.name + 'Aligned.sortedByCoord.out.bam'){
        exec """
        rm -rf ${branch.name}_STARtmp ;
        $time STAR --genomeDir $star_index
           --readFilesCommand zcat
           --readFilesIn $inputs
           --outSAMtype BAM SortedByCoordinate
           --outFileNamePrefix $branch.name
           --runThreadN $cores
           --limitBAMsortRAM $genome_mem ;
        $time samtools index $output
        """, 'star_align'
   }
}

dexseq_count = {
    output.dir = 'exon_counts'
    paired = inputs.gz.split().size() == 2 ? 'yes' : 'no'

    produce(branch.name + '_count.txt') {
        exec """
        $time $python $dexseq/dexseq_count.py \
            -p $paired -r pos -s no -f bam $input.gff \
            $input.bam $output
        """, 'dexseq_count'
    }
}

featurecounts_count = {
    output.dir = 'fc_counts'

    produce(branch.name + '_count.txt') {
        exec """
        $time $featurecounts -T $cores -t exon -g exon_id -a $input.gtf -o $output $input.bam
        """, 'featurecounts_count'
    }
}

make_ec_matrix = {
    def sample_names = inputs.txt.split().collect { it.split('/')[-4] }
    sample_names = sample_names.join(',')

    produce('ec_count_matrix.txt') {
        exec """
        $time $python $pipe_dir/create_salmon_ec_count_matrix.py $inputs $sample_names $output
        """, 'create_ec_count_matrix'
    }
}

run_dtu = {
    def dat = input
    def outname = feature
    if (feature == 'tx') {
        dat = 'quant/'
    } else if (feature == 'ex') {
        dat = 'exon_counts/'
    } else if (feature == 'ex_fc') {
        dat = 'fc_counts/'
    }

    def optargs = ''
    if (tx_ref != '$tx_ref' && tx_ref != '') {
        optargs = optargs + ' --ref ' + tx_ref
    }
    if (tx_lookup != '$tx_lookup' && tx_lookup != '') {
        optargs = optargs + ' --lookup ' + tx_lookup
    }
    if (bmart_dset != '$bmart_dset' && bmart_dset != '') {
        optargs = optargs + '--biomart ' + bmart_dset
    }

    produce(outname + '_results.RData') {
        exec """
        $time Rscript $pipe_dir/run_dtu.R $feature $dat $group $sample_regex $output $optargs --cores $cores
        """, 'run_dtu'
    }
}
