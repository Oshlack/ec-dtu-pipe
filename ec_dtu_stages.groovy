salmon='salmon'
code_base='/group/bioi1/marekc/20180202_ec_dtu/ec_dtu_pipeline'
threads=8

make_salmon_index = {
    output.dir = txome_index

    produce('hash.bin'){
        exec """
        $salmon index -t $txome -i $txome_index
        """
    }
}

run_salmon = {
    def (rf1, rf2) = inputs.fq.gz.split().collect { $it }
    def base_outdir = branch.name + '/salmon_out'
    output.dir = branch.name + '/salmon_out/aux_info'

    produce('eq_classes.txt'){
        exec """
        $salmon quant --dumpEq --seqBias -i $txome_index -l A -r $rf1 $rf2 -p $threads -o $base_outdir
        """, 'run_salmon'
    }
}

make_ec_matrix = {
    def sample_names=inputs.split().collect { it.split('/')[-4] }
    sample_names = sample_names.join(',')
   
    produce('ec_count_matrix.txt'){
        exec """
        python $code_base/create_salmon_ec_count_matrix.py $inputs $sample_names $output
        """, 'create_ec_count_matrix'
    }
}

run_dtu = {
    produce('eq_class_comp_diffsplice.txt'){
        exec """
        Rscript $code_base/run_dtu.R $feature $input.txt $group $sample_regex $tx_ref $tx_lookup ;
        """, 'run_dtu'
    }
}

