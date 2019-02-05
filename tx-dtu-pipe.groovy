load 'ec-dtu-stages.groovy'

run { make_salmon_index +
      fastqFormat * [ run_salmon.using(feature:'tx') + link_quant_file ] +
      run_dtu.using(feature:'tx')
}
