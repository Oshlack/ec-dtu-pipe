load 'ec_dtu_stages.groovy'

run { make_salmon_index +
      fastqFormat * [ run_salmon ] +
      run_dtu.using(feature:'tx')
}
