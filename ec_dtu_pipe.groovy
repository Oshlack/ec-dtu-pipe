load 'ec_dtu_stages.groovy'

run { make_salmon_index + 
      fastqFormat * [ run_salmon ] +
      make_ec_matrix +
      run_dtu
}
