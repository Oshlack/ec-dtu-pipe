load 'ec-dtu-stages.groovy'

run { make_salmon_index +
      fastqFormat * [ run_salmon.using(feature:'ec') ] +
      make_ec_matrix +
      run_dtu.using(feature:'ec')
}
