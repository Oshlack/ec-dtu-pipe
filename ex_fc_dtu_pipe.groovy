load 'ec_dtu_stages.groovy'

run { flatten_gtf + 
      make_star_index +
      fastqFormat * [ star_align + featurecounts_count ] +
      run_dtu.using(feature:'ex')
}
