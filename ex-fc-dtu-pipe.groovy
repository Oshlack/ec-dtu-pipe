load 'ec-dtu-stages.groovy'

run { flatten_gtf_featurecounts +
      make_star_index +
      fastqFormat * [ star_align + featurecounts_count ] +
      run_dtu.using(feature:'ex_fc')
}
