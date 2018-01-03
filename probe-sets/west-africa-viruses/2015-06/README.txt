20-virus probe set for viruses commonly circulating in West Africa

  - all/ contains probes designed for all viruses (run across
    combinations of parameters)
  - all-with-adapters/ contains a probe set for each virus (one
    selected choice of parameters per virus, as specified in
    params.90000.txt, and copied from all/) with adapters added
    onto the end of each probe; this (all.fasta) is the final
    20-virus, 90000-probe probe set submitted to CustomArray for
    synthesis
  - limited-with-adapters/ has a probe set (all.fasta) with just
    three viruses (based on one selected choice of parameters per
    virus, as specified in params.12000.txt, and coped from all/);
    I don't believe this probe set was ever synthesized or tested
  - num_probes.tsv has the number of probes for each dataset for
    each combination of parameters (mismatches and cover extension)
    as taken from the fasta files in all/
