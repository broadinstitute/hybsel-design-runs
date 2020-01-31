This is for the design of probe sets targeting Influenza A, B, and C viruses.

Running CATCH used commit 32147ef602b2d7618ff24eca21b966ed2c4ed58e of the CATCH
repo. (In fact, some was run from commits going back to 93e9f821a29e644f9f6da7df8826a82c077c32af,
but none of the changes between these two commits would affect results. All
calls to CATCH that used --expand-n, which was altered in 32147ef, were run using
32147ef.)

As input, I also used all data from commit 32147ef602b2d7618ff24eca21b966ed2c4ed58e
of the CATCH repo. This is the data that was committed in commit
b12cc141c87e898260ee51a7469a1ac3c890fe93 of the CATCH repo.

This treats each segment as a dataset. Influenza A virus has 8 segments, B has
8 segments, and C has 7 segments. Thus, there are 23 datasets in total.

To generate a probe set, I did the following:
  1) `./run.sh param-exploration-make-commands`
  2) `./run.sh param-exploration-parallel-run NJOBS` where NJOBS is the number of
     jobs to run in parallel
  3) `./run.sh param-exploration-combine-counts`
  4) `./run.sh find-optimal-pooled-params NPROBES` where NPROBES is the number of
     probes in the design
  5) `./run.sh make-probe-set-from-optimal-pooled-params NPROBES NJOBS` where
     NPROBES is the number of probes in the design and NJOBS is the number of jobs
     to run in parallel

I designed a probe set with 6,000 probes and another with 45,000 probes, so I
ran steps (4) and (5) above twice (once with NPROBES=6000 and again with
NPROBES=45000). Note that, since step (5) adds reverse complements to the
probes (via `--add-reverse-complements`), the total size of the probe set is
twice the design size (12,000 and 90,000).

Note that steps (1) and (2) do not in fact output the designed probes; they only
output the number of needed probes (not including reverse complements). They could
have output the probes for each choice of parameters, but I chose instead to
re-run `design.py` in CATCH (in step (5)) on the final choice of parameters and
save probes on these runs. Unlike the design in steps (1) and (2), the design
in step (5) adds adapters and reverse complements; it also adds `--expand-n 0`. The
size of the probe set in step (5) might be slightly different than the number
determined in step (2) (besides also being twice as large, due to adding
reverse complements) because of randomness due to `--filter-with-lsh-hamming`.

The final probe sets are in pooled/design-NPROBES/probes.fasta.gz.

Note that, even though run.sh uses `--expand-n 0` in the design (in the final
step), this was ineffective and some output probes contain N. This was due
to a bug in CATCH, and fixed in commit 074ac1c0.
