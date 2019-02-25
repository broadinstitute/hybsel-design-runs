This is for the design of probe sets targeting 588 viral species known to infect humans.

Running CATCH used commit 80641b28df7a17e84b003de550409b0cfcc6ace8 of the CATCH
repo.

As input, I also used all data from commit 80641b28df7a17e84b003de550409b0cfcc6ace8
of the CATCH repo. This is the data that was committed in commit
b12cc141c87e898260ee51a7469a1ac3c890fe93 of the CATCH repo.

This treats each species as a dataset. For Influenza A, B, and C viruses, it treats
each segment as a dataset. Influenza A virus has 8 segments, B has 8 segments,
and C has 7 segments. Thus, there are 608 datasets in total. These are listed
in datasets.txt.

To generate a probe set, I did the following:
  1) `./run.sh param-exploration-make-commands`
  2) `./run.sh param-exploration-parallel-run NJOBS` where NJOBS is the number of
     jobs to run in parallel
  3) `./run.sh param-exploration-combine-counts`
  4) `./run.sh find-optimal-pooled-params NPROBES NJOBS` where NPROBES is the number of
     probes in the design and NJOBS is the number of jobs to run in parallel
  5) `./run.sh make-probe-set-from-optimal-pooled-params NPROBES` where NPROBES is
     the number of probes in the design

I designed probe sets with 250,000 and 350,000 and 700,000 probes, so I ran steps
(4) and (5) three times (each with the corresponding value for NPROBES).

Note that steps (1) and (2) do indeed output designed probes for each choice of
parameter values. Step (3) reads the number of probes designed from the output,
and step (4) searches for an optimal combination of parameter values given those
numbers of probes. Step (5) combines probe outputs from step (2) based on the
optimal combination of parameter values, as determined in step (4).

The final probe sets are in pooled/design-NPROBES/probes.fasta.gz.
