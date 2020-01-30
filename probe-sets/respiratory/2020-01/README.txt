This is for the design of probe sets targeting respiratory-related viruses,
including 2019-nCoV. Note that this does *not* include influenza viruses; there
is a separate probe set for those, which could be combined with this one.

Running CATCH used commit cd0bfa2e3ca918d0079c8f97c55262d2c663958a of the CATCH
repo.

As input, I used 'download:' to downloaded the latest accessions for each tax id;
these are written. The exception is for 2019-nCoV, for which I used a fasta
as input.
 
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
save probes on these runs. Step (5) also adds `--expand-n 0`; it could also add
reverse complements and adapters (but does not here).
The size of the probe set in step (5) might be slightly different than the number
determined in step (2) because of randomness due to `--filter-with-lsh-hamming`.

The final probe sets are in pooled/design-NPROBES/probes.fasta.gz.
