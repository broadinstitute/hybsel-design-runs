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

The final probe sets are in pooled/design-NPROBES/probes.fasta.gz.

Note that the taxonomic ID for hpiv_4 is 11224, which is Human parainfluenza virus 4a
(not 1979161, which is Human rubulavirus 4). This is because the NCBI seems to have
some issues, where 1979161 does not return any results but redirects to 11224. 11224
in fact includes neighbors for both Human parainfluenza 4a and 4b (both major strains),
so using it is fine; this can be verified by looking under the taxid-acc/ for hpiv_4.
Likewise, the ID for human_mastadenovirus_g is 310540 (Simian adenovirus 1) rather than
536079 for this same reason; but 310540 includes neighbors for all the major lineages
of Human mastadenovirus G.
