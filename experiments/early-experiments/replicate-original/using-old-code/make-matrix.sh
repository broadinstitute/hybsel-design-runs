#!/bin/bash
for s in $(seq 0 2); do
  for m in $(seq 0 2); do
    java -jar RemoveSimilarSequences.jar OFFSET_ALLOWED=${s} MISMATCHES_ALLOWED=${m} I=candidate_probes.fasta O=final_probes.${s}.${m}.fasta
    grep -v '^>' final_probes.${s}.${m}.fasta | sort > final_probes.${s}.${m}.noheader.sorted.fasta
  done
done
