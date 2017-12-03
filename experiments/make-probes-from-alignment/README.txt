The naive tiling approach yields many probes in large part because unaligned sequences are offset from each other. I was curious how it would do when the input was an alignment of sequences.

To test this, I randomly sampled 100, 500, and 1000 HIV-1 sequences and aligned them using 'mafft --maxiterate 1000 --retree 2'. Then, I replaced gaps ('-') with 'N' so that hybseldesign would not skip over the gaps, allowing the probes to be tiled across the alignment with a fixed stride.

I edited hybseldesign/filter/candidate_probes.py so that the value of min_n_string_length used by make_candidate_probes_from_sequence(..) would be 1000 -- i.e., large enough so that it would design probes with however many Ns were in the alignment. This ensures that it would create candidate probes actually tiled across the alignment at a fixed stride, rather than resetting when it encountered a stretch of Ns.

Then I run ./run_naive_tiling.sh to get the number of probes based on naive tiling from the aligned input and naive tiling from unaligned input. Its output is in naive_tiling.out.

As expected, it appears that using the aligned input results in fewer probes than using unaligned input.
N=100:  aligned input needs 75% as many probes as unaligned input
N=500:  aligned input needs 65% as many probes as unaligned input
N=1000: aligned input needs 62% as many probes as unaligned input

So, in short: it does seem that using aligned input scales better and yields fewer probes, but not *too* many fewer probes (on this data at N=1000, it still needs >50% of the probes).
