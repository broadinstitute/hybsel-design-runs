# Convert a file that lists each sequence (and lineage) for each
# dataset (e.g., seqs-in-201709-pull.txt) into one that lists
# all the lineages for each dataset, separated by ' & '.
# Args:
#  1: input file where col 1 has dataset and col 4 has lineage

for d in $(awk -F'\t' '{print $1}' $1 | sort | uniq); do
    echo -ne "$d\t"; awk -F'\t' -v d="$d" '$1==d {print $4}' $1 | sort | uniq | paste -d'&' -s | sed 's/&/ & /g';
done
