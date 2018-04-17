#!/usr/bin/env Rscript

# Plot number of probes, time (in sec), and memory (in kb) for
# different datasets over different sample sizes.
#
# Args:
#   1: TSV file; col 1 gives name of an approach/data (dataset and fixed
#      parameter values); col 2 gives sample size (number of input genomes);
#      col 3 gives what LSH family was used (or 'nolsh' if LSH is not used);
#      col 4 gives number of probes; col 5 gives user runtime (in sec);
#      col6 gives memory (MRSS) required in kb
#   2: output PDF file
#
# By Hayden Metsky <hayden@mit.edu>

library(reshape2)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

in.tsv <- args[1]
out.pdf <- args[2]


#################
# Read input data
#################
input.table <- read.table(in.tsv, header=TRUE, row.names=NULL, sep='\t')

# Melt the data so we can show facets
stats <- melt(input.table, id=c("data", "sample.size", "lsh"))

###############
# Make the plot
###############
p <- ggplot(stats, aes(x=sample.size, y=value))
p <- p + geom_line(aes(color=lsh))

# Facet so that each row corresponds to data and each column corresponds
# to a variable
p <- p + facet_wrap(data ~ variable, ncol=3, scales="free")

p <- p + theme(strip.text=element_text(size=3),
               text=element_text(size=3))

# Save to PDF
p <- p + ggsave(out.pdf, width=4, height=12, useDingbats=FALSE)
