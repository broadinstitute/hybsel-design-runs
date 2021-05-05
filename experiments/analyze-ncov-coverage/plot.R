#!/usr/bin/env Rscript

require(ggplot2)
require(ggpubr)

args <- commandArgs(trailingOnly=TRUE)

in.tsv <- args[1]
out.pdf <- args[2]

# Read table
input.table <- read.table(gzfile(in.tsv), sep="\t", header=FALSE)
colnames(input.table) <- c("genome", "pos", "depth")

# Remove reverse complement (rc) rows
input.table <- input.table[!grepl("(rc)", input.table$genome),]

# Make the plot
p <- ggplot(input.table, aes(pos, depth))
p <- p + geom_area()
#p <- p + geom_hline(aes(yintercept=0), color="gray")
p <- p + scale_y_continuous(expand=c(0, 0), limits=c(0, NA))    # force the y-axis to start at 0
p <- p + xlab("SARS-CoV-2 genome position") + ylab("Number of probes")
p <- p + theme_pubr()

ggsave(out.pdf, p, width=8, height=4, useDingbats=FALSE)
