#!/usr/bin/env Rscript

# Produce a bubble plot showing each set of parameters in probe
# design (mismatches, cover extension), with the size of each
# bubble indicating the number of datasets that have that choice
# of parameters in the design.
#
# Args:
#  1: TSV file giving number of mismatches, cover extension, and
#     the number of datasets in the probe design with that pair
#     of parameters
#  2: output PDF file
#  3: breaks on the x-axis, separated by commas
#  4: breaks on the y-axis, separated by commas
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

in.tsv <- args[1]
out.pdf <- args[2]
x.breaks <- as.numeric(strsplit(args[3], ",")[[1]])
y.breaks <- as.numeric(strsplit(args[4], ",")[[1]])

# Read table
params <- read.table(in.tsv, header=FALSE)
colnames(params) <- c("mismatches", "cover.extension", "num.datasets")

# Initiate the plot
p <- ggplot(params, aes(x=mismatches, y=cover.extension, size=num.datasets)) +
     geom_point()

# Set the axis breaks as provided
p <- p + scale_x_continuous(breaks=x.breaks, limits=c(min(x.breaks), max(x.breaks))) +
         scale_y_continuous(breaks=y.breaks, limits=c(min(y.breaks), max(y.breaks)))

# Scale by area with a larger than default max point size
p <- p + scale_size_area(max_size=20)

# Put text labels on each bubble that give the number of datasets
# with the choice of parameters
p <- p + geom_text(data=params,
                   aes(x=mismatches, y=cover.extension, label=num.datasets),
                   color="white", size=2)

# Specify axis labels
p <- p + xlab("Mismatches") + ylab("Cover extension")

# Hide the legend
p <- p + guides(size=FALSE)

# Leave out usual ggplot2 background and minor grid; put
# in axis lines
p <- p + theme(panel.background=element_blank(),
               panel.border=element_blank(),
               axis.line = element_line(color="black"))

# Save the plot
p <- p + ggsave(out.pdf, width=8, height=4, useDingbats=FALSE)
