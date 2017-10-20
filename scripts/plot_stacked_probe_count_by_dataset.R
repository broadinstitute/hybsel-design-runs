#!/usr/bin/env Rscript

# Produce a stacked bar showing the probe count for each dataset in a
# probe design.
#
# Args:
#  1: TSV file giving dataset name and number of probes for it
#  2: output PDF file
#  3: [optional] stacked plot limits and tick increment; specified
#     by 'min,max,increment'
#  4: [optional] list of datasets, separated by ',', to label
#     on the stacked chart; if not set, this labels all datasets
#
# By Hayden Metsky <hayden@mit.edu>

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

in.tsv <- args[1]
out.pdf <- args[2]
if (length(args) >= 3) {
    plot.breaks <- as.numeric(strsplit(args[3], ",")[[1]])
    plot.breaks <- seq(plot.breaks[1], plot.breaks[2], plot.breaks[3])
} else {
    plot.breaks <- c('')
}
if (length(args) >= 4) {
    datasets.to.include <- strsplit(args[4], ",")[[1]]
} else {
    datasets.to.include <- c('')
}

# Read counts
counts <- read.table(in.tsv, header=FALSE)
colnames(counts) <- c("dataset", "count")

# The order of the datasets in the stacked chart is given by the
# order of the levels in the factor counts$dataset; to order
# such that the smallest dataset is on top and the largest is
# on bottom of the stack, re-factor counts$dataset so that the
# levels are ordered by counts$count
counts$dataset <- factor(counts$dataset, levels=counts$dataset[order(counts$count)])

# To determine heights of each dataset in the stacked bar,
# first sort in reverse order by count (this sort isn't necessary
# for ggplot2 to correctly create the stacked bar), then use
# a cumulative sum on counts to determine the height (y-value)
# of the center of each dataset's piece in the stacked bar. This
# height will be used to position labels.
counts <- counts[order(counts$count, decreasing=TRUE), ]
counts$stacked.height <- cumsum(counts$count) - 0.5*counts$count

# To just alternate between 2 colors in the stack, add a column
# called color that alternates between them
if (nrow(counts) %% 2 == 0) {
    # since there are an even number of rows, to get #686868 on
    # bottom of the stack (end of the data frame), start with #aaaaaa 
    first.color <- "#aaaaaa"
    second.color <- "#686868"
} else {
    # start with #686868
    first.color <- "#686868"
    second.color <- "#aaaaaa"
}
counts$color <- rep(first.color, nrow(counts))
counts$color[seq(2, nrow(counts), 2)] <- second.color

# If a list of datasets to include was not provided, include all
if (length(datasets.to.include) == 1 && datasets.to.include[1] == '') {
    datasets.to.include <- counts$dataset
}

# Initiate the plot; use stat="identity" so that heights represent
# values
p <- ggplot(counts, aes(x=1, y=count, fill=dataset)) +
     geom_bar(stat="identity")

# Set manual (alternating) colors
p <- p + scale_fill_manual(values=counts$color)

# Hide the legend
p <- p + guides(fill=FALSE)

# Specify the y-label
p <- p + ylab("Number of probes")

# Leave out x-axis text and ticks
# Also leave out the usual ggplot2 background
p <- p + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               panel.background=element_blank(),
               panel.border=element_blank(),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               plot.background=element_blank())

# Add text for the labels
# The label is the dataset name, and this only outputs labels that
# are in datasets.to.include
# This places the label at a y-position given by stacked.height,
# which is the middle of the region corresponding to a dataset in
# the stacked chart
p <- p + geom_text(data=subset(counts, dataset %in% datasets.to.include),
                   aes(x=1, y=stacked.height, label=dataset),
                   size=3)

# If explicit breaks to use in the plot were provided, use those
if (length(plot.breaks) > 1) {
    p <- p + scale_y_continuous(breaks=plot.breaks)
}

# Save the plot
p <- p + ggsave(out.pdf, width=8, height=8, useDingbats=FALSE)
