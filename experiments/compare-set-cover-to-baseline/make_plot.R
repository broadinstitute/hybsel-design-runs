#!/usr/bin/env Rscript

# Plot probe count as a function of number of input
# genomes, for various approaches.
#
# Args:
#   1: TSV file; col 1 gives the name of an approach (e.g.,
#      'Naive approach'); col 2 gives the number of genomes used
#      as input when running the approach; col 3 gives the
#      number of output probes (the same approach and number of
#      genomes can show multiple times, e.g. if there is randomly
#      sampled input)
#   2: output PDF file
#
# By Hayden Metsky <hayden@mit.edu>

library(reshape2)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

in.tsv <- args[1]
out.pdf <- args[2]

##################
# Helper functions
##################

## A helper function from:
##   http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
## Gives count, mean, standard deviation, standard error of the mean, and
## confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be
##     summariezed
##   groupvars: a vector containing names of columns that contain grouping
##     variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is
##     95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count
    # them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

#################
# Read input data
#################
input.table <- read.table(in.tsv, header=FALSE, sep='\t')
colnames(input.table) <- c("approach", "num.genomes", "num.probes")

# Summarize each combination of approach and num.genomes (they
# may appear multiple times, each time run on a different randomly
# sampled input) with summary statistics on num.probes, like the
# mean num.probes across the replicates and a confidence interval
data <- summarySE(input.table, measurevar="num.probes",
                  groupvars=c("approach", "num.genomes"))
# Note that, in data, num.probes is now the mean across the
# replicates for each combination of approach and num.genomes

###############
# Make the plot
###############
p <- ggplot(data, aes(x=num.genomes, y=num.probes, group=approach))
p <- p + geom_line(aes(color=approach))

# Can use geom_errorbar(..) to show error bars at each plotted
# x-value; alternatively, geom_ribbon(..) to show a continuous
# interval (i.e., confidence band) around each line. Note that
# this is a 95% pointwise confidence band, NOT a simultaneous
# confidence band.
p <- p + geom_ribbon(aes(ymin=num.probes-ci, ymax=num.probes+ci,
                     fill=approach), alpha=0.2)

p <- p + labs(title="Scaling probe count with genomes",
              x="Number of genomes",
              y="Number of probes")
#p <- p + scale_y_log10()

# Save to PDF
p <- p + ggsave(out.pdf, width=8, height=8, useDingbats=FALSE)

