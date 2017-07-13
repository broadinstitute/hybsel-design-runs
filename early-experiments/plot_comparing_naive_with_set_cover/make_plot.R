library(ggplot2)

data <- as.data.frame(read.csv("data.csv", header=TRUE))
print(data)

pdf("plot.pdf")

labels <- c("Set cover approach", "Naive approach")
colors <- c("Set cover approach"="#377eb8", "Naive approach"="#656565")
sizes <- c("Set cover approach"=2.0, "Naive approach"=2.0)

ggplot(data, aes(num_genomes)) +
    geom_line(aes(y=set_cover_probes, color="Set cover approach",
        size="Set cover approach")) +
    geom_line(aes(y=naive_probes, color="Naive approach",
        size="Naive approach")) +
    scale_color_manual("", values=colors, breaks=labels) +
    scale_size_manual("", values=sizes, breaks=labels) +
    theme(legend.position=c(0.15, 0.9)) +
    labs(title="Scaling probe count with genomes", y="Number of probes",
        x="Number of Ebola Zaire genomes") +
    scale_x_continuous(breaks=seq(0, 140, 20)) +
    scale_y_continuous(breaks=seq(0,30000,5000), trans="log")

dev.off()
