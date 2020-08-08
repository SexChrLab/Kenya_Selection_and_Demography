#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

data = read.csv(args[1], sep = "\t")

more_than_10_perc = subset(data, data$F_MISS>0.1)

print(nrow(more_than_10_perc)/nrow(data)*100)