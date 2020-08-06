# This script returns the sum DP and the number of sites for each chromosome
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

chr_data = read.table(args[1], header = TRUE)

chr_data_rmNA = subset(chr_data, chr_data$DP != "NA")

out = data.frame(sum = sum(chr_data_rmNA$DP), num_sites = nrow(chr_data_rmNA))

write.table(out, args[2], quote = F, row.names = F, col.names = F, sep = "\t")