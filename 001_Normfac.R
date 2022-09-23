####Determine normalization factor per donor for SIBERG####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
dirname(rstudioapi::getActiveDocumentContext()$path)

logcounts <- read.csv("log_counts_DESEq_normalized.csv", header = T, row.names = 1)
countsn <- exp(logcounts) -0.5
counts <- read.csv("counts_filtered.csv", header =T, row.names = 1)

nf <- counts/countsn

write.csv(nf, "nf.csv")
