#Computation of bimodal index#

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dirname(rstudioapi::getActiveDocumentContext()$path)

library(SIBERG)
library(dplyr)

#load unnormalized un-log counts
counts <- read.csv("counts_filtered.csv", header =T, row.names = 1) 
#load meta data
meta <- read.csv("Metadata.csv", header = T)
#load vector with normalization factors created in 001_Normfac
nf <- read.csv("nf.csv", header = T, row.names = 1)

#get nf per donor
nf[nf == 0] <- NA
nf <- nf %>% filter(complete.cases(.))
nf <- nf[1,]

#prepare data
rownames(meta)      <- meta[,1]
meta                <- meta[,2:15]
meta$patient        <- rep((1:156),each = 3)

meta_0 <- meta[meta$Dose =="0 Gy",]
meta_0.05<- meta[meta$Dose=="0.05 Gy",]
meta_2 <- meta[meta$Dose=="2 Gy",]


counts_0Gy            <- as.matrix(counts[,colnames(counts) %in% rownames(meta_0)])
counts_0.05Gy            <- as.matrix(counts[,colnames(counts) %in% rownames(meta_0.05)])
counts_2Gy            <- as.matrix(counts[,colnames(counts) %in% rownames(meta_2)])
nf0                   <- t(nf[,colnames(nf) %in% rownames(meta_0)])
nf0.05                   <- t(nf[,colnames(nf) %in% rownames(meta_0.05)])
nf2                   <- t(nf[,colnames(nf) %in% rownames(meta_2)])


#Calculate bimodal index for data of each radiation dose 
hdir <- apply(counts_2Gy, 1, SIBER,d=1/nf2, model= 'LN', zeroPercentThr = 0.2, base=exp(1), eps=10)
ldir <- apply(counts_0.05Gy, 1, SIBER,d=1/nf0.05, model= 'LN', zeroPercentThr = 0.2, base=exp(1), eps=10)
bl <- apply(counts_0Gy, 1, SIBER,d=1/nf0, model= 'LN', zeroPercentThr = 0.2, base=exp(1), eps=10)

hdir <- as.data.frame(t(hdir))
ldir <- as.data.frame(t(ldir))
bl <- as.data.frame(t(bl))

#extract numbers for manuscript
hdirf <- hdir[hdir$BI>= 1.1,]
ldirf <- ldir[ldir$BI>= 1.1,]
blf <- ldir[bl$BI>= 1.1,]

#store all results for each radiation dose
write.csv(hdir,"bimodal_test_hdir.csv")
write.csv(ldir,"bimodal_test_ldir.csv")
write.csv(bl,"bimodal_test_bl.csv")



