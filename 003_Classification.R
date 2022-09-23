#Compute expression variability and subsequently classifiy#


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dirname(rstudioapi::getActiveDocumentContext()$path)
library(ggplot2)
library(remotes)
set.seed(123)

#Bootstrap function taken from https://rdrr.io/github/joranE/statskier2/man/bs_mads.html

bs_mads <- function(vec,B = 1000){
  n <- length(vec)
  vec <- vec[!is.na(vec)]
  if (n == 1) return(list(orig = NA,bs = rep(NA,B),n = n))
  orig <- mad(vec)
  mat <- matrix(vec[sample(n,n*B,replace = TRUE)],nrow = n,ncol = B,byrow = FALSE)
  bs_vals <- matrixStats::colMads(mat)
  list(orig = orig,bs = bs_vals,n = n)
}


#load expression data
counts <- read.csv("log_counts_DESEq_normalized.csv", header = T, row.names = 1)
#load metadata
meta <- read.csv("Metadata.csv", header = T)
#unlog
counts <- exp(counts) -0.5
#load bimodal index per gene
BI <- read.csv("bimodal_test_bl.csv", header = T, row.names = 1)
#create filter of bimodally expressed genes
BIf <- BI[BI$BI>=1.1,]

#prepare data
rownames(meta)      <- meta[,1]
meta                <- meta[,2:15]
meta$patient        <- rep((1:156),each = 3)
meta_0 <- meta[meta$Dose=="0 Gy",]
meta_0.05<- meta[meta$Dose=="0.05 Gy",]
meta_2 <- meta[meta$Dose=="2 Gy",]
counts_0Gy            <- counts[,colnames(counts) %in% rownames(meta_0)]
counts_0.05Gy            <- counts[,colnames(counts) %in% rownames(meta_0.05)]
counts_2Gy            <- counts[,colnames(counts) %in% rownames(meta_2)]
counts_0Gy <- data.frame(t(counts_0Gy))
counts_0.05Gy <- data.frame(t(counts_0.05Gy))
counts_2Gy <- data.frame(t(counts_2Gy))
counts_0Gy$Group <- meta_0$Group[match(rownames(counts_0Gy),rownames(meta_0))]
counts_0.05Gy$Group <- meta_0.05$Group[match(rownames(counts_0.05Gy),rownames(meta_0.05))]
counts_2Gy$Group <- meta_2$Group[match(rownames(counts_2Gy),rownames(meta_2))]
rownames(counts_0Gy) <- counts_0Gy$patient
rownames(counts_0.05Gy) <- counts_0.05Gy$patient
rownames(counts_2Gy) <- counts_2Gy$patient

####Examplatory application for 0Gy####
#CO = cancer-free controls
countsCO <- counts_0Gy[counts_0Gy$Group == "CO",]
countsCO <- countsCO[,-ncol(countsCO)]
#SPN = second primary neoplasm
countsSPN <- counts_0Gy[counts_0Gy$Group == "SPN",]
countsSPN <- countsSPN[,-ncol(countsSPN)]
#FPN = first primary neoplasm
countsFPN <- counts_0Gy[counts_0Gy$Group == "FPN",]
countsFPN <- countsFPN[,-ncol(countsFPN)]

#filter bimodally expressed genes from expression data
ug_co <- countsCO[!colnames(countsCO) %in% rownames(BIf)]
###SPN
ug_sn <- countsSPN[!colnames(countsSPN) %in% rownames(BIf)]
###FPN
ug_fn <- countsFPN[!colnames(countsFPN) %in% rownames(BIf)]
#combine to list
list <- list(ug_co, ug_sn, ug_fn)
names(list) <- c("CO_0Gy","SPN_0Gy","FPN_0Gy")
#creat empty EV df
EV <- as.data.frame(colnames(ug_co))
rownames(EV) <- EV$`colnames(ug_co)`
EV <- EV[,0]

for (j in 1:length(list)){
  ##Application of BS-Function (Bootstrap MAD)
  n <- names(list)[j]
  k <- list[[j]]
  tk <- as.data.frame(t(k))
  bs <- apply(k, 2, bs_mads, B=1000)
  tk$mad <- apply(k, 2, mad)
  tk$med <- apply(k,2, median)
  ##Accessing bs-list + median of BS-MAD+
  mad_bs <- as.data.frame(apply(sapply(bs, "[[","bs"), 2, median))
  tk <- cbind(tk, mad_bs)
  #model loess function
  loessmod <- loess(mad~med, tk)
  #predict loess
  tk$loesspred <- predict(loessmod)
  #EV is calculated as the difference between the 
  #observed MAD (bootstrap) and the expected MAD values (Loess).
  tk$EV <- tk$`apply(sapply(bs, "[[", "bs"), 2, median)`- tk$loesspred
  EV <- cbind(EV, tk$EV)
  colnames(EV)[max(ncol(EV))] <- paste(n,sep="")
}


####Classification based on CO EV
#Classify genes x(EV) +/- 3x MAD(bs)
#x(EV) = median of EV per group
#MAD(bs) = median absolute deviation of EV from 1000 bs-replicates
#plot(density(tbig$EV))
#x(EV)
x <- median(EV$CO_0Gy)
#apply MAD-bootstrap to EV
ev_bs <- bs_mads(EV$CO_0Gy)
#median bootstrapped MAD(EV) 
mad_ev_bs <- median(ev_bs$bs)
minrange <- x - (3*mad_ev_bs)
maxrange <- x + (3*mad_ev_bs)
EV$classCO_0Gy  <- as.factor(ifelse(EV$CO_0Gy < minrange, 1, ifelse(EV$CO_0Gy > maxrange,3,2)))
EV$classSPN_0Gy <- as.factor(ifelse(EV$SPN_0Gy < minrange, 1, ifelse(EV$SPN_0Gy > maxrange,3,2)))
EV$classFPN_0Gy <- as.factor(ifelse(EV$FPN_0Gy < minrange, 1, ifelse(EV$FPN_0Gy > maxrange,3,2)))
write.csv(EV,"EV_class_0Gy.csv")

###Empty df for cv
EVcv <- as.data.frame(colnames(ug_co))
rownames(EVcv) <- EVcv$`colnames(ug_co)`
EVcv <- EVcv[,0]
#######Cross validation 50/50, 10 iterations#######
for (i in 1:10){
  # randomly split data in r
  s = sample(seq_len(nrow(ug_co)),size = nrow(ug_co)/2)
  s1co =ug_co[s,]
  s2co =ug_co[-s,]
  s1sn =ug_sn[s,]
  s2sn =ug_sn[-s,]
  s1fn =ug_fn[s,]
  s2fn =ug_fn[-s,]
  list <- list(s1co, s2co, s1sn, s2sn, s1fn, s2fn)
  names(list) <- c("s1co", "s2co", "s1sn", "s2sn", "s1fn", "s2fn")
  for (j in 1:length(list)){
    ##Application of BS-Function (Bootstrap MAD)
    n <- names(list)[j]
    k <- list[[j]]
    tk <- as.data.frame(t(k))
    bs <- apply(k, 2, bs_mads, B=1000)
    tk$mad <- apply(k, 2, mad)
    tk$med <- apply(k,2, median)
    ##Accessing bs-list + median of BS-MAD+
    mad_bs <- as.data.frame(apply(sapply(bs, "[[","bs"), 2, median))
    tk <- cbind(tk, mad_bs)
    #model loess function
    loessmod <- loess(mad~med, tk)
    #predict loess
    tk$loesspred <- predict(loessmod)
    #EV is calculated as the difference between the 
    #observed MAD (bootstrap) and the expected MAD values (Loess).
    tk$EVcv <- tk$`apply(sapply(bs, "[[", "bs"), 2, median)`- tk$loesspred
    EVcv <- cbind(EVcv, tk$EVcv)
    colnames(EVcv)[max(ncol(EVcv))] <- paste(n,"_",i,sep="")
  }
  
}

#store results
write.csv(EVcv,"EV_cv_all_genes.csv")
#EVcv <- read.csv("EV_cv_all_genes.csv",row.names = 1, header = T)
  EVcvco <- EVcv [,grepl("co", names(EVcv))]
  EVcvsn <- EVcv [,grepl("sn", names(EVcv))]
  EVcvfn <- EVcv [,grepl("fn", names(EVcv))]
  
  ####Classification based on CV-data
  #Classify genes x(EV) +/- 3x MAD(bs)
  #x(EV) = median EV 
  #MAD(bs) = median absolute deviation from EV after 1000 bootstrap-replicates
  #Median EV per cv-sample (n=200)
  xcv <- apply(EVcvco, 2, median)
  #MAD-bootstrap to each EV
  ev_bs_cv <- apply(EVcvco, 2, bs_mads)
  #median bootstrapped MAD(EV) 
  mad_ev_bs_cv <- as.data.frame(apply(sapply(ev_bs_cv, "[[","bs"), 2, median))
  minrange_cv <- xcv - (3*mad_ev_bs_cv)
  maxrange_cv <- xcv + (3*mad_ev_bs_cv)
  #classification
  cf <- function(vec){
    as.factor(ifelse((vec) < minrange_cv, 1, ifelse((vec) > maxrange_cv,3,2)))}
  classcocv <- as.data.frame(t(apply(EVcvco, 1, cf)))
  classsncv <- as.data.frame(t(apply(EVcvsn, 1, cf)))
  classfncv <- as.data.frame(t(apply(EVcvfn, 1, cf)))
  colnames(classcocv) <- colnames(EVcvco)
  colnames(classsncv) <- colnames(EVcvsn)
  colnames(classfncv) <- colnames(EVcvfn)
  
  write.csv(classcocv, "class_co_bl.csv")
  write.csv(classsncv, "class_sn_bl.csv")
  write.csv(classfncv, "class_fn_bl.csv")
  
  