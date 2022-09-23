#Binomial test to test results of split-sets and complete data for significant concordance

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dirname(rstudioapi::getActiveDocumentContext()$path)
set.seed(123)

cvco <- read.csv("class_co_bl.csv", header = T, row.names = 1)
cvsn <- read.csv("class_sn_bl.csv",header = T, row.names = 1)
cvfn <- read.csv("class_fn_bl.csv",header = T, row.names = 1)
all <- read.csv("EV_class_0Gy.csv",header = T, row.names = 1)

#Bashkeel:
#"To this end, a binomial test was used with a probability of success greater than 0.5. 
#For each split-retest replicate, an accurate EV classification, or "success", is defined as a consistent 
#EV classification of a gene between the original data set and both split data sets. 
#Gene classifications were considered significant with a p-value < 0.05 using the binomial test."
#Compare computed EVs from 10*2 split sets to full sample using a binomial test  

#####CO####:

cvcos1 <- cvco[,grepl("s1", colnames(cvco))]
cvcos2 <- cvco[,grepl("s2", colnames(cvco))]

co <- as.data.frame(all$classCO_0Gy)
rownames(co) <- rownames(all)
sf <- function(vec){
  ifelse((vec) == cvcos2 & co, 1,0)}

CO <- as.data.frame(apply(cvcos1, 2, sf))
CO$success <- rowSums(CO)
CO <-as.data.frame((CO[,11]))
rownames(CO) <- rownames(all)
colnames(CO)[1] <- "success"

btco <- lapply(CO$success, binom.test, 10, p= 0.5, alternative = "greater", conf.level= 0.95)
names(btco) <- rownames(CO)

a <- as.data.frame(sapply(btco, "[[","p.value"))
CO <- cbind(CO, a)
colnames(CO)[2] <- "binom.test.p.val"
CO$FDR <- p.adjust(CO$binom.test.p.val, method = "BH")
COf <- CO[CO$FDR <= 0.05,]

write.csv(CO, "co_binom_test.csv")
write.csv(COf, "co_binom_test_fdr.csv")

cof          <- as.data.frame(co[rownames(co) %in% rownames(COf),])
rownames(cof) <- rownames(COf)
colnames(cof)[1] <-"class"
cof$class <- as.factor(cof$class)


COres <- merge(co, cof, by= "row.names",all = T)
colnames(COres)[2] <- "Class.CO"
colnames(COres)[3] <- "Class.CO.cv"

####FPN:####

cvfn1 <- cvfn[,grepl("s1", colnames(cvfn))]
cvfn2 <- cvfn[,grepl("s2", colnames(cvfn))]

fn <- as.data.frame(all$classFPN_0Gy)
rownames(fn) <- rownames(all)
sf <- function(vec){
  ifelse((vec) == cvfn2 & fn, 1,0)}

FPN <- as.data.frame(apply(cvfn1, 2, sf))
FPN$success <- rowSums(FPN)
FPN <-as.data.frame((FPN[,11]))
rownames(FPN) <- rownames(all)
colnames(FPN)[1] <- "success"

bt<- lapply(FPN$success, binom.test, 10, p= 0.5, alternative = "greater", conf.level= 0.95)
names(bt) <- rownames(FPN)

a <- as.data.frame(sapply(bt, "[[","p.value"))
FPN <- cbind(FPN, a)
colnames(FPN)[2] <- "binom.test.p.val"
FPN$FDR <- p.adjust(FPN$binom.test.p.val, method = "BH")
FPNf <- FPN[FPN$FDR <= 0.05,]

write.csv(FPN, "fpn_binom_test.csv")
write.csv(FPNf, "fpn_binom_test_fdr.csv")

fpnf          <- as.data.frame(fn[rownames(fn) %in% rownames(FPNf),])
rownames(fpnf) <- rownames(FPNf)
colnames(fpnf)[1] <-"class"
fpnf$class <- as.factor(fpnf$class)
FPNres <- merge(fn, fpnf, by= "row.names",all = T)
colnames(FPNres)[2] <- "Class.FPN"
colnames(FPNres)[3] <- "Class.FPN.cv"


####SPN:####

cvsn1 <- cvsn[,grepl("s1", colnames(cvsn))]
cvsn2 <- cvsn[,grepl("s2", colnames(cvsn))]

sn <- as.data.frame(all$classSPN_0Gy)
rownames(sn) <- rownames(all)
sf <- function(vec){
  ifelse((vec) == cvsn2 & sn, 1,0)}

SPN <- as.data.frame(apply(cvsn1, 2, sf))
SPN$success <- rowSums(SPN)
SPN <-as.data.frame((SPN[,11]))
rownames(SPN) <- rownames(all)
colnames(SPN)[1] <- "success"

bt<- lapply(SPN$success, binom.test, 10, p= 0.5, alternative = "greater", conf.level= 0.95)
names(bt) <- rownames(SPN)

a <- as.data.frame(sapply(bt, "[[","p.value"))
SPN <- cbind(SPN, a)
colnames(SPN)[2] <- "binom.test.p.val"
SPN$FDR <- p.adjust(SPN$binom.test.p.val, method = "BH")
SPNf <- SPN[SPN$FDR <= 0.05,]

write.csv(SPN, "spn_binom_test.csv")
write.csv(SPNf, "spn_binom_test_fdr.csv")

spnf          <- as.data.frame(sn[rownames(sn) %in% rownames(SPNf),])
rownames(spnf) <- rownames(SPNf)
colnames(spnf)[1] <-"class"
spnf$class <- as.factor(spnf$class)

SPNres <- merge(sn, spnf, by= "row.names",all = T)
colnames(SPNres)[2] <- "Class.SPN"
colnames(SPNres)[3] <- "Class.SPN.cv"

ALL <- cbind(COres, FPNres, SPNres)
rownames(ALL) <- ALL$Row.names
ALL <- ALL[,c(2,3,5,6,8,9)]
write.csv(ALL, "Results_class_cv_0Gy_all.csv")
