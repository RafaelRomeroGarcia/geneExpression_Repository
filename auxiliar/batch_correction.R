source("https://bioconductor.org/biocLite.R")
biocLite("bladderbatch")
biocLite("pamr")
biocLite("limma")
library(sva)
library(pamr)
library(limma)

edata_RR=read.csv("AIBS_map/downloaded/genes_samples_gen.csv")
batch_RR=read.csv("AIBS_map/downloaded/batch.csv")
batch_RR_r=batch_RR$batch
modcombat_RR = model.matrix(~1, data=batch_RR)
combat_edata_RR=ComBat(dat=edata_RR,batch=batch_RR_r,mod=modcombat_RR,par.prior=TRUE,prior.plots=TRUE)

write.table(combat_edata_RR, file = "AIBS_map/downloaded/genes_samples_gen.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
