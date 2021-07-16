library(Seurat)
library(BayesSpace)
set.seed(100)
library(BPSC)

setwd("/home/yangyi/gmrf_analysis/realdata/CRC3_out");
CRC3_sce = readRDS(paste0("CRC3_K_8_beta_2.5.rds"))

data_dir <- "/home/yangyi/gmrf_analysis/realdata/Data/CRC3/"
sce <- readVisium(data_dir)
d = 15
sce <- spatialPreprocess(sce, platform="Visium", n.PCs = d, n.HVGs=2000)
dec <- scran::modelGeneVar(sce)
top <- scran::getTopHVGs(dec, n = 2000)

idx = match(top, rownames(sce))
data = logcounts(sce)[idx,]

table(colData(CRC3_sce)$x_gmrf)

sum(rownames(colData(CRC3_sce)) == colnames(data))

#Create a data set by merging the control group and the treated group
bp.mat = as.matrix(data)


for (i in 1:8){
  group=c(rep(1,sum(colData(CRC3_sce)$x_gmrf==i)),rep(2,sum(colData(CRC3_sce)$x_gmrf!=i)))
  
  #First, choose IDs of all cells of the control group for estimating parameters of BP models
  controlIds=which(group==2)
  
  #Create a design matrix including the group labels. All batch effects can be also added here if they are available
  design=model.matrix(~group) 
  #Select the column in the design matrix corresponding to the coefficient (the group label) for the GLM model testing
  coef=2 
  
  #Run BPglm for differential expression analysis
  res=BPglm(data=bp.mat, controlIds=controlIds, design=design, coef=coef, estIntPar=FALSE, useParallel=FALSE) 
  #In this function, user can also set estIntPar=TRUE to have better estimated beta-Poisson models for the generalized linear model. However, a longer computational time is required.
  
  #Plot the p-value distribution
  hist(res$PVAL, breaks=20)
  
  #Summarize the results
  ss=summary(res)
  #Compare the discovered DE genes and the true DE genes predefined beforeward
  fdr=p.adjust(res$PVAL, method="BH")
  bpglm.DE.ids=which(fdr<=0.05)
  #Print the indices of the DE genes discovered by BPglm:
  cat(sort(bpglm.DE.ids))
  
  
  write.csv(summary(res)$topTable[summary(res)$topTable[,3]<0.05,], file = paste0("CRC3_betapossion_full_gmrf_cluster_", i, ".csv"))
}




for (i in 1:8){
  group=c(rep(1,sum(colData(CRC3_sce)$x_bs==i)),rep(2,sum(colData(CRC3_sce)$x_bs!=i)))
  
  #First, choose IDs of all cells of the control group for estimating parameters of BP models
  controlIds=which(group==2)
  
  #Create a design matrix including the group labels. All batch effects can be also added here if they are available
  design=model.matrix(~group) 
  #Select the column in the design matrix corresponding to the coefficient (the group label) for the GLM model testing
  coef=2 
  
  #Run BPglm for differential expression analysis
  res=BPglm(data=bp.mat, controlIds=controlIds, design=design, coef=coef, estIntPar=FALSE, useParallel=FALSE) 
  #In this function, user can also set estIntPar=TRUE to have better estimated beta-Poisson models for the generalized linear model. However, a longer computational time is required.
  
  #Plot the p-value distribution
  hist(res$PVAL, breaks=20)
  
  #Summarize the results
  ss=summary(res)
  #Compare the discovered DE genes and the true DE genes predefined beforeward
  fdr=p.adjust(res$PVAL, method="BH")
  bpglm.DE.ids=which(fdr<=0.05)
  #Print the indices of the DE genes discovered by BPglm:
  cat(sort(bpglm.DE.ids))
  
  
  write.csv(summary(res)$topTable[summary(res)$topTable[,3]<0.05,], file = paste0("CRC3_betapossion_full_bs_cluster_", i, ".csv"))
}

