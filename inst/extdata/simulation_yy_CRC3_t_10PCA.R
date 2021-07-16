setwd("/home/yangyi/gmrf_analysis/simulation_CRC3")
library(Rcpp)
library(ggplot2) 
library(Matrix)
library(BayesSpace)
library(SingleCellExperiment)
library(mvtnorm)
library(GiRaF)
library(mritc)
library(gmrf)
library(Giotto)
set.seed(2013)
library(dplyr)


B <- 100
res <- matrix(0, ncol = 5, nrow = B)
#for(b in 1:B){
G <- 4
Bet <- 1
KK <- 8
p <- 10
mu <- matrix(c(rep(0,10),
               rep(0,10),
               rep(0,10),
               rep(0,10),
               rep(0,10),
               rep(0,10),
               rep(0,10),
               rep(0,10)), ncol = KK)

mu[1,1] = 4.75
mu[2,2] = 4.75
mu[3,3] = 4.75
mu[4,4] = 4.75
mu[5,5] = 4.75
mu[6,6] = 4.75
mu[7,7] = 4.75
mu[8:10,8] = 4.75

diagmat = array(0, dim = c(10, 10, 8))
diag(diagmat[,,1]) = c(1,rep(7,9))
diag(diagmat[,,2]) = c(7,1,rep(7,8))
diag(diagmat[,,3]) = c(rep(7,2),1, rep(7,7))
diag(diagmat[,,4]) = c(rep(7,3),1, rep(7,6))
diag(diagmat[,,5]) = c(rep(7,4),1, rep(7,5))
diag(diagmat[,,6]) = c(rep(7,5),1, rep(7,4))
diag(diagmat[,,7]) = c(rep(7,6),1, rep(7,3))
diag(diagmat[,,8]) = c(rep(7,7),1, rep(7,2))

n <- 2988 # # of cell in each indviduals 

replication = 50
ARI = matrix(0, replication, 20)
Time = matrix(0, replication, 6)
Time2 = matrix(0, replication, 6)
output = matrix(0, replication, 5)
colnames(output) = c("gmm_K", "beta", "gmrf_K", "giotto_K", "louvain_K")

data <- read.csv("tissue_positions_list.csv", header = F)
colnames(data)[1:4] = c("barcode", "tissue", "row", "col")
data2 = data %>% filter(tissue == 1)
pos = as.matrix(data2[,3:4])

data_dir = "/home/yangyi/gmrf_analysis/realdata/CRC3_out/"
CRC3 = readRDS(paste0(data_dir, "CRC3_K_8_beta_2.5.rds"))
clu1 = colData(CRC3)$x_gmrf

for (iter in 1:replication){
  print(iter)
  print("###### simulation yy CRC3 8 cluster 10PCA#####")
  
  x = clu1
  
  y <- matrix(0, nrow = n, ncol = p)
  
  for(i in 1:n)	{ # cell
    mu_i <- mu[, x[i]] 
    df_i <- ((x[i]==1)*diagmat[,,1] + (x[i]==2)*diagmat[,,2] + (x[i]==3)*diagmat[,,3] + 
                  (x[i]==4)*diagmat[,,4] + (x[i]==5)*diagmat[,,5] + (x[i]==6)*diagmat[,,6] + 
                  (x[i]==7)*diagmat[,,7] + (x[i]==8)*diagmat[,,8])
    for (j in 1:p){
      y[i, j] <- rt(1, df_i[j, j]) + mu_i[j]
    }
  }

  
  # -------------------------------------------------
  # make BayesSpace metadata used in BayesSpace
  counts <- t(y)
  rownames(counts) <- paste0("gene_", seq_len(p))
  colnames(counts) <- paste0("spot_", seq_len(n))
  
  ## Make array coordinates - filled rectangle
  cdata <- list()
  cdata$row <- pos[,1]
  cdata$col <- pos[,2]
  cdata <- as.data.frame(do.call(cbind, cdata))
  ## Scale and jitter image coordinates
  #scale.factor <- rnorm(1, 8);  n_spots <- n
  #cdata$imagerow <- scale.factor * cdata$row + rnorm(n_spots)
  #cdata$imagecol <- scale.factor * cdata$col + rnorm(n_spots)
  cdata$imagerow <- cdata$row
  cdata$imagecol <- cdata$col 
  ## Make SCE
  ## note: scater::runPCA throws warning on our small sim data, so use prcomp
  sce <- SingleCellExperiment(assays=list(counts=counts), colData=cdata)
  reducedDim(sce, "PCA") <- y
  sce$spatial.cluster <- floor(runif(ncol(sce), 1, 3))
  
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  # -------------------------------------------------
  library(mclust, quietly=TRUE)
  
  tic <- proc.time()
  fit_int = Mclust(y, G = 7:14)
  toc <- proc.time()
  Time[iter, 1] = (toc[1] - tic[1])/8
  Time2[iter, 1] = (toc[3] - tic[3])/8
  mu_int <- unname(fit_int$parameter$mean)
  sigma_int <- unname(fit_int$parameter$variance$sigma) 
  x_gmm <- fit_int$classification
  ARI[iter,1] = adjustedRandIndex(x, x_gmm)
  output[iter, 1] = fit_int$G
  
  # -------------------------------------------------
  for (K in (KK-1):(KK+1)){
    fit_int = Mclust(y, G = K)
    x_gmm <- fit_int$classification
    tic <- proc.time()
    scc <- spatialCluster(sce, q=K, d=p, platform='Visium', init=x_gmm,
                          nrep=10000, gamma=3)
    toc <- proc.time()
    x_bs <- colData(scc)$spatial.cluster
    ARI[iter,K] = adjustedRandIndex(x, x_bs)
  }
  Time[iter, 2] = (toc[1] - tic[1])
  Time2[iter, 2] = (toc[3] - tic[3])
    # -------------------------------------------------
    # markov random fields
    D <- getPairDist(pos)
    source("/home/yangyi/gmrf_analysis/realdata/find_neighbors.R")
    ij <- find_neighbors(sce, "Visium")
    library(Matrix)
    Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1, dims = dim(D))
    
    y_dist <- getPairDist(y)
    Icut <- 2
    
    BIC = matrix(0, 8, 26)
    bet <- seq(0,5,0.2)
    time_gmrf = matrix(0, 8, 2)
    for (K in 7:14){
      fit_int = Mclust(y, G = K)
      x_gmm <- fit_int$classification
      mu_int <- unname(fit_int$parameter$mean)
      sigma_int <- unname(fit_int$parameter$variance$sigma) 
      alpha <- -log(fit_int$parameter$pro)*0
      tic <- proc.time()
      for (b in 1:length(bet)){
        fit_sp <- gmrfICMEM(y, x_gmm, Adj_sp, mu_int, sigma_int, alpha, bet[b], TRUE, 10, 50)
        dr <- K * dim(y)[2] + dim(y)[2]*(dim(y)[2] + 1)/2 * K
        BIC[K-6,b] <- -2 * fit_sp$ell - log(dim(y)[1]) * dr
      }
      toc <- proc.time()
      time_gmrf[K-6,1] = (toc - tic)[1]
      time_gmrf[K-6,2] = (toc - tic)[3]
    }
    Time[iter, 3] = sum(time_gmrf[,1])
    Time2[iter, 3] = sum(time_gmrf[,2])
    
    BIC2 = apply(BIC, 1,  max)
    best_K = which.max(BIC2) + 6
    idx = which.max(BIC[which.max(BIC2),])
    best_beta = seq(0,5,0.2)[idx]
    output[iter, c(2,3)] = c(best_beta, best_K)
    
    fit_int = Mclust(y, G = best_K)
    x_gmm <- fit_int$classification
    mu_int <- unname(fit_int$parameter$mean)
    sigma_int <- unname(fit_int$parameter$variance$sigma) 
    alpha <- -log(fit_int$parameter$pro)*0
    fit_sp <- gmrfICMEM(y, x_gmm, Adj_sp, mu_int, sigma_int, alpha, best_beta, TRUE, 10, 50)
    x_gmrf <- fit_sp$x
    ARI[iter,10] = adjustedRandIndex(x, x_gmrf)
    
  # -------------------------------------------------
  # Giotto
  ## http://spatialgiotto.rc.fas.harvard.edu/giotto.visium.brain.html
  workdir="/home/yangyi/out" #where results and plots will be saved
  myinst=createGiottoInstructions(save_plot=T, show_plot=F, save_dir=workdir, python_path="/home/yangyi/miniconda3/bin/python")
  
  gobject = createGiottoObject(
    raw_exprs = counts,
    spatial_locs = pos,
    norm_expr = NULL,
    norm_scaled_expr = NULL,
    custom_expr = NULL,
    cell_metadata = NULL,
    gene_metadata = NULL,
    spatial_network = NULL,
    spatial_network_name = NULL,
    spatial_grid = NULL,
    spatial_grid_name = NULL,
    spatial_enrichment = NULL,
    spatial_enrichment_name = NULL,
    dimension_reduction = NULL,
    nn_network = NULL,
    images = NULL,
    offset_file = NULL,
    instructions = myinst,
    cores = NA
  )
  PC <- y
  rownames(PC) <- paste0("spot_", seq_len(n))
  colnames(PC) <- paste0("PC_", seq_len(p))
  gobject@dimension_reduction$cells$pca$pca$coordinates = PC
  
  tic <- proc.time()
  ## sNN network (default)
  gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:p, k = 30)
  ## Leiden clustering
  gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)
  toc <- proc.time()
  Time[iter, 4] = toc[1] - tic[1]
  Time2[iter, 4] = toc[3] - tic[3]
  
  x_giotto = gobject@cell_metadata$leiden_clus
  ARI[iter, 11] = adjustedRandIndex(x, x_giotto)
  output[iter, 4] = length(unique(x_giotto))
  
  # -------------------------------------------------
  # louvain
  # https://edward130603.github.io/BayesSpace/articles/maynard_DLPFC.html
  tic <- proc.time()
  g.jaccard = scran::buildSNNGraph(sce, use.dimred="PCA", type="jaccard")
  x_louvain <- igraph::cluster_louvain(g.jaccard)$membership
  ARI[iter, 12] = adjustedRandIndex(x, x_louvain)
  toc <- proc.time()
  
  Time[iter, 5] = toc[1] - tic[1]
  Time2[iter, 5] = toc[3] - tic[3]
  output[iter, 5] = length(unique(x_louvain))
  
  ##--------------------------------------------------
  ## K means
  tic <- proc.time()
  for (j in (KK-1):(KK+1)){
    x_kmeans = kmeans(y, centers=j)$cluster
    ARI[iter, 12 + j - 1] = adjustedRandIndex(x, x_kmeans)
  }
  toc <- proc.time()
  Time[iter, 6] = toc[1] - tic[1]
  Time2[iter, 6] = toc[3] - tic[3]

  
  write.table(ARI, file = paste0("ARI_simulation_CRC3_t_", KK, "cluster_", p, "PCA_48.txt"), quote = F, col.names = F, row.names = F)
  write.table(Time, file = paste0("Time_simulation_CRC3_t_", KK, "cluster_", p, "PCA_48.txt"), quote = F, col.names = F, row.names = F)
  write.table(Time2, file = paste0("Time2_simulation_CRC3_t_", KK, "cluster_", p, "PCA_48.txt"), quote = F,  col.names = F, row.names = F)
  write.table(output, file = paste0("output_simulation_CRC3_t_", KK, "cluster_", p, "PCA_48.txt"), quote = F,  row.names = F)
}
