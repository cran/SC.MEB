library(RColorBrewer)
library(Seurat)
library(patchwork)
library(dplyr)
library(scCATCH)
library(hdf5r)
library(ggplot2)
library(ggpubr)


##############################################################################################
library(BayesSpace)
set.seed(100)
data_dir <- "/home/yangyi/gmrf_analysis/realdata/Data/CRC3"
sce <- readVisium(data_dir)

d = 15
sce <- spatialPreprocess(sce, platform="Visium", n.PCs = d, n.HVGs=2000)
position <- colData(sce)[, c("row", "col", "imagerow", "imagecol")]

y <- reducedDim(sce, "PCA")[, seq_len(d)]
all(rownames(poisson) == rownames(y))

tim <- matrix(0, nrow = 5, ncol = 3)
k <- 1

Bet = seq(0, 10, 0.1)
B = length(Bet)
BIC <- matrix(0, 10, length(seq(0,10,0.1)))
dr <- numeric(10)
ell <- matrix(0, 10, length(seq(0,10,0.1)))
cluster = matrix(0, dim(y)[1], length(seq(0,10,0.1)))


for(K in 2:10) {
  library(mclust, quietly=TRUE)
  tic <- proc.time()
  #fit_int = Mclust(y, G=K)
  fit_int = Mclust(y, G=K);
  plot(fit_int$BIC)
  toc <- proc.time()
  x_gmm <- fit_int$classification
  tim[k, 1] <- (toc - tic)[1]


  # -------------------------------------------------
  # GMRF
  library(SC.MEB)
  source("/home/yangyi/gmrf_analysis/realdata/find_neighbors.R")
  ij <- find_neighbors(sce, "Visium")
  nrow(ij)/length(x_gmm)

  library(Matrix)
  n <- nrow(y)
  Adj_sp <- sparseMatrix(ij[,1], ij[,2], x = 1, dims = c(n, n))

  alpha <- -log(fit_int$parameter$pro) * 0
  mu_int <- unname(fit_int$parameter$mean)
  sigma_int <- unname(fit_int$parameter$variance$sigma)
  x_int <- x_gmm
  tic <- proc.time()
  for (bet in 1:B)  {
    tic <- proc.time()
    fit_sp <- gmrfICMEM(y, x_int, Adj_sp, mu_int, sigma_int, alpha, Bet[bet], TRUE, 10, 50)
    toc <- proc.time()
    toc - tic
    x_int <- fit_sp$x
    mu_int <- fit_sp$mu
    sigma_int <- fit_sp$sigma
    x_gmrf <- fit_sp$x

    dr[K] <- K * dim(y)[2] + dim(y)[2]*(dim(y)[2] + 1)/2 * K
    ell[K,bet] <- fit_sp$ell
    BIC[K,bet] <- -2 * fit_sp$ell - log(dim(y)[1]) * dr[K]
    cluster[,bet] = x_gmrf

  }
  save(BIC, file = paste0("CRC3_BIC.RData"))
  write.table(cluster, file = paste0("CRC3_gmrf_K_", K, ".txt"), quote = F, row.names = F, col.names = F)

  toc <- proc.time()
}




