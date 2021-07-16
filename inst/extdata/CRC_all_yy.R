library(Seurat)
library(BayesSpace)
set.seed(100)
data_dir <- "/home/yangyi/gmrf_analysis/realdata/Data/CRC3/"
sce <- readVisium(data_dir)

d = 15
sce <- spatialPreprocess(sce, platform="Visium", n.PCs = d, n.HVGs=2000)
position <- colData(sce)[, c("row", "col", "imagerow", "imagecol")]

y <- reducedDim(sce, "PCA")[, seq_len(d)]
all(rownames(poisson) == rownames(y))
all(rownames(poisson) == rownames(colData(sce)[, c("col", "row")]))


tim <- matrix(0, nrow = 1, ncol = 5)
tim2 <- matrix(0, nrow = 1, ncol = 5)
k <- 1
K = 8
library(mclust, quietly=TRUE)
tic <- proc.time()
fit_int = Mclust(y, G=K)
plot(fit_int$BIC)
toc <- proc.time()
x_gmm <- fit_int$classification
tim[k, 1] <- (toc - tic)[1]
tim2[k, 1] <- (toc - tic)[3]
K <- fit_int$G
# -------------------------------------------------
# -------------------------------------------------
# make BayesSpace metadata used in BayesSpace
d = 15
library(BayesSpace)
tic <- proc.time()
scc <- spatialCluster(sce, q=K, d=d, platform='Visium')
toc <- proc.time()

tim[k, 2] <- (toc - tic)[1]
tim2[k, 2] <- (toc - tic)[3]

x_bs <- colData(scc)$spatial.cluster
colData(scc)$x_bs <- x_bs
colData(scc)$x_gmm <- x_gmm
colData(scc)$r_gmm <- fit_int$z


# -------------------------------------------------
# SC.MEB
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

bet = 2.5
tic <- proc.time()
fit_sp <- gmrfICMEM(y, x_int, Adj_sp, mu_int, sigma_int, alpha, bet, TRUE, 10, 50)
toc <- proc.time()

tim[k, 3] <- (toc - tic)[1]
tim2[k, 3] <- (toc - tic)[3]

x_gmrf <- fit_sp$x
colData(scc)$x_gmrf <- x_gmrf
colData(scc)$r_gmrf <- fit_sp$gam


# -------------------------------------------------
# Giotto
## http://spatialgiotto.rc.fas.harvard.edu/giotto.visium.brain.html
library(Giotto)
workdir="/home/yangyi/out" #where results and plots will be saved
myinst=createGiottoInstructions(save_plot=T, show_plot=F, save_dir=workdir, python_path="/home/yangyi/miniconda3/bin/python")

counts = t(y)
pos <- as.matrix(colData(sce)[, c("row", "col")])
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
# rownames(PC) <- paste0("spot_", seq_len(n))
# colnames(PC) <- paste0("PC_", seq_len(p))
gobject@dimension_reduction$cells$pca$pca$coordinates = PC

tic <- proc.time()
## sNN network (default)
gobject <- createNearestNetwork(gobject = gobject, dimensions_to_use = 1:d, k = 30)
## Leiden clustering
gobject <- doLeidenCluster(gobject = gobject, resolution = 1, n_iterations = 1000)
toc <- proc.time()
tim[k, 4] <- (toc - tic)[1]
tim2[k, 4] <- (toc - tic)[3]

x_giotto = gobject@cell_metadata$leiden_clus
write.table(x_giotto, file = "CRC3_giotto.txt", quote = F, col.names = F, row.names = F)
colData(scc)$x_giotto <- x_giotto

# -------------------------------------------------
# louvain
# https://edward130603.github.io/BayesSpace/articles/maynard_DLPFC.html
tic <- proc.time()
g.jaccard = scran::buildSNNGraph(sce, use.dimred="PCA", type="jaccard")
x_louvain <- igraph::cluster_louvain(g.jaccard)$membership
toc <- proc.time()
write.table(x_louvain, file = "CRC3_louvain.txt", quote = F, col.names = F, row.names = F)
colData(scc)$x_louvain <- x_louvain
tim[k, 5] <- (toc - tic)[1]
tim2[k, 5] <- (toc - tic)[3]

setwd("/home/yangyi/gmrf_analysis/realdata/CRC3_out")
saveRDS(scc, paste0("CRC3_ALL_K_", K, "_beta_", bet, ".rds"))
tim_all = rbind(tim, tim2)
write.table(tim_all, file = "CRC3_ALL_tim.txt", quote = F, col.names = F, row.names = F)

