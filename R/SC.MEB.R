#' SC.MEB.
#' 
#' @description
#' SC.MEB implements the model SC-MEB, spatial clustering with hidden Markov random field using empirical Bayes. 
#'
#' @details SC.MEB implements the model SC-MEB, spatial clustering with hidden Markov random field using empirical Bayes. 
#' @param sce is a SingleCellExperiment object containing PCA and position informatin. 
#' @param d is a integer specifying the dimension of PCA. The default is 15. 
#' @param K is an integer vector specifying the numbers of mixture components (clusters) for which the BIC is to be calculated. The default is K = 2:9. 
#' @param platform is the name of spatial transcriptomic platform. Specify 'Visium' for hex lattice geometry or 'ST' for square lattice geometry. Specifying this parameter is optional as this information is included in their metadata. 
#' @param bet is a numeric vector specifying the smoothness of Random Markov Field. The default is seq(0,5,0.2). 
#' @param maxIter_ICM is the maximum iteration of ICM algorithm. The default is 10. 
#' @param maxIter is the maximum iteration of EM algorithm. The default is 50.
#' @return a list, We briefly explain the output of the SC.MEB. 
#' 
#' The item 'best_K' is the optimal K we choose according to BIC rule. 
#' 
#' The item 'best_beta' is also the optimal beta we choose according to BIC rule. 
#' 
#' The item 'best_cluster_label' is the optimal clustering result corresponding to optimal K and optimal beta. 
#' 
#' The item 'best_BIC' is the optimal BIC corresponding to optimal K and optimal beta. 
#' 
#' The item 'best_ell' is the optimal opposite log-likelihood corresponding to optimal K and optimal beta. 
#' 
#' The item 'best_mu' is the optimal mean for each component corresponding to optimal K and optimal beta. 
#' 
#' The item 'best_sigma' is the optimal variance for each component corresponding to optimal K nd optimal beta. 
#' 
#' The item 'best_gam' is the optimal posterior probability matrix corresponding to optimal K and optimal beta.
#'  
#' The item 'cluster_label' is 3-dimensional n-by-b-by-q matrix, storing all clustering results for each K and beta. n is the number of cells, b is the length of vector 'bet', q is the length of vector 'K'. 
#' 
#' The item 'BIC' contains all BIC value for each K and beta. 
#' 
#' The item 'ell' is the opposite log-likelihood for each beta and K. 
#' 
#' The item 'mu' is the mean of each component for each beta and K.
#' 
#' The item 'sigma' is the variance of each component for each beta and K.
#' 
#' @references Yang Y, Shi X, Zhou Q, et al. SC-MEB: spatial clustering with hidden Markov random field using empirical Bayes[J]. bioRxiv, 2021.
#' @import mclust
#' @examples
#' ## read the simulated data 
#' data(sce)
#' d = 5
#' K = 3
#' bet = 1
#' platform = "ST"
#' maxIter_ICM = 10
#' maxIter = 50
#' ## run the SC.MEB function
#' out = SC.MEB(sce, d = d, K = K, bet=bet, platform = platform)
#' @export
SC.MEB <- function(sce, d=15, K=2:9, bet=seq(0,5,0.2), platform = "Visium", maxIter_ICM = 10,  maxIter = 50) {
  # select K and select beta
  y <- SingleCellExperiment::reducedDim(sce, "PCA")[, seq_len(d)]
  pos <- as.matrix(SingleCellExperiment::colData(sce)[, c("row", "col")])
  Adj_sp <- find_neighbors2(sce, platform)

  BIC = matrix(0, length(K), length(bet))
  rownames(BIC) = paste0("K = ", K)
  colnames(BIC) = paste0("beta = ", bet)
  ell = matrix(0, length(K), length(bet))
  rownames(ell) = paste0("K = ", K)
  colnames(ell) = paste0("beta = ", bet)
  cluster_label = array(0, c(dim(y)[1], length(bet), length(K)))
  mu = list(1)
  sigma = list(1)
  gam = list(1)
  pxgn = list(1)
  pygx = list(1)

  for (k in K){
    fit_int = mclust::Mclust(y, G=k, verbose = FALSE)
    mu_int <- unname(fit_int$parameter$mean)
    sigma_int <- unname(fit_int$parameter$variance$sigma)
    x_int <- fit_int$classification
    alpha <- -log(fit_int$parameter$pro) * 0
    index_k = match(k, K)

    mu_k = array(0, c(dim(y)[2], k, length(bet)))
    sigma_k = array(0, c(dim(y)[2], dim(y)[2], k, length(bet)))
    gam_k = array(0, c(dim(y)[1], k, length(bet)))
    pxgn_k = array(0, c(dim(y)[1], k, length(bet)))
    pygx_k = array(0, c(dim(y)[1], k, length(bet)))

    for (b in bet)
    {
      fit_sp <- gmrfICMEM(y, x_int, Adj_sp, mu_int, sigma_int, alpha, b, TRUE, maxIter_ICM, maxIter)
      dr <- k * dim(y)[2] + dim(y)[2]*(dim(y)[2] + 1)/2 * k
      index_b = match(b, bet)
      BIC[index_k, index_b] <- -2 * fit_sp$ell - log(dim(y)[1]) * dr
      cluster_label[ , index_b, index_k] = fit_sp$x
      ell[index_k, index_b] = fit_sp$ell

      mu_k[ , ,index_b] = fit_sp$mu
      sigma_k[ , , ,index_b] = fit_sp$sigma
      gam_k[ , ,index_b] = fit_sp$gam
      pxgn_k[ , ,index_b] = fit_sp$pxgn
      pygx_k[ , ,index_b] = fit_sp$pygx
    }
    mu[[index_k]] = mu_k
    sigma[[index_k]] = sigma_k
    gam[[index_k]] = gam_k
    pxgn[[index_k]] = pxgn_k
    pygx[[index_k]] = pygx_k

  }

  BIC2 = apply(BIC, 1,  max)
  idx1 = which.max(BIC2)
  best_K = K[idx1]
  idx2 = which.max(BIC[which.max(BIC2),])
  best_beta = bet[idx2]

  out = list()
  out$best_K = best_K
  out$best_beta = best_beta
  out$best_cluster_label = cluster_label[, idx2, idx1]
  out$best_BIC = BIC[idx1, idx2]
  out$best_ell = ell[idx1, idx2]
  out$best_mu = mu[[idx1]][ , ,idx2]
  out$best_sigma = sigma[[idx1]][ , , ,idx2]
  out$best_gam = gam[[idx1]][ , ,idx2]
  out$cluster_label = cluster_label
  out$BIC = BIC
  out$ell = ell
  out$mu = mu
  out$sigma = sigma
  out$gam = gam
  out$pxgn = pxgn
  out$pygx = pygx
  out
}

#' gmrfICMEM.
#' 
#' @description
#' gmrfICMEM was used to conduct spatial clustering with hidden Markov random field for fixed beta and fixed number of clusters
#'
#' @details gmrfICMEM was used to conduct spatial clustering with hidden Markov random field for fixed beta and fixed number of clusters
#' @param y is a matrix of PCs containing gene expression.
#' @param x_int is a vector of initial cluster label.
#' @param Adj is a matrix containing neighborhood information generated by find_neighbors2.
#' @param mu_int is a initial mean vector. we often generated it by Gaussian mixture model. 
#' @param sigma_int is a initial co-variance matrix. we often generated it by Gaussian mixture model. 
#' @param alpha is a intercept.
#' @param beta is a smoothing parameter that can be specified by user.
#' @param PX is a logical value specifying the parameter expansion in EM algorithm.
#' @param maxIter_ICM is the maximum iteration of ICM algorithm.
#' @param maxIter is the maximum iteration of EM algorithm.
#' @examples
#' ## read the simulated data 
#' data(PC)
#' data(sce)
#' ## find the neighborhood of spot
#' platform = "ST"
#' Adj <- find_neighbors2(sce, platform)
#' ## set the initial value
#' fit_int = mclust::Mclust(PC, G=2, verbose = FALSE)
#' mu_int <- unname(fit_int$parameter$mean)
#' sigma_int <- unname(fit_int$parameter$variance$sigma)
#' x_int <- fit_int$classification
#' alpha <- -log(fit_int$parameter$pro) * 0
#' bet = 1
#' maxIter_ICM = 10
#' maxIter = 50
#' ## run the gmrfICMEM function
#' out = gmrfICMEM(PC, x_int, Adj, mu_int, sigma_int, alpha, beta = bet, PX = TRUE, 
#' maxIter_ICM = maxIter_ICM, maxIter = maxIter)
#' @return a list.
#' 
#' The item 'x' is the clustering result. 
#' 
#' The item 'gam' is the posterior probability matrix.
#'  
#' The item 'ell' is the opposite log-likelihood. 
#' 
#' The item 'mu' is the mean of each component.
#' 
#' The item 'sigma' is the variance of each component.
#' 
#' @export 
gmrfICMEM <- function(y, x_int, Adj, mu_int, sigma_int, alpha, beta, PX, maxIter_ICM, maxIter) {
    .Call(`_SC_MEB_gmrfICMEM`, y, x_int, Adj, mu_int, sigma_int, alpha, beta, PX, maxIter_ICM, maxIter)
}

#' find_neighbors2.
#' 
#' @description
#' find_neighbors2 was used to find the neighborhood of spot. 
#'
#' @details find_neighbors2 was used to find the neighborhood of spot. 
#' @param sce is a SingleCellExperiment object containing PCA and position informatin. 
#' @param platform is the name of spatial transcriptomic platform. Specify 'Visium' for hex lattice geometry or 'ST' for square lattice geometry. Specifying this parameter is optional as this information is included in their metadata. 
#' @return a matrix recording the information of neighborhood. 
#' @examples
#' ## read the simulated data 
#' data(sce)
#' platform = "ST"
#' out = find_neighbors2(sce, platform)
#' @export
find_neighbors2 <- function(sce, platform) {
  if (platform == "Visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
  } else if (platform == "ST") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
  } else {
    stop(".find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  ## Get array coordinates (and label by index of spot in SCE)
  spot.positions <- SingleCellExperiment::colData(sce)[, c("col", "row")]
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  
  ## Compute coordinates of each possible spot neighbor
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
  
  ## Select spots that exist at neighbor coordinates
  neighbors <- merge(as.data.frame(neighbor.positions), 
                     as.data.frame(spot.positions), 
                     by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                     suffixes=c(".primary", ".neighbor"),
                     all.x=TRUE)
  
  ## Shift to zero-indexing for C++
  #neighbors$spot.idx.neighbor <- neighbors$spot.idx.neighbor - 1
  
  ## Group neighbor indices by spot 
  ## (sort first for consistency with older implementation)
  neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)
  
  ## Discard neighboring spots without spot data
  ## This can be implemented by eliminating `all.x=TRUE` above, but
  ## this makes it easier to keep empty lists for spots with no neighbors
  ## (as expected by C++ code)
  df_j <- purrr::map(df_j, function(nbrs) purrr::discard(nbrs, function(x) is.na(x)))
  
  ## Log number of spots with neighbors
  n_with_neighbors <- length(purrr::keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ",
          ncol(sce), " spots.")
  
  n <- length(df_j) 
  ##D <- matrix(0,  nrow = n, ncol = n)
  D <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 0)
  for (i in 1:n) {
    if(length(df_j[[i]]) != 0)
      D[i, df_j[[i]]] <- 1
  }
  D
}

