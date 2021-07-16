
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
 
using namespace std;
using namespace Rcpp;
using namespace arma;


const double log2pi = std::log(2.0 * M_PI);


////////////////////////////////////////////////////////////////////////////////
arma::mat getPairDist(const arma::mat x)	{
	int N = x.n_rows;
	arma::mat D(N, N);
	for (int j = 0; j < N; ++j)
	{
		for (int i = j; i < N; ++i)
		{
			D(i, j)	= norm(x.row(i) - x.row(j), 2);
			D(j, i) = D(i, j);
		}
	}
	return D;
}

double elbo(arma::mat& U, arma::mat& gam)	{
   return accu(U % gam) + accu(gam % log(gam + (gam==0)));
}


arma::mat Smooth(int K, double beta)	{
	arma::mat S = beta * (ones<arma::mat>(K, K) - eye<arma::mat>(K, K));
	return S;
}

vec dmvnrm(const mat& x,  
           rowvec mean,  
           mat sigma, 
           bool logd = false,
           int cores = 1) { 
    omp_set_num_threads(cores);
    int n = x.n_rows;
    int xdim = x.n_cols;
    vec out(n);
    mat rooti = inv(chol(sigma, "lower"));
    double rootisum = sum(log(rooti.diag()));
    double constants = -(xdim/2) * log2pi;
    #pragma omp parallel for schedule(static) 
    for (int i=0; i < n; i++) {
        vec z = rooti * trans( x.row(i) - mean) ;    
        out(i)      = constants - 0.5 * sum(z%z) + rootisum;     
    }  
      
    if (logd==false) {
        out=exp(out);
    }
    return(out);
}



sp_mat get_spNbs(ivec x, const sp_mat& Adj) {   
// row is for pixel.
//output a sparse matrix, i-th row contains labels of neighbor_i. 
	// Make const iterator
    arma::sp_mat::const_iterator start = Adj.begin();
    //arma::sp_mat::const_iterator end   = Adj.end();

    // Calculate number of nonzero points
    //int n = std::distance(start, end);
    int n = Adj.n_nonzero;

	sp_mat spNbs(x.n_elem, x.n_elem);    //neiborhood state matrix, matched with Adj.

	
	arma::sp_mat::const_iterator it = start; 
    for(int i = 0; i < n; ++i)
    {
        //temp(0) = it.row();
        //temp(1) = it.col();
        spNbs(it.row(), it.col()) = x(it.col());
        ++it; // increment
    }

	return spNbs;
}

arma::rowvec calXenergy2D_i(arma::ivec x, int i, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
	arma::sp_mat spNbs = get_spNbs(x, Adj);
	arma::sp_mat spNbs_t = spNbs.t();  // transform spNbs to iterate by column.
	arma::rowvec Ux_i = zeros<arma::rowvec>(K);
	int k;
	for (k = 0; k < K; k++)
	{
		arma::sp_mat col(spNbs_t.col(i));
		double n_sameS = 0;
		int nn = col.n_nonzero;
		if (nn == 0)	{
			Ux_i(k) = alpha(k);
		} else {
			for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
				n_sameS += (*j) == (k+1);
			}
			Ux_i(k) = alpha(k) + beta * (nn - n_sameS)/2;
			//Ux_i(k) = alpha(k) + 2 * beta * (1 - n_sameS/nn);
		}
	}
	return Ux_i;
}
 


arma::mat calXenergy2D_sp(arma::ivec x, const arma::sp_mat& Adj, int K, const arma::vec alpha, const double beta)	{
	int n = x.n_rows;
	arma::sp_mat spNbs = get_spNbs(x, Adj);
	arma::sp_mat spNbs_t = spNbs.t();  // transform spNbs to iterate by column.
	arma::mat Ux(n, K);
	int i, k;
	for (k = 0; k < K; k++)
	{
		for (i = 0; i < n; i++)
		{
			arma::sp_mat col(spNbs_t.col(i));
			double n_sameS = 0;
			//double n_nb = 0;
			//int nn = std::distance(col.begin(), col.end());
			int nn = col.n_nonzero;
			for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
				n_sameS += (*j) == (k+1);
				//n_nb++;
			}
			Ux(i, k) = alpha(k) + beta * (nn - n_sameS)/2;
		}
	}
	return Ux;

}


 
Rcpp::List runICM_sp (const arma::mat &y, arma::ivec x, arma::mat mu, arma::cube sigma, const arma::sp_mat& Adj,
			arma::vec alpha, double beta, int maxIter_ICM)	{
	int n = y.n_rows, K = mu.n_cols;
	int iter, k;

	// energy of y
	arma::mat Uy(n, K);
	for (k = 0; k < K; k++)	{
		arma::vec mu_k = mu.col(k);
		arma::mat sigma_k = sigma.slice(k);
		Uy.col(k) = -dmvnrm(y, conv_to< rowvec >::from(mu_k), sigma_k, true);
	}

	arma::vec Energy(maxIter_ICM);
   	Energy(0) = INFINITY;
   	arma::mat Ux(n, K);
   	arma::mat U(n, K);
	//--------------------------------------------------------------------------------	
	// ICM algrithm
	//--------------------------------------------------------------------------------
	int Iteration = 1;
	for (iter = 1; iter < maxIter_ICM; iter ++ ) {
		
		Ux = calXenergy2D_sp(x, Adj, K, alpha, beta);

		U = Uy + Ux;
		arma::vec Umin = min(U, 1);
		arma::uvec x_u = index_min(U, 1);
		x = conv_to< ivec >::from(x_u) + 1;

		Energy(iter) = sum(Umin);
		if (Energy(iter) - Energy(iter - 1) > 1e-5) {
			Rprintf("diff Energy = %f\n", Energy(iter) - Energy(iter - 1));
			break;
		}

		if (Energy(iter-1) - Energy(iter) < 1e-5)
		{
			Rprintf("ICM Converged at Iteration = %d\n", iter);
			break;
		}
	}
	
	if (iter == maxIter_ICM) {
		Iteration = iter - 1;
	} else {
		Iteration = iter;
	}

	vec energy = Energy.subvec(1, Iteration);
	
	arma::mat pxgn = normalise(exp(-Ux), 1, 1); // set all rowSums to be ONE.

	List output = List::create(
		Rcpp::Named("x") = x,
		Rcpp::Named("U") = U,
		Rcpp::Named("Uy") = Uy,
		Rcpp::Named("Ux") = Ux,
		Rcpp::Named("pxgn") = pxgn,
		Rcpp::Named("energy") = energy);

	return output; 

}
 


vec ud_mu(const mat& y, vec gam_k)	{
	double N_k = sum(gam_k);
	int n = y.n_rows, p = y.n_cols;

	vec mu_k = zeros<vec>(p);
	for (int j = 0; j < p; j ++)	{
		double temp = 0;
		for (int i = 0; i < n; i ++)	{
		 	temp += gam_k(i) * y(i, j);
		}
		mu_k(j) = temp/N_k;
	}
	return mu_k;
} 


mat ud_sigma(const mat& y, vec mu_k, vec gam_k)	{
	double N_k = sum(gam_k);
	int n = y.n_rows, p = y.n_cols;

	mat sigma = eye<mat>(p, p);
	for (int j = 0; j < p; j++)
	{
		for (int k = j; k < p; k++)	{
			double temp = 0;
			for (int i = 0; i < n; i++)
			{
				temp += gam_k(i) * (y(i, j) - mu_k(j)) * (y(i, k) - mu_k(k));
			}
			sigma(j, k) = temp/N_k;
			sigma(k, j) = sigma(j, k);
		}
		sigma(j, j) += 1e-6;
	}
	return sigma;
}
 

mat runEstep(mat Uy, mat pxgn)	{
	int n = Uy.n_rows, K = Uy.n_cols;
	mat gam = zeros<mat>(n, K); 
	mat dy = exp(-Uy);
	mat gam1 = zeros<mat>(n, K); 
	int k;
	// compute pygn
	for (k = 0; k < K; k ++)	{
		gam1.col(k) = dy.col(k) % pxgn.col(k);
	}
	
	vec pygn = sum(gam1, 1);

	// compute gam
	for (k = 0; k < K; k ++)	{
		gam.col(k) = gam1.col(k)/pygn;
	}
	return gam;
}


 
// [[Rcpp::export]]
List	gmrfICMEM(const arma::mat& y, const arma::ivec& x_int, const arma::sp_mat& Adj,
			   const arma::mat& mu_int, const arma::cube& sigma_int, 
			   const arma::vec& alpha, 	const double&  beta,
	 	  	   const bool& PX,          const int& maxIter_ICM, const int& maxIter)
{
	int n = y.n_rows, K = mu_int.n_cols; //, p = y.n_cols;
	mat mu = mu_int;
	cube sigma = sigma_int;
	ivec x = x_int;		  // label	
	mat pxgn(n, K); // p(x_i | x_{N_i})
	pxgn.fill(1.0/K);  
	vec pygn = ones<vec>(n);		  // p(y_i | x_{N_i})
	mat gam = zeros<mat>(n, K);   // responsibility
	mat pygx = zeros<mat>(n, K);  
	//mat gam1 = zeros<mat>(n, K);  // p(y_i | x_i = k) * p(x_i = k | x_{N_i})

	// Parameter Expansion; double lam = 1.0;
	double ell = 0;
	vec LogLik(maxIter);
   	LogLik(0) = INFINITY;
   	int iter, k;
	//--------------------------------------------------------------------------------	
	// EM algrithm
	//--------------------------------------------------------------------------------
	int Iteration = 1;
	for (iter = 1; iter < maxIter; iter ++ ) {

		// ICM, List;
		// update x and pxgn
		List fitICM = runICM_sp(y, x, mu, sigma, Adj, alpha, beta, maxIter_ICM);
		ivec xHat = fitICM["x"]; 
		x = xHat;
		mat U = fitICM["U"]; 
		mat Uy = fitICM["Uy"]; 
		pygx = Uy; // save for latter use, will delete later.
		mat pxgnHat = fitICM["pxgn"];
		pxgn = pxgnHat;
		vec energy = fitICM["energy"];
		LogLik(iter) = min(energy);
		// E-step, update gamma.
		
		gam = runEstep(Uy, pxgn);

		// M-step 
		for (k = 0; k < K; k ++)	{
			vec gam_k = gam.col(k);
			mu.col(k) = ud_mu(y, gam_k);
			vec mu_k = mu.col(k);
			sigma.slice(k) = ud_sigma(y, mu_k, gam_k);
		}
		
        mat U2 = -log(pxgn) + Uy;
		ell = elbo(U2, gam);
        
		if ( LogLik(iter) - LogLik(iter - 1) > 1e-5 ){
			//perror("The energy failed to decrease!");
			break;
		}

		//if (abs(LogLik(iter) - LogLik(iter - 1)) < 1e-5 || rcond(sigma) < 1e-7) {
		if (abs(LogLik(iter) - LogLik(iter - 1)) < 1e-5) {
Rprintf("Converged at Iteration = %d\n", iter);
			break;
		}
	}
	

	if (iter == maxIter) {
		Iteration = iter - 1;
	} else {
		Iteration = iter;
	}

	vec loglik;
	loglik = LogLik.subvec(1, Iteration);
	// Reduce-step
	

	List output = List::create(
		Rcpp::Named("x") = x,
		Rcpp::Named("gam") = gam,
		Rcpp::Named("pxgn") = pxgn,
		Rcpp::Named("pygx") = pygx,
		Rcpp::Named("mu") = mu,
		Rcpp::Named("sigma") = sigma,
		Rcpp::Named("ell") = ell,
		Rcpp::Named("loglik") = loglik);

	return output; 

}


 

 


