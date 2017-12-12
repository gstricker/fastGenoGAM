#define ARMA_USE_SUPERLU 1

#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::sp_mat mult_sp(arma::sp_mat X, arma::sp_mat Y) {
  return X * Y;
}

// [[Rcpp::export]]
arma::sp_mat trans_(arma::sp_mat X) {
  return trans(X);
}

// [[Rcpp::export]]
arma::mat mult_spmat(arma::sp_mat X, arma::mat Y){
    return X * Y;
}

// [[Rcpp::export]]
double ll_pen_nb(arma::vec beta, arma::sp_mat X, arma::vec y, arma::vec offset,
		double theta, double lambda, arma::sp_mat S, double ll_factor,
		double lambda_factor, int n) {
  // theta must be a matrix to use lgamma on it
  arma::mat theta_mat(1,1);
  theta_mat.fill(theta);

  // mu = exp(eta) = exp(offset + X * beta)
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta);

  // some other parts
  arma::vec aux1 = theta + y;
  arma::vec aux2 = theta + mu;
  arma::vec lga = lgamma(aux1);
  arma::mat lgt = lgamma(theta_mat);
  double lgts = as_scalar(lgt);
  arma::vec lf = lgamma(y + 1);

  // first part of equation with the gamma and factorial fucntion
  double ls = sum(lga - (lf + lgts));

  // second part of the equation
  arma::vec rs = (trans(y) * eta) + (n * theta * log(theta)) - (trans(aux1) * log(aux2));
  double rss = as_scalar(rs);

  // penalization term
  arma::mat pen = (trans(beta) * S) * beta;
  double pens = as_scalar(pen);

  // all together normalized by a factor
  double res = (-1/ll_factor) * (ls + rss) + ((1/lambda_factor) * (lambda * pens));
  return res;
}

// [[Rcpp::export]]
arma::vec gr_ll_pen_nb(arma::vec beta, arma::sp_mat X, arma::vec y,
		       arma::vec offset, double theta,
		       double lambda, arma::sp_mat S) {
  
  // mu = exp(eta) = exp(offset + X * beta)
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta);
  arma::vec z = (y - mu)/(1 + (mu/theta));
  arma::vec pen = S * beta;
  arma::vec gr = trans(X) * z;
  arma::vec res = ((-1) * gr) + (2 * lambda * pen);
  return res;
}

// [[Rcpp::export]]
arma::vec gr_ll_pen_nb_alt(arma::vec beta, arma::sp_mat X, arma::sp_mat XT, arma::vec y,
			   arma::vec offset, double theta,
			   double lambda, arma::sp_mat S) {
  
  // mu = exp(eta) = exp(offset + X * beta)
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta);
  arma::vec z = (y - mu)/(1 + (mu/theta));
  arma::vec pen = S * beta;
  arma::vec gr = XT *z;
  arma::vec res = ((-1) * gr) + (2 * lambda * pen);
  return res;
}

// [[Rcpp::export]]
arma::vec gr_ll_pen_nb_mat(arma::mat beta, arma::sp_mat X, arma::mat y,
		       arma::mat offset, double theta,
		       double lambda, arma::sp_mat S) {
  
  // mu = exp(eta) = exp(offset + X * beta)
  arma::mat lin = X * beta;
  arma::mat eta = offset + lin;
  arma::mat mu = exp(eta);
  arma::mat z = (y - mu)/(1 + (mu/theta));
  arma::mat pen = S * beta;
  arma::mat gr = trans(X) * z;
  arma::mat res = ((-1) * gr) + (2 * lambda * pen);
  return res;
}

// [[Rcpp::export]]
arma::vec gr_ll_pen_nb_mat_alt(arma::mat beta, arma::sp_mat X, arma::sp_mat XT, arma::mat y,
			   arma::mat offset, double theta,
			   double lambda, arma::sp_mat S) {
  
  // mu = exp(eta) = exp(offset + X * beta)
  arma::mat lin = vectorise(X * beta);
  arma::mat eta = offset + lin;
  arma::mat mu = exp(eta);
  arma::mat z = (y - mu)/(1 + (mu/theta));
  arma::mat pen = S * beta;
  arma::mat gr = XT * z;
  arma::mat res = ((-1) * gr) + (2 * lambda * pen);
  return res;
}

// [[Rcpp::export]]
arma::vec hessian_nb_d(arma::vec eta, arma::vec y, double theta) {
  
  arma::vec mu = exp(eta);
  arma::vec ls = mu % (1 + y/theta);
  arma::vec rs = square(1 + mu/theta);
  arma::vec res = ls/rs;
  return res;
}

// [[Rcpp::export]]
arma::sp_mat compute_hessian_nb(arma::sp_mat XT, arma::sp_mat X, arma::sp_mat D,
				 arma::sp_mat S, double lambda) {
  arma::sp_mat H = XT * D * X;
  arma::sp_mat pen = 2 * lambda * S;
  arma::sp_mat res = H + pen;
  return res;
}

// [[Rcpp::export]]
arma::sp_mat compute_pen_hessian(arma::vec beta, arma::sp_mat X, arma::sp_mat XT, arma::vec offset,
				 arma::vec y, arma::sp_mat S, double lambda, double theta) {
  arma::vec lin = vectorise(X * beta);
  arma::vec eta = offset + lin;
  arma::vec mu = exp(eta);

  arma::vec ls = mu % (1 + y/theta);
  arma::vec rs = square(1 + mu/theta);
  arma::vec dd = ls/rs;

  arma::sp_mat D = diagmat(sp_mat(dd));

  arma::sp_mat H = XT * D * X;
  arma::sp_mat pen = 2 * lambda * S;
  arma::sp_mat res = H + pen;
  return res;
}

// [[Rcpp::export]]
arma::sp_mat standard_error(arma::sp_mat X, arma::sp_mat H) {
  arma::sp_mat V = (X * H) % X;
  arma::sp_mat Diag = sum(V, 1);
  arma::sp_mat res = sqrt(Diag);
  return res;
}

// [[Rcpp::export]]
arma::vec solve_hessian(arma::sp_mat H, arma::vec e) {
  arma::vec res = spsolve(H, e);
  return res;
}
