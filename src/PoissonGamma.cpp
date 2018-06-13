#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_beta_psi_mu2(arma::vec Y, arma::vec X, double beta, double psi, double mu2, int n){

  arma::vec gradient = arma::zeros<arma::vec>(3);

  double common_terms = log(psi) - R::digamma(psi);
  for(int i = 0; i < n; ++i){

    double e_term = exp(X(i)*beta + mu2);
    double e_term_psi = e_term + psi;
    double Y_i_psi = Y(i) + psi;
    double e_term_fac = Y_i_psi/e_term_psi;
    double f_term = Y(i) - e_term*e_term_fac;
    gradient(0) += f_term*X(i);

    double gradient_psi_i = - log(e_term_psi) - e_term_fac + R::digamma(Y_i_psi);
    gradient(1) += gradient_psi_i;

    gradient(2) += f_term; //mu2

  }
  gradient(1) += common_terms*n + n;

  return gradient;
}

// [[Rcpp::export]]

double log_factorial(int Y){
  double res = 0;
  for(int kk = 1; kk <= Y; ++kk){
    res += log(kk);
  }
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec log_factorial_calculated(int N){

  arma::vec values = arma::zeros<arma::vec>(N+1);

  for(int kk = 1; kk <= N; ++kk){
    values(kk) = values(kk-1) + log(kk);
  }

  return values;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double LogLikelihood_beta_psi_mu2(arma::vec Y, arma::vec X, double beta, double psi, double mu2, int n, double sum_log_factorial_Y){

  double Likelihood = 0;
  double common_term = -lgamma(psi) + psi*log(psi);

  for(int i = 0; i < n; ++i){
    double Xbeta_mu2 = X(i)*beta + mu2;
    double e_term = exp(Xbeta_mu2);
    double psi_Yi = psi+Y(i);
    double ll = lgamma(psi_Yi) - psi_Yi*log(psi + e_term) + Y(i)*(Xbeta_mu2);
    Likelihood += ll;
  }

  Likelihood += common_term*n - sum_log_factorial_Y;
  return Likelihood;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_beta(arma::vec Y, arma::vec X, double gra_beta, double ll, double beta, double psi, double mu2, int n, double gamma, double sum_log_factorial_Y, double down){

  double gra_beta2 = gra_beta*gra_beta*gamma;
//  double lb = LogLikelihood_beta_psi(Y, X, beta, psi, n);
  double start = 0.1;
  if(beta >= 0.001) start = sqrt(abs(beta/gra_beta))/2;

  double aa = start;
  double selected = beta;
  while(aa > 0){
    double aa2 = aa*aa;
    double beta_prime = beta + aa2*gra_beta;
    double lb_prime = LogLikelihood_beta_psi_mu2(Y, X, beta_prime, psi, mu2, n, sum_log_factorial_Y);
    if(lb_prime - ll - aa2*gra_beta2 > 0 ) { // abs(lb_prime - ll - aa2*gra_beta2) < 0.0001) {
      selected = beta_prime;
      break;
    }
    aa = aa - start*down;
  }

  return selected;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_psi(arma::vec Y, arma::vec X, double gra_psi, double ll, double beta, double psi, double mu2, int n, double gamma, double sum_log_factorial_Y, double down, double psi_cutoff){

  double gra_psi2 = gra_psi*gra_psi*gamma;
  double start = sqrt(abs(psi/gra_psi))/2;

  double aa = start;
  double selected = psi;
  while(aa > 0){
    double aa2 = aa*aa;
    double psi_prime = psi + aa2*gra_psi;
    if(psi_prime > 0){
      double lpsi_prime = LogLikelihood_beta_psi_mu2(Y, X, beta, psi_prime, mu2, n, sum_log_factorial_Y);
      if(lpsi_prime - ll - aa2*gra_psi2 > 0 ) { // abs(lpsi_prime - ll - aa2*gra_psi2) < 0.0001) {
        selected = psi_prime;
        break;
      }
    }
    aa = aa - start*down;
  }
  if(selected >= psi_cutoff) selected = psi_cutoff;
  return selected;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_mu2(arma::vec Y, arma::vec X, double gra_mu2, double ll, double beta, double psi, double mu2, int n, double gamma, double sum_log_factorial_Y, double down){

  double gra_mu2_2 = gra_mu2*gra_mu2*gamma;
  double start = sqrt(abs(mu2/gra_mu2))/2;

  double aa = start;
  double selected = mu2;
  while(aa > 0){
    double aa2 = aa*aa;
    double mu2_prime = mu2 + aa2*gra_mu2;
    double lmu2_prime = LogLikelihood_beta_psi_mu2(Y, X, beta, psi, mu2_prime, n, sum_log_factorial_Y);
    if(lmu2_prime - ll - aa2*gra_mu2_2 > 0 ) { // abs(lmu2_prime - ll - aa2*gra_mu2_2) < 0.0001) {
      selected = mu2_prime;
      break;
    }
    aa = aa - start*down;
  }

  return selected;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_descent_beta_psi_mu2(arma::vec Y, arma::vec X, double beta, double psi, double mu2, double gamma, int steps, double sum_log_factorial_Y, double down, double psi_cutoff){

  int n = Y.n_elem;
  arma::vec gradient = gradient_beta_psi_mu2(Y, X, beta, psi, mu2, n);
  double ll = LogLikelihood_beta_psi_mu2(Y, X, beta, psi, mu2, n, sum_log_factorial_Y);

  double beta_prime = 0; double psi_prime = 0; double mu2_prime = 0;
  for(int i = 0; i < steps; ++i){
    if(abs(gradient(0)) >= 0.0001){
      beta_prime = select_stepsize_for_beta(Y, X, gradient(0), ll, beta, psi, mu2, n, gamma, sum_log_factorial_Y, down);
    } else {
      beta_prime = beta;
    }
    if(abs(gradient(1)) >= 0.0001){
      psi_prime = select_stepsize_for_psi(Y, X, gradient(1), ll, beta, psi, mu2, n, gamma, sum_log_factorial_Y, down, psi_cutoff);
    } else {
      psi_prime = psi;
    }
    if(abs(gradient(2)) >= 0.0001){
      mu2_prime = select_stepsize_for_mu2(Y, X, gradient(2), ll, beta, psi, mu2, n, gamma, sum_log_factorial_Y, down);
    } else {
      mu2_prime = mu2;
    }
    beta = beta_prime; psi = psi_prime; mu2 = mu2_prime;
    gradient = gradient_beta_psi_mu2(Y, X, beta, psi, mu2, n);
    ll = LogLikelihood_beta_psi_mu2(Y, X, beta, psi, mu2, n, sum_log_factorial_Y);
  }

  arma::vec est = arma::zeros<arma::vec>(4);
  est(0) = beta;
  est(1) = psi;
  est(2) = mu2;
  est(3) = ll;

  return est;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_descent_alt(arma::vec Y, arma::vec X, double psi, double mu2, double gamma, int steps, double sum_log_factorial_Y, double down, double psi_cutoff){

  int n = Y.n_elem;
  arma::vec gradient = gradient_beta_psi_mu2(Y, X, 0, psi, mu2, n);
  double ll = LogLikelihood_beta_psi_mu2(Y, X, 0, psi, mu2, n, sum_log_factorial_Y);
  double psi_prime = 0; double mu2_prime = 0;

  for(int i = 0; i < steps; ++i){
    if(abs(gradient(1)) >= 0.0001){
      psi_prime = select_stepsize_for_psi(Y, X, gradient(1), ll, 0, psi, mu2, n, gamma, sum_log_factorial_Y, down, psi_cutoff);
    }  else {
      psi_prime = psi;
    }
    if(abs(gradient(2)) >= 0.0001){
      mu2_prime = select_stepsize_for_mu2(Y, X, gradient(2), ll, 0, psi, mu2, n, gamma, sum_log_factorial_Y, down);
    } else {
      mu2_prime = mu2;
    }
    psi = psi_prime; mu2 = mu2_prime;
    gradient = gradient_beta_psi_mu2(Y, X, 0, psi, mu2, n);
    ll = LogLikelihood_beta_psi_mu2(Y, X, 0, psi, mu2, n, sum_log_factorial_Y);
  }

  arma::vec est = arma::zeros<arma::vec>(4);
  est(1) = psi;
  est(2) = mu2;
  est(3) = ll;

  return est;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PoissionGamma(arma::vec Y, arma::vec X, double beta, double psi, double mu2, double gamma, int steps, double down, double psi_cutoff){

  arma::vec calculated_values = log_factorial_calculated(Y.max());

  int n = Y.n_elem;
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  double sum_log_factorial_Y = sum(log_factorial_Y);
  arma::vec est = gradient_descent_beta_psi_mu2(Y, X, beta, psi, mu2, gamma, steps, sum_log_factorial_Y, down, psi_cutoff);
  arma::vec est_alt = gradient_descent_alt(Y, X, psi, mu2, gamma, steps, sum_log_factorial_Y, down, psi_cutoff);

  double test_statistics = est(3) - est_alt(3);
  double p_value = 1 - R::pchisq(2*test_statistics, 1, TRUE, FALSE);

  return Rcpp::List::create(Rcpp::Named("beta") = est(0), Rcpp::Named("p_value") = p_value,
                            Rcpp::Named("psi") = est(1), Rcpp::Named("alt_psi") = est_alt(1),
                            Rcpp::Named("mu2") = est(2), Rcpp::Named("alt_mu2") = est_alt(2),
                            Rcpp::Named("likelihood") = est(3), Rcpp::Named("alt_likelihood") = est_alt(3));
}





// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat Fisher_information_one_beta(arma::vec Y, arma::vec X, double beta, double psi, double mu2, int n){

  arma::mat Fisher = arma::zeros<arma::mat>(3,  3);

  arma::vec e_term = exp(X*beta + mu2);
  arma::vec e_term_fac_prime = 1/(e_term + psi);
  arma::vec e_term_fac = e_term%e_term_fac_prime;
  arma::vec e_term_fac_2 = 1/pow(e_term + psi, 2);
  arma::vec e_term_square_inverse = e_term%e_term_fac_2;
  arma::vec Y_i_psi_e_term_psi_inv2 = (Y + psi)%e_term_square_inverse;
  arma::vec X_square = pow(X, 2);

  arma::vec psi2_vec = arma::zeros<arma::vec>(n);
  for(int j = 0; j < n; ++j){
    psi2_vec(j) = - R::trigamma(Y(j) + psi) - 1/psi - (Y(j) + psi)*e_term_fac_2(j) + R::trigamma(psi) + 2*e_term_fac_prime(j);
  }

  Fisher(0, 0) = sum(psi2_vec); //psi, psi
  Fisher(0, 1) = Fisher(1, 0) = sum(e_term_fac) - sum(Y_i_psi_e_term_psi_inv2); //psi, mu2
  Fisher(0, 2) = Fisher(2, 0) = sum(e_term_fac%X) - sum(X%Y_i_psi_e_term_psi_inv2); //psi, beta

  Fisher(1, 1) = sum(psi*Y_i_psi_e_term_psi_inv2); //mu2, mu2
  Fisher(1, 2) = Fisher(2, 1) = sum(psi*X%Y_i_psi_e_term_psi_inv2); //mu2, beta

  Fisher(2, 2) = sum(psi*Y_i_psi_e_term_psi_inv2%X_square); //beta, beta

  return Fisher;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PoissionGamma_FISHER(arma::vec Y, arma::vec X, double beta, double psi, double mu2, double gamma, int steps, double down, double psi_cutoff){

  arma::vec calculated_values = log_factorial_calculated(Y.max());

  int n = Y.n_elem;
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  double sum_log_factorial_Y = sum(log_factorial_Y);
  arma::vec est = gradient_descent_beta_psi_mu2(Y, X, beta, psi, mu2, gamma, steps, sum_log_factorial_Y, down, psi_cutoff);

  psi = est(1);
  mu2 = est(2);
  beta = est(0);

  arma::mat Fisher = Fisher_information_one_beta(Y, X, beta, psi, mu2, n);

  arma::mat est_cov = arma::inv(Fisher);

  arma::vec all_SE = arma::zeros<arma::vec>(3);
  arma::vec all_pvalue = arma::zeros<arma::vec>(3);

  for(int j = 0; j < 3; ++j){
    double SE = sqrt(abs(est_cov(j, j)));
    all_SE(j) = SE;
    if(j == 0){
      all_pvalue(0) = 2*R::pnorm(-abs(psi)/SE, 0, 1, TRUE, FALSE);
    }
    if(j == 1){
      all_pvalue(1) = 2*R::pnorm(-abs(mu2)/SE, 0, 1, TRUE, FALSE);
    }
    if(j == 2){
      all_pvalue(2) = 2*R::pnorm(-abs(beta)/SE, 0, 1, TRUE, FALSE);
    }
  }


  return Rcpp::List::create(Rcpp::Named("psi") = psi,
                            Rcpp::Named("mu2") = mu2,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("SE") = all_SE,
                            Rcpp::Named("p_value") = all_pvalue,
                            Rcpp::Named("likelihood") = est(3));

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List gradient_and_LogLikelihood_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, arma::vec WY, arma::vec WWY, arma::vec W3Y, arma::vec W4Y,
                                         arma::vec VY, arma::vec VVY, arma::vec V3Y, arma::vec V4Y,
                                         arma::vec WW, arma::vec W3, arma::vec W4, arma::vec VV, arma::vec V3, arma::vec V4,
                                         double a0, double a1, double a2, double a3, double a4,
                                         double b1, double b2, double b3, double b4,  int n, double sum_log_factorial_Y){

  arma::vec gradient = arma::zeros<arma::vec>(9);
  double Likelihood = 0;

  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    double exp_Ai = -exp(A_i);
    double ll =  exp_Ai*B_i;
    Likelihood += ll + Y(i)*(A_i + log(B_i));

    gradient(0) += ll + Y(i);
    gradient(1) += ll*W(i) + WY(i);
    gradient(2) += ll*WW(i) + WWY(i);
    gradient(3) += ll*W3(i) + W3Y(i);
    gradient(4) += ll*W4(i) + W4Y(i);

    double term2 = 1/B_i;
    gradient(5) += exp_Ai*V(i) + term2*VY(i);
    gradient(6) += exp_Ai*VV(i) + term2*VVY(i);
    gradient(7) += exp_Ai*V3(i) + term2*V3Y(i);
    gradient(8) += exp_Ai*V4(i) + term2*V4Y(i);

  }

  Likelihood -= sum_log_factorial_Y;
  return Rcpp::List::create(Rcpp::Named("gradient") = gradient,
                            Rcpp::Named("Likelihood") = Likelihood);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]


double LogLikelihood_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, arma::vec WW, arma::vec VV,
                                           arma::vec W3, arma::vec V3, arma::vec W4, arma::vec V4,
                                           double a0, double a1, double a2, double a3, double a4,
                                           double b1, double b2, double b3, double b4,
                                           int n, double sum_log_factorial_Y){

  double Likelihood = 0;

  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    double exp_Ai = -exp(A_i);
    double ll =  exp_Ai*B_i;
    Likelihood += ll + Y(i)*(A_i + log(B_i));
  }

  Likelihood -= sum_log_factorial_Y;
  return Likelihood;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double select_stepsize_for_a_parameter(arma::vec Y, arma::vec W, arma::vec V, arma::vec WW, arma::vec VV,
                                       arma::vec W3, arma::vec V3, arma::vec W4, arma::vec V4, double ll, double sum_log_factorial_Y,
                                       arma::vec gradient, arma::vec parameters, int ind, double gamma, int n, double down){

  //if(ind < 5) gamma = 0.75;
  double gra = gradient(ind);
  double gra_2 = gra*gra*gamma;
  double para = parameters(ind);
  double start = sqrt(abs(para/gra))/5;
  //if(ind < 5) start = start/10;
  double a0 = parameters(0);
  double a1 = parameters(1);
  double a2 = parameters(2);
  double a3 = parameters(3);
  double a4 = parameters(4);
  double b1 = parameters(5);
  double b2 = parameters(6);
  double b3 = parameters(7);
  double b4 = parameters(8);

  double aa = start;
  double selected = para;
  double ll_prime = ll;
  while(aa > 0){
    double aa2 = aa*aa;
    double para_prime = para + aa2*gra;
    if(ind == 0){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, para_prime, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(ind == 1){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, para_prime, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(ind == 2){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, para_prime, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(ind == 3){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, para_prime, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(ind == 4){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, para_prime, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(ind == 5){
    //  Rcpp::Rcout << "para" << para << std::endl;
    //  Rcpp::Rcout << "gra" << gra  << std::endl;
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, para_prime, b2, b3, b4, n, sum_log_factorial_Y);
    //  Rcpp::Rcout << "dis" << ll_prime - ll - aa2*gra_2  << std::endl;
    //  Rcpp::Rcout << "aa" << aa << std::endl;
    //  Rcpp::Rcout << "para_prime" << para_prime << std::endl;
    //  Rcpp::Rcout << "ll_prime" << ll_prime << std::endl;
    }
    if(ind == 6){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, para_prime, b3, b4, n, sum_log_factorial_Y);
    }
    if(ind == 7){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, para_prime, b4, n, sum_log_factorial_Y);
    }
    if(ind == 8){
      ll_prime = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, para_prime, n, sum_log_factorial_Y);
    }

    if(ll_prime - ll - aa2*gra_2 > 0 ) { //| abs(ll_prime - ll - aa2*gra_2) < 0.0001) {
      selected = para_prime;
      break;
    }
    aa = aa - start*down;
  }

  return selected;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List gradient_descent_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, double a0, double a1, double a2, double a3, double a4,
                                                  double b1, double b2, double b3, double b4, double gamma, int steps, double down){

  int n = Y.n_elem;

  arma::vec calculated_values = log_factorial_calculated(Y.max());
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  double sum_log_factorial_Y = sum(log_factorial_Y);

  arma::vec WW = W%W;
  arma::vec W3 = W%WW;
  arma::vec W4 = W%W3;
  arma::vec WY = W%Y;
  arma::vec WWY = W%WY;
  arma::vec W3Y = W%WWY;
  arma::vec W4Y = W%W3Y;

  arma::vec VV = V%V;
  arma::vec V3 = V%VV;
  arma::vec V4 = V%V3;
  arma::vec VY = V%Y;
  arma::vec VVY = V%VY;
  arma::vec V3Y = V%VVY;
  arma::vec V4Y = V%V3Y;
  Rcpp::List res = gradient_and_LogLikelihood_for_individual_sample(Y, W, V, WY, WWY, W3Y, W4Y, VY, VVY, V3Y, V4Y,
                                                      WW, W3, W4, VV, V3, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
  arma::vec gradient = res["gradient"];
  double ll = res["Likelihood"];

  arma::vec parameters = arma::zeros<arma::vec>(9);
  parameters(0) = a0;
  parameters(1) = a1;
  parameters(2) = a2;
  parameters(3) = a3;
  parameters(4) = a4;
  parameters(5) = b1;
  parameters(6) = b2;
  parameters(7) = b3;
  parameters(8) = b4;

  double a0_prime = 0; double a1_prime = 0; double a2_prime = 0; double a3_prime = 0; double a4_prime = 0;
  double b1_prime = 0; double b2_prime = 0; double b3_prime = 0; double b4_prime = 0;

  for(int i = 0; i < steps; ++i){
    if(abs(gradient(0)) >= 0.0001){
      a0_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 0, gamma, n, down);
    } else {
      a0_prime = a0;
    }
    if(abs(gradient(1)) >= 0.0001){
      a1_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 1, gamma, n, down);
    } else {
      a1_prime = a1;
    }
    if(abs(gradient(2)) >= 0.0001){
      a2_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y,  gradient, parameters, 2, gamma, n, down);
    } else {
      a2_prime = a2;
    }
    if(abs(gradient(3)) >= 0.0001){
      a3_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 3, gamma, n, down);
    } else {
      a3_prime = a3;
    }
    if(abs(gradient(4)) >= 0.0001){
      a4_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 4, gamma, n, down);
    } else {
      a4_prime = a4;
    }
    if(abs(gradient(5)) >= 0.0001){
      b1_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 5, gamma, n, down);
    } else {
      b1_prime = b1;
    }
    if(abs(gradient(6)) >= 0.0001){
      b2_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y,  gradient, parameters, 6, gamma, n, down);
    } else {
      b2_prime = b2;
    }
    if(abs(gradient(7)) >= 0.0001){
      b3_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 7, gamma, n, down);
    } else {
      b3_prime = b3;
    }
    if(abs(gradient(8)) >= 0.0001){
      b4_prime = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 8, gamma, n, down);
    } else {
      b4_prime = b4;
    }

    a0 = a0_prime; a1 = a1_prime; a2 = a2_prime; a3 = a3_prime; a4 = a4_prime;
    b1 = b1_prime; b2 = b2_prime; b3 = b3_prime; b4 = b4_prime;

    Rcpp::List res = gradient_and_LogLikelihood_for_individual_sample(Y, W, V, WY, WWY, W3Y, W4Y, VY, VVY, V3Y, V4Y,
                                                                      WW, W3, W4, VV, V3, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    arma::vec gradient_1 = res["gradient"];
    gradient = gradient_1;
    ll = res["Likelihood"];

    parameters(0) = a0;
    parameters(1) = a1;
    parameters(2) = a2;
    parameters(3) = a3;
    parameters(4) = a4;
    parameters(5) = b1;
    parameters(6) = b2;
    parameters(7) = b3;
    parameters(8) = b4;

   // Rcpp::Rcout << gradient << std::endl;
  }

  arma::vec corrected = Y;
  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    corrected(i) = B_i*exp(A_i);
  }

  return Rcpp::List::create(Rcpp::Named("parameters") = parameters,
                          Rcpp::Named("corrected") = corrected);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List coordinate_descent_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, double a0, double a1, double a2, double a3, double a4,
                                                  double b1, double b2, double b3, double b4, double gamma, int steps, double down){

  int n = Y.n_elem;

  arma::vec calculated_values = log_factorial_calculated(Y.max());
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  double sum_log_factorial_Y = sum(log_factorial_Y);

  arma::vec WW = W%W;
  arma::vec W3 = W%WW;
  arma::vec W4 = W%W3;
  arma::vec WY = W%Y;
  arma::vec WWY = W%WY;
  arma::vec W3Y = W%WWY;
  arma::vec W4Y = W%W3Y;

  arma::vec VV = V%V;
  arma::vec V3 = V%VV;
  arma::vec V4 = V%V3;
  arma::vec VY = V%Y;
  arma::vec VVY = V%VY;
  arma::vec V3Y = V%VVY;
  arma::vec V4Y = V%V3Y;
  Rcpp::List res = gradient_and_LogLikelihood_for_individual_sample(Y, W, V, WY, WWY, W3Y, W4Y, VY, VVY, V3Y, V4Y,
                                                                    WW, W3, W4, VV, V3, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
  arma::vec gradient = res["gradient"];
  double ll = res["Likelihood"];

  arma::vec parameters = arma::zeros<arma::vec>(9);
  parameters(0) = a0;
  parameters(1) = a1;
  parameters(2) = a2;
  parameters(3) = a3;
  parameters(4) = a4;
  parameters(5) = b1;
  parameters(6) = b2;
  parameters(7) = b3;
  parameters(8) = b4;

  for(int i = 0; i < steps; ++i){
    if(abs(gradient(0)) >= 0.0001){
      a0 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 0, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(1)) >= 0.0001){
      a1 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 1, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(2)) >= 0.0001){
      a2 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y,  gradient, parameters, 2, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(3)) >= 0.0001){
      a3 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 3, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(4)) >= 0.0001){
      a4 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 4, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(5)) >= 0.0001){
      b1 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 5, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(6)) >= 0.0001){
      b2 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y,  gradient, parameters, 6, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(7)) >= 0.0001){
      b3 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 7, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }
    if(abs(gradient(8)) >= 0.0001){
      b4 = select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 8, gamma, n, down);
      ll = LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    }

    Rcpp::List res = gradient_and_LogLikelihood_for_individual_sample(Y, W, V, WY, WWY, W3Y, W4Y, VY, VVY, V3Y, V4Y,
                                                                      WW, W3, W4, VV, V3, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, n, sum_log_factorial_Y);
    arma::vec gradient_1 = res["gradient"];
    gradient = gradient_1;
    ll = res["Likelihood"];

    parameters(0) = a0;
    parameters(1) = a1;
    parameters(2) = a2;
    parameters(3) = a3;
    parameters(4) = a4;
    parameters(5) = b1;
    parameters(6) = b2;
    parameters(7) = b3;
    parameters(8) = b4;

   // Rcpp::Rcout << gradient << std::endl;
  }

  arma::vec corrected = Y;
  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    corrected(i) = B_i*exp(A_i);
  }

  return Rcpp::List::create(Rcpp::Named("parameters") = parameters,
                            Rcpp::Named("corrected") = corrected);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_multiple_beta_psi_mu2(arma::vec Y, arma::mat X, arma::vec Xbeta, double psi, double mu2, int n, int k){

  arma::vec gradient = arma::zeros<arma::vec>(2 + k);

  double common_terms = log(psi) - R::digamma(psi);

  for(int i = 0; i < n; ++i){

    double e_term = exp(Xbeta(i) + mu2);
    double e_term_psi = e_term + psi;
    double Y_i_psi = Y(i) + psi;
    double e_term_fac = Y_i_psi/e_term_psi;
    double f_term = Y(i) - e_term*e_term_fac;
    double gradient_psi_i = -log(e_term_psi) - e_term_fac + R::digamma(Y_i_psi);
    gradient(0) += gradient_psi_i; //psi
    gradient(1) += f_term; //mu2
    for(int j = 0; j < k; ++j){
      gradient(j + 2) += f_term*X(i, j);
    }

  }
  gradient(0) += common_terms*n + n; //psi


  return gradient;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double LogLikelihood_multiple_beta_psi_mu2(arma::vec Y, arma::vec Xbeta, double psi, double mu2, int n, double sum_log_factorial_Y){

  double Likelihood = 0;
  double common_term = -lgamma(psi) + psi*log(psi);

  for(int i = 0; i < n; ++i){
    double Xbeta_mu2 = Xbeta(i) + mu2;
    double e_term = exp(Xbeta_mu2);
    double psi_Yi = psi+Y(i);
    double ll = lgamma(psi_Yi) - psi_Yi*log(psi + e_term) + Y(i)*(Xbeta_mu2);
    Likelihood += ll;
  }
  Likelihood += common_term*n - sum_log_factorial_Y;
  return Likelihood;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_multiple_beta(arma::vec Y, arma::mat X, arma::vec betas, arma::vec Xbeta, int ind,
                                         arma::vec m_gra_beta, double ll, double psi, double mu2,
                                         int n, double gamma, double sum_log_factorial_Y, double down){

  double gra_beta = m_gra_beta(ind);
  double gra_beta2 = gra_beta*gra_beta*gamma;
  double start = 0.1;
  double beta = betas(ind);
  arma::vec Xbeta_remain = Xbeta - X.col(ind)*beta;
  if(beta >= 0.001) start = sqrt(abs(beta/gra_beta))/2;
  double aa = start;
  double selected = beta;
  while(aa > 0){
    double aa2 = aa*aa;
    double beta_prime = beta + aa2*gra_beta;
    arma::vec X_beta_prime = Xbeta_remain + X.col(ind)*beta_prime;
    double lb_prime = LogLikelihood_multiple_beta_psi_mu2(Y, X_beta_prime, psi, mu2, n, sum_log_factorial_Y);
    if(lb_prime - ll - aa2*gra_beta2 > 0 ){ // | abs(lb_prime - ll - aa2*gra_beta2) < 0.0001) {
      selected = beta_prime;
      break;
    }
    aa = aa - start*down;
  }
 // if(selected >= 20) selected = beta;
  return selected;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_psi_with_multiple_beta(arma::vec Y, arma::vec Xbeta, double gra_psi, double ll, double psi, double mu2, int n, double gamma, double sum_log_factorial_Y, double down, double psi_cutoff){

  double gra_psi2 = gra_psi*gra_psi*gamma;
  double start = sqrt(abs(psi/gra_psi))/2;

  double aa = start;
  double selected = psi;
  while(aa > 0){
    double aa2 = aa*aa;
    double psi_prime = psi + aa2*gra_psi;
    if(psi_prime > 0){
      double lpsi_prime = LogLikelihood_multiple_beta_psi_mu2(Y, Xbeta, psi_prime, mu2, n, sum_log_factorial_Y);
      if(lpsi_prime - ll - aa2*gra_psi2 > 0) { // | abs(lpsi_prime - ll - aa2*gra_psi2) < 0.0001) {
        selected = psi_prime;
        break;
      }
    }
    aa = aa - start*down;
  }

  if(selected >= psi_cutoff) selected = psi_cutoff;
  return selected;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

double select_stepsize_for_mu2_with_multiple_beta(arma::vec Y, arma::vec Xbeta, double gra_mu2, double ll, double psi, double mu2, int n, double gamma, double sum_log_factorial_Y, double down){

  double gra_mu2_2 = gra_mu2*gra_mu2*gamma;
  double start = sqrt(abs(mu2/gra_mu2))/2;

  double aa = start;
  double selected = mu2;
  while(aa > 0){
    double aa2 = aa*aa;
    double mu2_prime = mu2 + aa2*gra_mu2;
    double lmu2_prime = LogLikelihood_multiple_beta_psi_mu2(Y, Xbeta, psi, mu2_prime, n, sum_log_factorial_Y);
    if(lmu2_prime - ll - aa2*gra_mu2_2 > 0 ){ // | abs(lmu2_prime - ll - aa2*gra_mu2_2) < 0.0001) {
      selected = mu2_prime;
      break;
    }
    aa = aa - start*down;
  }

  return selected;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec gradient_descent_multiple_beta_psi_mu2(arma::vec Y, arma::mat X, arma::vec beta, double psi, double mu2, double gamma,
                                                 int steps, double sum_log_factorial_Y, int n, int k, double down, double psi_cutoff){

  arma::vec Xbeta = X*beta;
  arma::vec gradient = gradient_multiple_beta_psi_mu2(Y, X, Xbeta, psi, mu2, n, k);
  arma::vec m_gra_beta = gradient(arma::span(2, 1 + k));
  double ll = LogLikelihood_multiple_beta_psi_mu2(Y, Xbeta, psi, mu2, n, sum_log_factorial_Y);

  arma::vec beta_prime = arma::zeros<arma::vec>(k); double psi_prime = 0; double mu2_prime = 0;
  for(int i = 0; i < steps; ++i){
    if(abs(gradient(0)) >= 0.0001){
      psi_prime = select_stepsize_for_psi_with_multiple_beta(Y, Xbeta, gradient(0), ll, psi, mu2, n, gamma, sum_log_factorial_Y, down, psi_cutoff);
    } else {
      psi_prime = psi;
    }
    if(abs(gradient(1)) >= 0.0001){
      mu2_prime = select_stepsize_for_mu2_with_multiple_beta(Y, Xbeta, gradient(1), ll, psi, mu2, n, gamma, sum_log_factorial_Y, down);
    } else {
      mu2_prime = mu2;
    }
    for(int j = 0; j < k; ++j){
      if(abs(gradient(j + 2)) >= 0.0001){
        beta_prime(j) = select_stepsize_for_multiple_beta(Y, X, beta, Xbeta, j, m_gra_beta, ll, psi, mu2, n, gamma, sum_log_factorial_Y, down);
      } else {
        beta_prime(j) = beta(j);
      }
    }

    beta = beta_prime; psi = psi_prime; mu2 = mu2_prime;
    Xbeta = X*beta;
    gradient = gradient_multiple_beta_psi_mu2(Y, X, Xbeta, psi, mu2, n, k);
    ll = LogLikelihood_multiple_beta_psi_mu2(Y, Xbeta, psi, mu2, n, sum_log_factorial_Y);

  }

  arma::vec est = arma::zeros<arma::vec>(3+k);
  est(0) = psi;
  est(1) = mu2;
  est(2) = ll;
  est(arma::span(3, 2 + k)) = beta;

  return est;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat Fisher_information(arma::vec Y, arma::mat X, arma::vec Xbeta, double psi, double mu2, int n, int k){

  arma::mat Fisher = arma::zeros<arma::mat>(2 + k,  2 + k);

  arma::vec e_term = exp(Xbeta + mu2);
  arma::vec e_term_fac_prime = 1/(e_term + psi);
  arma::vec e_term_fac = e_term%e_term_fac_prime;
  arma::vec e_term_fac_2 = 1/pow(e_term + psi, 2);
  arma::vec e_term_square_inverse = e_term%e_term_fac_2;
  arma::vec Y_i_psi_e_term_psi_inv2 = (Y + psi)%e_term_square_inverse;
  arma::mat X_square = pow(X, 2);

  arma::vec psi2_vec = arma::zeros<arma::vec>(n);
  for(int j = 0; j < n; ++j){
    psi2_vec(j) = - R::trigamma(Y(j) + psi) - 1/psi - (Y(j) + psi)*e_term_fac_2(j) + R::trigamma(psi) + 2*e_term_fac_prime(j);
  }

  Fisher(0, 0) = sum(psi2_vec); //psi, psi
  Fisher(0, 1) = Fisher(1, 0) = sum(e_term_fac) - sum(Y_i_psi_e_term_psi_inv2); //psi, mu2
  for(int j = 0; j < k; ++j){
    Fisher(0, j + 2) = Fisher(j + 2, 0) = sum(e_term_fac%X.col(j)) - sum(X.col(j)%Y_i_psi_e_term_psi_inv2); //psi, beta_k
  }
  Fisher(1, 1) = sum(psi*Y_i_psi_e_term_psi_inv2); //mu2, mu2
  for(int j = 0; j < k; ++j){
    Fisher(1, j + 2) = Fisher(j + 2, 1) = sum(psi*X.col(j)%Y_i_psi_e_term_psi_inv2); //mu2, beta_k
  }

  for(int j = 0; j < k; ++j){
    Fisher(j + 2, j + 2) = sum(psi*Y_i_psi_e_term_psi_inv2%X_square.col(j)); //beta_k, beta_k
  }

  for(int j = 0; j < k; ++j){
    for(int j_prime = j + 1; j_prime < k; ++j_prime){
      Fisher(j + 2, j_prime + 2) = Fisher(j_prime + 2, j + 2) = sum(psi*Y_i_psi_e_term_psi_inv2%X.col(j)%X.col(j_prime)); //beta_k, beta_k_prime
    }
  }

  return Fisher;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List PoissionGamma_multiple_beta(arma::vec Y, arma::mat X, arma::vec beta, double psi, double mu2, double gamma, int steps, double down, double psi_cutoff){

  arma::vec calculated_values = log_factorial_calculated(Y.max());

  int n = Y.n_elem;
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  double sum_log_factorial_Y = sum(log_factorial_Y);

  int k = beta.n_elem;
//  Rcpp::Rcout << "1" << std::endl;
  arma::vec est = gradient_descent_multiple_beta_psi_mu2(Y, X, beta, psi, mu2, gamma, steps, sum_log_factorial_Y, n, k, down, psi_cutoff);

//  Rcpp::Rcout << "2" << std::endl;
  psi = est(0);
  mu2 = est(1);
  beta = est(arma::span(3, 2 + k));
  arma::vec Xbeta = X*beta;
  arma::mat Fisher = Fisher_information(Y, X, Xbeta, psi, mu2, n, k);
 // Rcpp::Rcout << Fisher << std::endl;
  arma::mat est_cov = arma::inv(Fisher);

  arma::vec all_SE = arma::zeros<arma::vec>(2 + k);
  arma::vec all_pvalue = arma::zeros<arma::vec>(2 + k);

  for(int j = 0; j < k + 2; ++j){

    double SE = sqrt(abs(est_cov(j, j)));
    all_SE(j) = SE;
    if(j == 0){
      all_pvalue(j) = 2*R::pnorm(-abs(psi)/SE, 0, 1, TRUE, FALSE);
    }
    if(j == 1){
      all_pvalue(j) = 2*R::pnorm(-abs(mu2)/SE, 0, 1, TRUE, FALSE);
    }
    if(j > 1){
      all_pvalue(j) = 2*R::pnorm(-abs(beta(j - 2))/SE, 0, 1, TRUE, FALSE);
    }
  }


   return Rcpp::List::create(Rcpp::Named("psi") = psi,
                             Rcpp::Named("mu2") = mu2,
                             Rcpp::Named("beta") = beta,
                             Rcpp::Named("SE") = all_SE,
                             Rcpp::Named("p_value") = all_pvalue,
                             Rcpp::Named("likelihood") = est(2));
}


