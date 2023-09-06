// Author: Zheng Li
// Date: 2023-07-28
// Poisson/Binomial mixed model using nearest neighbor Gaussian process

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

class nnPQL
{
private:
  struct Data
  {
    string model; // PMM or BMM
    int n; // number of samples
    int c; // number of covariates
    vec y; // n by 1 vector of count values
    vec lib_size; // n by 1 vector of library size
                  // total Bernoulli trials for Binomial models
    mat X;  // n by (c + 2) matrix (1, W, x) including
            // the intercept, covariates, and predictor
    mat K; // sample relatedness matrix
    sp_mat nn_mtx; // nearest neighbor matrix
    vec offset;
  } dat;
  
  struct Paras
  {
    vec alpha_beta; // (alpha^T, beta)^T
    vec taus; // tau1 and tau2
    vec alpha_beta_prev; // parameter estimates at the previous iteration
    vec taus_prev;
  } paras;
  
  struct Algorithm
  {
    vec d; // diagonal elements of D
    mat H;
    mat Hinv; // H^-1
    mat HinvX; // H^-1 * X
    mat XtHinvX; // X^T * H^-1 * X
    mat XtHinvX_inv; // (X^T * H^-1 * X)^-1
    mat P;
    vec y_tilde; // pseudo-data
    vec eta;
    vec mu;
  } alg;
  
  struct Control
  {
    bool nngp; // whether to use the nearest neighbor Gaussian process
    int maxIter; // maximum iterations of the iterative algorithm
    double tol; // tolerance used to declare convergence
    bool is_converged;
  } control;
  
  // save time spent at each updating step
  struct Time
  {
    wall_clock timer;
    double prev_time = 0;
    double curr_time = 0;
    vec elapsed_time;
  } time;
  
  int iter = 0;
  
public:
  void load_data(const string model, const vec &y, const mat &X, const mat &K, 
    const vec &lib_size, const sp_mat &nn_mtx)
  {
    dat.model = model;
    dat.y = y;
    dat.X = X;
    dat.n = dat.X.n_rows;
    dat.c = dat.X.n_cols - 2; // exclude intercept and predictor
    dat.K = K;
    dat.lib_size = lib_size;
    dat.nn_mtx = nn_mtx;
    if(dat.model == "PMM")
    {
      dat.offset = log(dat.lib_size);
    }
    else if(dat.model == "BMM")
    {
      dat.offset = zeros(dat.n);
    }
  }
  
  void set_control(const bool nngp, const int maxIter, const double tol)
  {
    control.nngp = nngp;
    control.maxIter = maxIter;
    control.tol = tol;
    control.is_converged = false;
  }
  
  void init_paras(const vec &init_alpha_beta)
  {
    paras.alpha_beta = init_alpha_beta;
    // initialize eta without random effects
    alg.eta = dat.X * paras.alpha_beta + dat.offset;
    update_d();
    update_y_tilde();
    paras.taus = ones(2) * min(0.9, var(alg.y_tilde) / 2.0);
    time.elapsed_time = vec(8, fill::zeros);
  }
  
  void update_eta()
  {
    alg.eta = alg.y_tilde - alg.Hinv * (alg.y_tilde - 
      dat.X * paras.alpha_beta) / alg.d + dat.offset;  
  }
  
  void update_d()
  {
    if(dat.model == "PMM")
    {
      alg.mu = exp(alg.eta);
      alg.d = alg.mu;
    }
    else if(dat.model == "BMM")
    {
      alg.mu = 1.0 / (1 + exp(-alg.eta)) % dat.lib_size;
      alg.d = alg.mu % (dat.lib_size - alg.mu) / dat.lib_size;
    }
  }
  
  void update_y_tilde()
  {
    alg.y_tilde = alg.eta - dat.offset + (dat.y - alg.mu) / alg.d;
  }
  
  void update_H(bool nngp)
  {
    if(nngp)
    {
      uvec idx_nbr; // neighbors indeces of each sample
      uvec idx_i(1);
      vec k_Nii; // K_{N(i), i}
      mat K_NiNi; // K_{N(i), i}
      double fii;
      vec bi;
      
      alg.Hinv = mat(dat.n, dat.n, fill::zeros);
      alg.Hinv(0,0) = 1.0 / (paras.taus(0) * dat.K(0,0) + paras.taus(1) + 
        1.0 / alg.d(0));
      for(int i = 1; i < dat.n; i++)
      {
        idx_i(0) = i;
        idx_nbr = find(dat.nn_mtx.row(i));
        k_Nii = dat.K(idx_nbr, idx_i);
        K_NiNi = dat.K(idx_nbr, idx_nbr);
        
        mat temp = paras.taus(0) * K_NiNi;
        temp.diag() += paras.taus(1) + 1.0 / alg.d(idx_nbr);
        temp = inv(temp);
        
        fii = paras.taus(0) * dat.K(i,i) + paras.taus(1) + 1.0 / alg.d(i);
        fii -= as_scalar(k_Nii.t() * temp  * k_Nii) * paras.taus(0) * 
          paras.taus(0);
        bi = paras.taus(0) * temp * k_Nii;
        
        alg.Hinv(idx_nbr, idx_nbr) += bi * bi.t() / fii;
        alg.Hinv(i,i) += 1.0 / fii;
        alg.Hinv(idx_i, idx_nbr) -= bi.t() / fii;
        alg.Hinv(idx_nbr, idx_i) -= bi / fii;
      }
    }
    else
    {
      alg.H = paras.taus(0) * dat.K;
      alg.H.diag() += 1.0 / alg.d + paras.taus(1);
      mat I = eye(dat.n, dat.n);
      alg.Hinv = solve(alg.H, I);
    }
  }
  
  void update_P()
  {
    alg.HinvX = alg.Hinv * dat.X;
    alg.XtHinvX = dat.X.t() * alg.HinvX;
    alg.XtHinvX_inv = inv(alg.XtHinvX);
    alg.P = alg.Hinv - alg.HinvX * alg.XtHinvX_inv * alg.HinvX.t();
  }
  
  void update_alpha_beta()
  {
    paras.alpha_beta_prev = paras.alpha_beta;
    paras.alpha_beta = alg.XtHinvX_inv * alg.HinvX.t() * alg.y_tilde;
  }
  
  void update_taus()
  {
    paras.taus_prev = paras.taus;
    vec scores(paras.taus.n_elem);
    mat AI(paras.taus.n_elem, paras.taus.n_elem);
    
    vec Py_tilde = alg.P * alg.y_tilde;
    vec PKPy_tilde = alg.P * dat.K * Py_tilde;
    mat PPy_tilde = alg.P * Py_tilde;
    AI(0,0) = accu(Py_tilde % (dat.K * PKPy_tilde));
    AI(1,1) = accu(Py_tilde % PPy_tilde);
    AI(0,1) = accu(PKPy_tilde % Py_tilde);
    AI(1,0) = AI(0,1);
    scores(0) = accu(PKPy_tilde % alg.y_tilde) - accu(alg.P % dat.K);
    scores(1) = accu(PPy_tilde % alg.y_tilde) - accu(alg.P.diag());
    scores /= 2.0;
    
    AI.diag() += 0.01; // to encourage positive definiteness
    paras.taus += solve(AI, scores);
    paras.taus.elem(find(paras.taus < 0)).zeros();
  }
  
  void record_time(int i_step)
  {
    time.curr_time = time.timer.toc();
    time.elapsed_time(i_step) += time.curr_time - time.prev_time;
    time.prev_time = time.curr_time;
  }
  
  void check_converge()
  {
    double diff1 = max(abs(paras.alpha_beta - paras.alpha_beta_prev) / 
      (abs(paras.alpha_beta) + abs(paras.alpha_beta_prev) + control.tol));
    double diff2 = max(abs(paras.taus - paras.taus_prev) / 
      (abs(paras.taus) + abs(paras.taus_prev) + control.tol));
    control.is_converged = (2 * max(diff1, diff2)) < control.tol;
  }
  
  void run(bool verbose)
  {
    time.timer.tic();
    while(iter < control.maxIter)
    {
      if(verbose)
      {
        cout << "Iter: " << iter + 1 << endl;
      }
      update_H(control.nngp);
      record_time(0);
      update_P();
      record_time(1);
      update_alpha_beta();
      record_time(2);
      update_taus();
      record_time(3);
      update_eta();
      record_time(4);
      update_d();
      record_time(5);
      update_y_tilde();
      record_time(6);
      check_converge();
      record_time(7);
      if(control.is_converged)
      {
        break;
      }
      iter += 1;
    }
  }
  
  List get_output()
  {
    List output;
    output = List::create(
      _["n"] = dat.n,
      _["tau1"] = paras.taus(0),
      _["tau2"] = paras.taus(1),
      _["h2"] = paras.taus(0) / sum(paras.taus),
      _["sigma2"] = sum(paras.taus),
      _["alpha"] = paras.alpha_beta.subvec(0, dat.c),
      _["beta"] = paras.alpha_beta.tail(1),
      _["se_beta"] = sqrt(alg.XtHinvX_inv(dat.c+1, dat.c+1)),
      _["converged"] = control.is_converged,
      _["niter"] = iter,
      _["elapsed_time"] = time.elapsed_time
    );
    return output;
  }
};


// [[Rcpp::export]]
List run_nnpql(const arma::vec &y, const arma::mat &X, const arma::mat &K,
  const arma::vec &init_alpha_beta, const std::string model, const int maxIter, 
  const double tol, const arma::vec &lib_size, const bool nngp, 
  const arma::sp_mat &nn_mtx, const bool verbose)
{
  wall_clock timer;
  timer.tic();
  nnPQL nnpql_model;
  
  nnpql_model.load_data(model, y, X, K, lib_size, nn_mtx);
  nnpql_model.set_control(nngp, maxIter, tol);
  nnpql_model.init_paras(init_alpha_beta);
  nnpql_model.run(verbose);
  
  List output = nnpql_model.get_output();
  double elapsed = timer.toc();
  output.push_back(elapsed, "elapsed_time");
  return output;
}

