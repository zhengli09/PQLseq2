// Author: Zheng Li
// Date: 2023-07-28
// Poisson/Binomial mixed model using nearest neighbor Gaussian process

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <float.h>
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
    field<mat> Ks; // sample relatedness/kernel matrices
    sp_mat nn_mtx; // nearest neighbor matrix
    mat In; // n by n identity matrix
    vec offset;
  } dat;
  
  struct Paras
  {
    vec alpha_beta; // (alpha^T, beta)^T
    vec alpha_beta_prev; // parameter estimates at the previous iteration
    vec taus; // tau1 and tau2
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
    uvec est_idx; // whether to estimate or fix the parameters
    bool nngp; // whether to use the nearest neighbor Gaussian process
    int maxIter; // maximum iterations of the iterative algorithm
    double tol; // tolerance used to declare convergence
    bool is_converged;
    bool is_failed;
  } control;
  
  struct IntEst
  {
    mat taus;
  } int_est;
  
  // save time spent at each updating step
  struct Time
  {
    wall_clock timer;
    double prev_time = 0;
    double curr_time = 0;
    vec step_time;
  } time;
  
  int iter;
  double step_size;
  
public:
  void load_data(const string model, const vec &y, const mat &X, 
    const field<mat> &Ks, const vec &lib_size, const sp_mat &nn_mtx)
  {
    dat.model = model;
    dat.y = y;
    dat.X = X;
    dat.n = dat.X.n_rows;
    dat.c = dat.X.n_cols - 2; // exclude intercept and predictor
    dat.Ks = Ks;
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
    dat.In = eye(dat.n, dat.n);
  }
  
  void set_control(const bool nngp, const int maxIter, const double tol)
  {
    control.nngp = nngp;
    control.maxIter = maxIter;
    control.tol = tol;
    control.is_converged = false;
    control.is_failed = false;
  }
  
  void set_control(const uvec est_idx)
  {
    control.est_idx = est_idx;
    control.is_converged = false;
    control.is_failed = false;
  }
  
  void update_step_size()
  {
    if((iter+1) % 10 == 0)
    {
      step_size = 0.9 * step_size;
    }
  }
  
  void init_paras(const vec &init_alpha_beta, const vec &init_taus)
  {
    paras.alpha_beta = init_alpha_beta;
    paras.taus = init_taus;
    // initialize eta without random effects
    alg.eta = dat.X * paras.alpha_beta + dat.offset;
    update_d();
    update_y_tilde();
    // further initialize taus
    paras.taus(control.est_idx) = min(0.9, var(alg.y_tilde) / 
      control.est_idx.n_elem) * ones(control.est_idx.n_elem);
    update_H(control.nngp);
    update_P();
    vec Py_tilde = alg.P * alg.y_tilde;
    for(int i = 0; i < control.est_idx.n_elem; i++)
    {
      int l = control.est_idx(i);
      paras.taus(l) += (accu(Py_tilde % (dat.Ks(l) * Py_tilde)) - 
        accu(alg.P % dat.Ks(l))) * paras.taus(l) * paras.taus(l);
      paras.taus(l) = max(0.0, paras.taus(l) / dat.n);
    }
    time.step_time = vec(7, fill::zeros);
    iter = 0;
    step_size = 1;
  }
  
  void init_int_est()
  {
    int_est.taus.zeros(paras.taus.n_elem, control.maxIter);
  }
  
  void update_eta()
  {
    alg.eta = alg.y_tilde - (alg.Hinv * (alg.y_tilde - 
      dat.X * paras.alpha_beta)) / alg.d + dat.offset;
    // follow glm family$mu.eta function to handle edge cases
    for(int i = 0; i < dat.n; i++)
    {
      if(alg.eta(i) < -30 || alg.eta(i) > 30)
      {
        alg.eta(i) = DBL_EPSILON;
      }
    }
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
      uvec idx_nbr; // neighbors indexes of each sample
      uvec idx_i(1);
      double k_ii;
      vec K_Nii; // K_{N(i), i}
      mat K_NiNi; // K_{N(i), i}
      mat K_NiNi_inv;
      double fii;
      vec bi;

      alg.Hinv = mat(dat.n, dat.n, fill::zeros);
      for(int i = 0; i < control.est_idx.n_elem; i++)
      {
        int l = control.est_idx(i);
        alg.Hinv(0,0) += paras.taus(l) * dat.Ks(l)(0,0);
      }
      alg.Hinv(0,0) = 1.0 / (1.0 / alg.d(0) + alg.Hinv(0,0));
      for(int i = 1; i < dat.n; i++)
      {
        idx_i(0) = i;
        idx_nbr = find(dat.nn_mtx.row(i));
        k_ii = 1.0 / alg.d(i);
        K_Nii = vec(idx_nbr.n_elem, fill::zeros);
        K_NiNi = mat(idx_nbr.n_elem, idx_nbr.n_elem, fill::zeros);
        K_NiNi.diag() = 1.0 / alg.d(idx_nbr);
        for(int j = 0; j < control.est_idx.n_elem; j++)
        {
          int l = control.est_idx(j);
          k_ii += paras.taus(l) * dat.Ks(l)(i,i);
          K_Nii += paras.taus(l) * dat.Ks(l)(idx_nbr, idx_i);
          K_NiNi += paras.taus(l) * dat.Ks(l)(idx_nbr, idx_nbr);
        }
        K_NiNi_inv = inv(K_NiNi);

        fii = as_scalar(k_ii - K_Nii.t() * K_NiNi_inv * K_Nii);
        bi = K_NiNi_inv * K_Nii;

        alg.Hinv(idx_nbr, idx_nbr) += bi * bi.t() / fii;
        alg.Hinv(i,i) += 1.0 / fii;
        alg.Hinv(idx_i, idx_nbr) -= bi.t() / fii;
        alg.Hinv(idx_nbr, idx_i) -= bi / fii;
      }
    }
    else
    {
      alg.H = diagmat(1.0 / alg.d);
      for(int i = 0; i < control.est_idx.n_elem; i++)
      {
        int l = control.est_idx(i);
        alg.H += paras.taus(l) * dat.Ks(l);
      }
      alg.Hinv = solve(alg.H, dat.In);
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
    vec scores(control.est_idx.n_elem);
    mat AI(control.est_idx.n_elem, control.est_idx.n_elem);
    
    vec Py_tilde = alg.P * alg.y_tilde;
    for(int i = 0; i < control.est_idx.n_elem; i++)
    {
      int l = control.est_idx(i);
      vec KlPy_tilde = dat.Ks(l) * Py_tilde;
      scores(i) = accu(Py_tilde % KlPy_tilde) - accu(alg.P % dat.Ks(l));
      for(int j = 0; j <= i; j++)
      {
        int m = control.est_idx(j);
        vec KmPy_tilde = dat.Ks(m) * Py_tilde;
        AI(i,j) = accu(KlPy_tilde % (alg.P * KmPy_tilde));
        if(i != j)
        {
          AI(j,i) = AI(i,j);
        }
      }
    }
    paras.taus(control.est_idx) += step_size * solve(AI, scores);
    
    // handle tau less than zero case
    paras.taus.elem(find((paras.taus_prev < control.tol) %
      (paras.taus < control.tol))).zeros();
    double ss = step_size;
    while(any(paras.taus < 0.0))
    {
      ss = ss * 0.5;
      paras.taus(control.est_idx) = ss * solve(AI, scores) +
        paras.taus_prev(control.est_idx);
      paras.taus.elem(find((paras.taus_prev < control.tol) %
        (paras.taus < control.tol))).zeros();
    }
    paras.taus.elem(find(paras.taus < control.tol)).zeros();
    int_est.taus.col(iter) = paras.taus;
  }
  
  void record_time(int i_step)
  {
    time.curr_time = time.timer.toc();
    time.step_time(i_step) += time.curr_time - time.prev_time;
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
  
  void check_anomaly()
  {
    if(iter >= 9)
    {
      vec check_zeros =  sum(int_est.taus.cols(iter-9,iter), 1);
      int n_zeros = accu(check_zeros == 0);
      if(n_zeros > 0 && n_zeros < control.est_idx.n_elem)
      {
        control.is_failed = true;
      }
    }
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
      update_eta();
      record_time(3);
      update_taus();
      record_time(4);
      update_d();
      record_time(5);
      update_y_tilde();
      record_time(6);
      check_converge();
      check_anomaly();
      if(control.is_converged || control.is_failed)
      {
        break;
      }
      update_step_size();
      iter += 1;
    }
  }
  
  List get_output()
  {
    List output;
    vec u = alg.eta - dat.X * paras.alpha_beta;
    output = List::create(
      _["n"] = dat.n,
      _["taus"] = paras.taus,
      _["tau1"] = paras.taus(0),
      _["tau2"] = paras.taus(1),
      _["h2"] = paras.taus(0) / sum(paras.taus),
      _["sigma2"] = sum(paras.taus),
      _["alpha"] = paras.alpha_beta.subvec(0, dat.c),
      _["beta"] = paras.alpha_beta.tail(1),
      _["se_beta"] = sqrt(alg.XtHinvX_inv(dat.c+1, dat.c+1)),
      _["converged"] = control.is_converged,
      _["niter"] = iter,
      _["step_time"] = time.step_time,
      _["u"] = u,
      _["eta"] = alg.eta,
      _["mu"] = alg.mu,
      _["H"] = alg.H
    );
    return output;
  }
};


// [[Rcpp::export]]
List run_nnpql(const arma::vec &y, const arma::mat &X, 
  const arma::field<arma::mat> &Ks, const arma::vec &init_alpha_beta, 
  const std::string model, const int maxIter, const double tol, 
  const arma::vec &lib_size, const bool nngp, const arma::sp_mat &nn_mtx, 
  const bool fix_h2eq1, const bool verbose)
{
  wall_clock timer;
  timer.tic();
  List output;
  uvec est_idx_curr;
  uvec est_idx_prev;
  vec init_taus;
  
  est_idx_curr = regspace<uvec>(0, Ks.n_elem-1);
  if(fix_h2eq1)
  {
    est_idx_curr = est_idx_curr.subvec(0, Ks.n_elem-2);
  }
  init_taus = vec(Ks.n_elem, fill::zeros);
  
  nnPQL nnpql_model;
  nnpql_model.load_data(model, y, X, Ks, lib_size, nn_mtx);
  nnpql_model.set_control(nngp, maxIter, tol);
  do
  {
    nnpql_model.set_control(est_idx_curr);
    nnpql_model.init_paras(init_alpha_beta, init_taus);
    nnpql_model.init_int_est();
    nnpql_model.run(verbose);
    output = nnpql_model.get_output();
    init_taus = as<vec>(output["taus"]);
    est_idx_prev = est_idx_curr;
    est_idx_curr = find(init_taus > tol);
  } while(est_idx_curr.n_elem != est_idx_prev.n_elem || 
    any(est_idx_curr != est_idx_prev));

  double elapsed = timer.toc();
  output.push_back(elapsed, "elapsed_time");
  return output;
}

