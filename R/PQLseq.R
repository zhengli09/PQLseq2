# Author: Zheng Li
# Date: 2023-08-12
# A faster re-implementation of the Poisson/Binomial mixed model

#' Poisson/Binomial linear mixed model
#'
#' @param Y q by n count matrix for q genes/sites and n individuals 
#' @param x n by 1 vector of the predicting variable
#' @param K n by n relatedness matrix 
#' @param W n by c matrix of covariates for n individuals and c covariates
#' @param lib_size q by n matrix of total read count
#' @param model either "PMM" for poisson model or "BMM" for binomial model
#' @param maxIter maximum number of iterations for fitting the model
#' @param tol a threshold parameter to declare convergence
#' @param ncores number of cpus to use for parallel computing
#' @param filter whether to filter out genes that do not have at least two
#'   individuals having read counts greater than 5
#' @param check_K whether to check if the relatedness matrix K is positive
#'   definite or not and map it to the nearest PD matrix if not.
#' @param nngp whether to use NNGP to invert the H matrix (under testing)
#' @param k maximum number of nearest neighbors used in NNGP
#' @param fix_h2eq1 whether to fix h2 to be 1
#' @param outfile also output results to a file on the disk if specified
#' @param verbose whether to print intermediate information for testing
#' 
#' @return \code{pqlseq2} returns a data frame where each row stores results
#' for a gene/site. The results include:
#' \item{outcome}{name of the analyzed gene/site}
#' \item{n}{number of individuals}
#' \item{intercept}{intercept estimate}
#' \item{se_intercept}{standard error of the intercept estimate}
#' \item{beta}{fixed effect estimate for the predicting variable of interest}
#' \item{se_beta}{standard error of the beta estimate}
#' \item{pvalue}{Wald test p value for testing H0:beta = 0}
#' \item{h2}{heritability estimate}
#' \item{sigma2}{total variance component}
#' \item{converged}{whether the algorithm converged}
#' \item{elapsed_time}{time (seconds) for analyzing the gene/site}
#' 
#' @export
pqlseq2 <- function(Y, x, K, W = NULL, lib_size = NULL, model = c("PMM", "BMM"), 
  maxIter = 500, tol = 1e-5, ncores = 1, filter = TRUE, check_K = FALSE, 
  nngp = FALSE, k = 10, fix_h2eq1 = FALSE, outfile = NULL, verbose = FALSE)
{
  # 1.check inputs
  model <- match.arg(model)
  if(ncores > parallel::detectCores()){
    warning("the specified number of cores is larger than that detected")
    ncores <- parallel::detectCores() - 1
  }
  if(model == "BMM"){
    if(is.null(lib_size)){
      stop("please input the library size (lib_size)")
    }
    if(!all(dim(lib_size) == dim(Y))){
      stop(paste("the dimension of the count matrix and library size matrix",
        "does not match"))
    }
  }
  
  # 2.filter genes/sites
  if(filter & model == "PMM"){
    idx_keep <- apply(Y, MARGIN = 1, function(y){
      sum(y > 5) >= 2
    })
    Y <- Y[idx_keep, , drop = FALSE]
  }
  if(model == "BMM"){
    # exclude sites that have count values exceed the library size
    p_succ <- Y / lib_size
    p_succ[is.na(p_succ)] <- 0
    idx_rm <- apply(p_succ > 1, MARGIN = 1, sum) > 0
    Y <- Y[!idx_rm, , drop = FALSE]
    lib_size <- lib_size[!idx_rm, , drop = FALSE]
  }
  Y <- t(as.matrix(Y))
  lib_size <- t(as.matrix(lib_size))
  x <- as.matrix(x)
  n <- nrow(Y)
  q <- ncol(Y)
  
  # 3.process the covariate matrix W
  if(!is.null(W)){
    # remove the intercept
    if(all(W[, 1] == 1)){
      W <- W[, -1, drop = FALSE]
    }
    W <- as.matrix(W)
    c <- ncol(W)
  } else{
    c <- 0
  }
  
  cat("## number of individuals:", n, "\n")
  cat("## number of genes/sites:", q, "\n")
  cat("## number of covariates:", c, "\n")
  
  # 4.process the sample relatedness matrix
  K <- as.matrix(K)
  if(check_K)
  {
    eigvals <- eigen(K, only.values = TRUE)$values
    if(any(eigvals < 1e-10)){ 
      warning(paste("K is singular, approximate K with its",
        "nearest positive definite matrix"))
      K <- as.matrix(Matrix::nearPD(K, doSym = TRUE, corr = TRUE)$mat)	
    } 
  }
  
  # 5.construct a adjacency matrix for NNGP
  # note1: samples with more relatives are better to have smaller index
  # note2: need dgCMatrix
  if(nngp){
    nbr_deg <- apply(K, MARGIN = 1, function(x){sum(abs(x))})
    idx_order <- order(nbr_deg, decreasing = TRUE)
    Y <- Y[idx_order, , drop = FALSE]
    lib_size <- lib_size[idx_order, , drop = FALSE]
    x <- x[idx_order]
    W <- W[idx_order, , drop = FALSE]
    K <- K[idx_order, idx_order, drop = FALSE]
    
    nn_mtx <- const_nnMtx(K, k)
  } else{
    nn_mtx <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = T, doDiag = F)
    nn_mtx <- as(nn_mtx, "generalMatrix")
  }

  # 6.run algorithm
  header <- c("outcome", "n", "intercept", "se_intercept", "beta", "se_beta", 
    "pvalue", "h2", "sigma2", "converged", "elapsed_time")
  if(!is.null(outfile)){
    write.table(t(header), file = outfile, col.names = F, row.names = F,
      quote = F, sep = ',')
  }
  # 6.1.Binomial Mixed Model
  if(model == "BMM"){
    cat("# fitting binomial mixed model ... \n")
    res <- pbmcapply::pbmclapply(1:q, function(i){
      if(verbose){
        cat("NO. Gene/Site =", i, "\n")
      }
      
      p_succ <- Y[, i] / lib_size[, i]
      mod0_formula <- "p_succ ~ x"
      if(c != 0){
        mod0_formula <- "p_succ ~ W + x"
      }
      # fit model under g = e = 0 to obtain initials 
      mod0 <- glm(as.formula(mod0_formula), family = binomial(link = "logit"),
        weights = lib_size[, i], na.action = na.omit)
      init_alpha_beta <- mod0$coefficients
      
      # focus on the complete set of data
      idx_keep <- match(
        rownames(model.frame(as.formula(mod0_formula), na.action = na.omit)),
        rownames(model.frame(as.formula(mod0_formula), na.action = na.pass)))
      # after excluding NA cases, find the covariates that are homogeneous
      Wx <- model.matrix(mod0)
      covs_homo <- sapply(2:ncol(Wx), function(j){
        length(unique(Wx[, j])) == 1
      })
      if(any(covs_homo)){
        warning(paste0("the ", paste0(which(covs_homo), collapse = ","), "-th ",
          "column of covariates are the same for gene/site ", colnames(Y)[i]))
      }
      
      # fit model
      if(!any(covs_homo)){
        res <- tryCatch({
          Ks <- list(K[idx_keep, idx_keep], diag(length(idx_keep)))
          run_pql(Y[idx_keep, i], Wx, Ks, init_alpha_beta, model, maxIter, 
            tol, lib_size[idx_keep, i], nngp, nn_mtx, fix_h2eq1, verbose)
        }, error = function(e){
          print(e)
          NA
        })
        
        # organized results into an output object
        if(is.list(res)){
          beta <- res$beta
          if(verbose){
            se_beta <- sqrt(res$cov[nrow(res$cov), nrow(res$cov)])
            se_intercept <- sqrt(res$cov[1, 1])
          } else{
            se_beta <- res$se_beta
            se_intercept <- res$se_intercept
          }
          pvalue <- pchisq((beta / se_beta)^2, df = 1, lower.tail = F)
          
          # handle cases when both variance components are zeros
          if(fix_h2eq1 & all(res$taus == 0)){
            res$h2 <- 1
          }
          
          ests <- data.frame(outcome = colnames(Y)[i], n = res$n, 
            intercept = res$intercept, se_intercept = se_intercept,
            beta = beta, se_beta = se_beta, pvalue = pvalue, h2 = res$h2, 
            sigma2 = res$sigma2, converged = res$converged, 
            elapsed_time = res$elapsed_time)
          if(verbose){
            others <- list(taus = res$taus, P = res$P, cov = res$cov, 
              alpha = res$alpha, se_alpha = sqrt(diag(res$cov))[-nrow(res$cov)], 
              eta = res$eta, H = res$H, intercept = res$intercept, 
              se_intercept = sqrt(res$cov[1, 1]), beta = res$beta, 
              se_beta = sqrt(res$cov[nrow(res$cov), nrow(res$cov)]),
              pvalue = pvalue, tau1 = res$tau1, tau2 = res$tau2, 
              sigma2 = res$sigma2, h2 = res$h2, iter = res$niter, mu = res$mu, 
              Y = Y, y = Y[idx_keep, i], res = Y[idx_keep, i] - res$mu, 
              converged = res$converged, X = Wx, K = K, 
              V = res$tau1 * Ks[[1]] + res$tau2 * Ks[[2]], D = diag(res$d[, 1]), 
              family = mod0$family, numTotal = lib_size[idx_keep, i], u = res$u, 
              LL = res$LL, AIC = res$AIC)
          }
        } else{
          ests <- c(colnames(Y)[i], rep(NA, 10))
          ests <- data.frame(t(ests))
          colnames(ests) <- header
          others <- list()
        }
        
        if(verbose){
          output <- list(estimates = ests, others = others)
        } else{
          output <- ests
        }

        if(!is.null(outfile)){
          write.table(ests, file = outfile, col.names = F, row.names = F, 
            quote = F, sep = ',', append = T)
        }
      }
      return(output)
    }, mc.cores = ncores)
    
    if(!is.null(names(res))){
      if(names(res)[2] == "warning"){
        res <- res[["value"]]
      }
    }
    if(verbose){
      res_ests <- do.call(rbind, lapply(res, `[[`, "estimates"))
      res_othrs <- lapply(res, `[[`, "others")
      names(res_othrs) <- res_ests$outcome
      res <- new(Class = "output", estimates = res_ests, others = res_othrs)
    } else{
      res <- do.call(rbind, res)
    }
    res
  }
}


const_nnMtx <- function(K, k)
{
  n <- nrow(K)
  nn_mtx <- Matrix::Matrix(0, nrow = n, ncol = n, sparse = T, doDiag = F)
  for(i in 1:n)
  {
    idx_topk <- head(order(abs(K[i, 1:(i-1)]), decreasing = T), k)
    nn_mtx[i, idx_topk] <- 1
  }
  nn_mtx[1, 1] <- 0
  return(nn_mtx)  
}


setClass("output", slots = list(
  estimates = "data.frame",
  others = "list"  
))

