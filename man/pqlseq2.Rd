% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PQLseq.R
\name{pqlseq2}
\alias{pqlseq2}
\title{Poisson/Binomial linear mixed model}
\usage{
pqlseq2(
  Y,
  x,
  K,
  W = NULL,
  lib_size = NULL,
  model = c("PMM", "BMM"),
  maxIter = 500,
  tol = 1e-05,
  ncores = 1,
  filter = TRUE,
  check_K = FALSE,
  nngp = FALSE,
  k = 10,
  fix_h2eq1 = FALSE,
  outfile = NULL,
  verbose = FALSE
)
}
\arguments{
\item{Y}{q by n count matrix for q genes/sites and n individuals}

\item{x}{n by 1 vector of the predicting variable}

\item{K}{n by n relatedness matrix}

\item{W}{n by c matrix of covariates for n individuals and c covariates}

\item{lib_size}{q by n matrix of total read count}

\item{model}{either "PMM" for poisson model or "BMM" for binomial model}

\item{maxIter}{maximum number of iterations for fitting the model}

\item{tol}{a threshold parameter to declare convergence}

\item{ncores}{number of cpus to use for parallel computing}

\item{filter}{whether to filter out genes that do not have at least two
individuals having read counts greater than 5}

\item{check_K}{whether to check if the relatedness matrix K is positive
definite or not and map it to the nearest PD matrix if not.}

\item{nngp}{whether to use NNGP to invert the H matrix (under testing)}

\item{k}{maximum number of nearest neighbors used in NNGP}

\item{fix_h2eq1}{whether to fix h2 to be 1}

\item{outfile}{also output results to a file on the disk if specified}

\item{verbose}{whether to print intermediate information for testing}
}
\value{
\code{pqlseq2} returns a data frame where each row stores results
for a gene/site. The results include:
\item{outcome}{name of the analyzed gene/site}
\item{n}{number of individuals}
\item{intercept}{intercept estimate}
\item{se_intercept}{standard error of the intercept estimate}
\item{beta}{fixed effect estimate for the predicting variable of interest}
\item{se_beta}{standard error of the beta estimate}
\item{pvalue}{Wald test p value for testing H0:beta = 0}
\item{h2}{heritability estimate}
\item{sigma2}{total variance component}
\item{converged}{whether the algorithm converged}
\item{elapsed_time}{time (seconds) for analyzing the gene/site}
}
\description{
Poisson/Binomial linear mixed model
}
