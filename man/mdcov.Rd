\name{mdcov}
\alias{mdcov}
\title{
Computation of the multidimensional distance covariance statistic for mutual independence using
characteristic functions.
}
\description{
Computation of the multidimensional distance covariance statistic for mutual independence using
characteristic functions. Compute the eigenvalues associated with the empirical covariance of the
limiting Gaussian procces. Compute the \eqn{p}-value associated with the test statistic, using the Imhof procedure.
}
\usage{
mdcov(X, vecd, a = 1, weight.choice = 1, N = 200, cubature = FALSE, K =
100, epsrel = 10 ^ -6, norming = TRUE, thresh.eigen = 10 ^ -8, estim.a =
FALSE, Cpp = TRUE, pval.comp = TRUE)
}
\arguments{
  \item{X}{
Data.frame or matrix with observations corresponding to rows and variables to columns.
}
  \item{vecd}{
a vector giving the sizes of each subvector.
}
  \item{a}{
parameter for the weight function.
}
  \item{weight.choice}{
Integer value in 1, 2, 3, 4, 5 corresponding to the choice in our paper.
}
  \item{N}{
Number of Monte-Carlo samples.
}
  \item{cubature}{
Logical. If \code{FALSE}, a Monte-Carlo approach is used. If \code{TRUE}, a cubature approach
is used.
}
  \item{K}{
Number of eigenvalues to compute.
}
\item{epsrel}{
  relative accuracy requested for the Imhof procedure.
}
\item{norming}{
  Logical. Should we normalize the test statistic with \eqn{H_n}.
}
\item{thresh.eigen}{
We will not compute eigenvalues (involved in the limiting distribution) below that threshold.
}
\item{estim.a}{
Logical. Should we automatically estimate the value of \eqn{a}.
}
\item{Cpp}{
Logical. If \code{TRUE} computations will be done using a fast C
code. The use of \code{FALSE} is only useful to compare the results with
the one given by the C code.
}
\item{pval.comp}{Logical. If \code{FALSE} do not compute the p-values
  and lambdas.}
}
\value{
  A list with the following components:
  \item{mdcov}{the value of the statistic \eqn{nT_n(w)}{n * Tn(w)} (this value
has been normed if \code{norming = TRUE})}
\item{Hn}{the denominator of \eqn{nT_n(w)}{nTn(w)}, namely \eqn{H_n}}
\item{pvalue}{the \eqn{p}-value of the test}
\item{lambdas}{the vector of eigenvalues computed (they have not been divided by
their sum)}
}
\author{
Lafaye de Micheaux P.
}

\references{Fan Y., Lafaye de Micheaux P., Penev S. and Salopek D. (2017). Multivariate nonparametric test of independence, Journal of Multivariate Analysis, 153, 189-210.
  }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\donttest{
a <- 1

# 4.1 Dependence among four discrete variables
set.seed(1)
n <- 100
w1 <- rpois(n, 1)
w3 <- rpois(n, 1)
w4 <- rpois(n, 1)
w6 <- rpois(n, 1)
w2 <- rpois(n, 3)
w5 <- rpois(n, 3)
x1 <- w1 + w2
x2 <- w2 + w3
x3 <- w4 + w5
x4 <- w5 + w6

X <- cbind(x1, x2, x3, x4)
mdcov(X, vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)
mdcov(X, vecd = c(1, 1, 1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)

X <- cbind(x1, x2)
mdcov(X, vecd = c(1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)

X <- cbind(x3, x4)
mdcov(X, vecd = c(1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)

# 4.2 Dependence between three bivariate vectors
set.seed(2)
n <- 200
Sigma <- matrix(c(
1 , 0 ,  0 ,  0 ,  0 ,  0 ,
0 , 1 ,  0 ,  0 ,  0 ,  0 ,
0 , 0 ,  1 ,  0 , .4 , .5 ,
0 , 0 ,  0 ,  1 , .1 , .2 ,
0 , 0 , .4 , .1 ,  1 ,  0 ,
0 , 0 , .5 , .2 ,  0 ,  1 ) ,
nrow = 6 , ncol = 6)   
W <- mvrnorm(n = n, mu = rep(0,6), Sigma = Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
mdcov(W, vecd = c(2, 2, 2), a, weight.choice = 1, N = 100, cubature = TRUE, epsrel = 10 ^ -7)

 # X^{(1)} with X^{(2)}^2
mdcov(W[,1:4], vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)

# X^{(2)} with X^{(3)}^2
mdcov(W[,2:6], vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)

# X^{(1)} with X^{(3)}^2
mdcov(W[,c(1:2, 4:6)], vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)

# 4.3 Four-dependent variables which are 2-independent and 3-independent
set.seed(3)
n <- 300
W <- sample(1:8, n, replace = TRUE)
X1 <- W \%in\% c(1, 2, 3, 5)
X2 <- W \%in\% c(1, 2, 4, 6)
X3 <- W \%in\% c(1, 3, 4, 7)
X4 <- W \%in\% c(2, 3, 4, 8)
X <- cbind(X1, X2, X3, X4)
   # pairwise independence
mdcov(X[,c(1, 2)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
mdcov(X[,c(1, 3)], vecd = c(1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)
mdcov(X[,c(1, 4)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
mdcov(X[,c(2, 3)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
mdcov(X[,c(2, 4)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
mdcov(X[,c(3, 4)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
   # 3-independence
mdcov(X[,c(1, 2, 3)], vecd = c(1, 1, 1), a, weight.choice =
                  1, cubature = TRUE)
mdcov(X[,c(1, 2, 4)], vecd = c(1, 1, 1), a, weight.choice =
                  1, cubature = TRUE)
mdcov(X[,c(1, 3, 4)], vecd = c(1, 1, 1), a, weight.choice =
                  1, cubature = TRUE)
mdcov(X[,c(2, 3, 4)], vecd = c(1, 1, 1), a, weight.choice =
                  1, cubature = TRUE)
   # 4-dependence
mdcov(X, vecd = c(1, 1, 1, 1), a, weight.choice = 1, cubature = TRUE)
}
}
