\name{dependogram}
\alias{dependogram}
\title{Nonparametric tests of independence between random vectors}
\description{This function can be used for the following two problems:
  1) testing mutual independence between many
  numerical random vectors, and 2) testing for serial independence of a multivariate
  stationary numerical time series. The proposed test does not assume
  continuous marginals. It is valid for any probability distribution. It
is also invariant with respect to the affine general linear group of transformations on
the vectors. This test is based on a
characterization of mutual independence defined from probabilities
of half-spaces in a combinatorial formula of Mobius. As such,
it is a natural generalization of tests of independence between
univariate random variables using the empirical distribution
function. Without the assumption that each vector is
one-dimensional with a continuous cumulative distribution
function, any test of independence can not be distribution free.
The critical values of the proposed test are thus computed with
the bootstrap which was shown to be consistent in this context.}
\usage{
dependogram(X,vecd.ou.p,N=10,B=2000,alpha=0.05,display=TRUE,graphics=TRUE,nbclus=1)
}
\arguments{
  \item{X}{Data.frame or matrix
    with observations corresponding to rows and variables to columns.}
  \item{vecd.ou.p}{For the mutual independence problem 1), a vector
    giving the sizes of each subvector. For the serial independence
    problem 2), an integer indicating the number of consecutive observations.}
  \item{N}{Integer. Number of points of the discretization to obtain
    directions on the sphere in order to 
    evaluate the value of the test statistic.}
  \item{B}{Integer. Number of bootstrap replications. Note that B can be slightly modified if \code{nbclus}>1}
  \item{alpha}{Double. Level of the test.}
  \item{display}{Logical. TRUE to display the values of the A-dependence statistic.}
\item{graphics}{Logical. TRUE to plot the dependogram.}
\item{nbclus}{Integer. Number of nodes in the cluster. Used only for parallel computations.}
}
\value{
 A list with the following components:

  In the mutual independence case:
  
  \item{norm.RnA}{ ... should be completed ... }
  \item{Rn}{ ... should be completed ... }
  \item{rA}{ ... should be completed ... }
  \item{r}{ ... should be completed ... }
  \item{RnAsstar}{ ... should be completed ... }

  In the serial case:

  \item{norm.SnA}{ ... should be completed ... }
  \item{Sn}{ ... should be completed ... }
  \item{sA}{ ... should be completed ... }
  \item{s}{ ... should be completed ... }
  \item{RnAsstar}{ ... should be completed ... }
  
}
\author{M. Bilodeau, P. Lafaye de Micheaux}

\examples{
%% All these examples are taken from the paper
%% "Nonparametric tests of independence between random vectors"
%%  R. Beran, M. Bilodeau, P. Lafaye de Micheaux
%%  Journal of Multivariate Analysis 98 (2007) 1805-1824

% NOTE: In real applications, B should be increased to at least 1000

% Example 4.1: Dependence among four discrete variables
% NOTE: with B=1000, this one took 65s with nbclus=1 and 15s with nbclus=7 on my computer
n <- 100
W1 <- rpois(n,1)
W3 <- rpois(n,1)
W4 <- rpois(n,1)
W6 <- rpois(n,1)
W2 <- rpois(n,3)
W5 <- rpois(n,3)
X1 <- W1 + W2
X2 <- W2 + W3
X3 <- W4 + W5
X4 <- W5 + W6
X <- cbind(X1,X2,X3,X4)
dependogram(X,vecd.ou.p=c(1,1,1,1),N=10,B=20,alpha=0.05,display=TRUE,graphics=TRUE)

% Example 4.2: Dependence between three bivariate vectors
n <- 50
Sigma <- matrix(c(1, 0, 0, 0, 0, 0,
                  0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0,.4,.5,
                  0, 0, 0, 1,.1,.2,
                  0, 0,.4,.1, 1, 0,
                  0, 0,.5,.2, 0, 1),nrow=6,ncol=6,byrow=TRUE)
W <- chol(Sigma) %*% matrix(rnorm(6*n),nrow=6,ncol=n)
X1 <- cbind(W[1,],W[2,])
X2 <- cbind(W[3,],W[4,])
X3 <- cbind(W[5,],W[6,])
X <- cbind(X1,X2,X3)
\donttest{
dependogram(X,vecd.ou.p=c(2,2,2),N=10,B=20,alpha=0.05,display=TRUE,graphics=TRUE)
}

% Example 4.3: Four-dependent variables which are 2-independent and 3-independent
n <- 100
W <- sample(x=1:8,size=n,TRUE)
X1 <- W %in% c(1,2,3,5)
X2 <- W %in% c(1,2,4,6)
X3 <- W %in% c(1,3,4,7)
X4 <- W %in% c(2,3,4,8)
X <- cbind(X1,X2,X3,X4)
dependogram(X,vecd.ou.p=c(1,1,1,1),N=10,B=20,alpha=0.05,display=TRUE,graphics=TRUE)

% Example 4.4: Serial independence of a binary sequence of zeros and ones
n <- 100
W <- rbinom(n,1,0.8)
Y <- W[1:(n-3)]*W[4:n]
dependogram(W,vecd.ou.p=4,N=10,B=20,alpha=0.05,display=TRUE,graphics=TRUE)
dependogram(Y,vecd.ou.p=4,N=10,B=20,alpha=0.05,display=TRUE,graphics=TRUE)

% Example 4.5: Serial independence of directional data
n <- 75
U <- matrix(rnorm(2*n),nrow=n,ncol=2)
W <- U[1:(n-1),] + sqrt(2)*U[2:n,]
Y <- W/apply(W,MARGIN=1,FUN=function(x) sqrt(x[1]^2+x[2]^2))
\donttest{
dependogram(Y,vecd.ou.p=3,N=10,B=20,alpha=0.05,display=TRUE,graphics=TRUE)
}
}
