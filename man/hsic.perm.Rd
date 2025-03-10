% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kpc.r
\name{hsic.perm}
\alias{hsic.perm}
\title{Hilber Schmidt Independence Criterion permutation test}
\usage{
hsic.perm(x, y, sig = 1, p = 100, numCol = 50)
}
\arguments{
\item{x}{data of first sample}

\item{y}{data of second sample}

\item{sig}{Gaussian kernel width for HSIC tests. Default is 1}

\item{p}{Number of permutations. Default is 100}

\item{numCol}{maximum number of columns that we use for the incomplete Cholesky decomposition}
}
\value{
hsic.perm() returns a list with class htest containing
\item{method}{description of test}
\item{statistic}{observed value of the test statistic}
\item{estimate}{HSIC(x,y)}
\item{estimates}{a vector: \link{HSIC(x,y), mean of HSIC(x,y), variance of HSIC(x,y)}}
\item{replicates}{replicates of the test statistic}
\item{p.value}{approximate p-value of the test}
\item{data.name}{desciption of data}
}
\description{
Test to check the independence between two variables x and y using HSIC.
The hsic.perm() function, uses Hilbert-Schmidt independence criterion to test for independence between
random variables.
}
\details{
Let x and y be two samples of length n. Gram matrices K and L are defined as: \eqn{K_{i,j} = \exp\frac{(x_i-x_j)^2}{\sigma^2}}{K_{i,j} =exp((x_i-x_j)^2/sig^2)} and \eqn{L_{i,j} = \exp\frac{(y_i-y_j)^2}{\sigma^2}}{L_{i,j} =exp((y_i-y_j)^2/sig^2)}. \eqn{H_{i,j} = \delta_{i,j} - \frac{1}{n}}{H_{i,j} = \delta_{i,j} - 1/n}. Let \eqn{A=HKH} and \eqn{B=HLH}, then \eqn{HSIC(x,y)=\frac{1}{n^2}Tr(AB)}{HSIC(x,y)=Tr(AB)/n^2}. Permutation test permutes y p times to get \eqn{y_{(p)}} and calculates HSIC(x,y_{(p)}). \eqn{pval = \frac{1(HSIC(x,y)>HSIC(x,y_{(p)}))}{p}}{(HSIC(x,y)>HSIC(x,y_{(p)}))/p}.
}
\examples{
library(energy)
set.seed(10)
#independence
x <- runif(300)
y <- runif(300)

hsic.gamma(x,y)
hsic.perm(x,y)
dcov.gamma(x,y)
dcov.test(x,y)

#uncorelated but not dependent
z <- 10*(runif(300)-0.5)
w <- z^2 + 10*runif(300)

cor(z,w)
hsic.gamma(z,w)
hsic.perm(z,w)
dcov.gamma(z,w)
dcov.test(z,w)
}
\references{
A. Gretton et al. (2005). Kernel Methods for Measuring Independence. JMLR 6 (2005) 2075-2129.
}
\seealso{
\link{hsic.gamma}, \link{hsic.clust}, \link{kernelCItest}
}
\author{
Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
}
