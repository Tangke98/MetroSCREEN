#'
#' Test to check the independence between two variables x and y using the Distance Covariance.
#' The dcov.gamma() function, uses Distance Covariance independence criterion with gamma approximation to test for independence between two random variables.
#'
#' @details Let x and y be two samples of length n. Gram matrices K and L are defined as: \eqn{K_{i,j} = \| x_i-x_j \|^s}{K_{i,j} =|x_i-x_j|^s} and \eqn{L_{i,j} = \| y_i-y_j \|^s}{L_{i,j} =|y_i-y_j|^s}, where 0<s<2. \eqn{H_{i,j} = \delta_{i,j} - \frac{1}{n}}{H_{i,j} = delta_{i,j} - 1/n}. Let A=HKH and B=HLH, then \eqn{nV^2=\frac{1}{n^2}\sum A_{i,j} B_{i,j}}{nV^2 = \sum A_{i,j} B_{i,j} \n^2}. For more detail: \link{dcov.test} in package energy. Gamma test compares \eqn{nV^2_n(x,y)} with the \eqn{\alpha}{alpha} quantile of the gamma distribution with mean and variance same as \eqn{nV^2_n} under independence hypothesis.
#' @param x data of first sample
#' @param y data of second sample
#' @param index exponent on Euclidean distance, in (0,2]
#' @param numCol Number of columns used in incomplete Singular Value Decomposition
#'
#' @references A. Gretton et al. (2005). Kernel Methods for Measuring Independence. JMLR 6 (2005) 2075-2129.
#' @references G. Szekely, M. Rizzo and N. Bakirov (2007). Measuring and Testing Dependence by Correlation of Distances. The Annals of Statistics 2007, Vol. 35, No. 6, 2769-2794.
#'
#' @importFrom RSpectra eigs
#' @importFrom stats var pgamma dist
#' @import methods
#'
#' @export
#' @return dcov.gamma() returns a list with class htest containing
#' \item{method}{description of test}
#' \item{statistic}{observed value of the test statistic}
#' \item{estimate}{nV^2(x,y)}
#' \item{estimates}{a vector: [nV^2(x,y), mean of nV^2(x,y), variance of nV^2(x,y)]}
#' \item{replicates}{replicates of the test statistic}
#' \item{p.value}{approximate p-value of the test}
#' \item{data.name}{desciption of data}
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
#' @seealso \link{hsic.perm}, \link{hsic.clust}, \link{hsic.gamma}, \link{dcov.test}, \link{kernelCItest}
#' @examples
#' library(energy)
#' set.seed(10)
#' #independence
#' x <- runif(300)
#' y <- runif(300)
#'
#' hsic.gamma(x,y)
#' hsic.perm(x,y)
#' dcov.gamma(x,y)
#' dcov.test(x,y)
#'
#' #uncorelated but not dependent
#' z <- 10*(runif(300)-0.5)
#' w <- z^2 + 10*runif(300)
#'
#' cor(z,w)
#' hsic.gamma(z,w)
#' hsic.perm(z,w)
#' dcov.gamma(z,w)
#' dcov.test(z,w)

dcov.gamma <- function(x, # first variable
                       y, # second variable
                       index=1, #type of distance
                       numCol=100 #number of columns used in incomplete SVD decomposition
){
  n <- length(x)
  m <- length(y)
  if (index < 0 || index > 2) {
    warning("index must be in [0,2), using default index=1")
    index = 1
  }
  if (n != m)
    stop("Sample sizes must agree")
  if (!(all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  H <- diag(n)-matrix(1/n,n,n)
  matx <- as.matrix(dist(x))
  maty <- as.matrix(dist(y))
  ##
  P <- eigs(matx,numCol)
  Q <- eigs(maty,numCol)
  Ux <- P$vectors
  Sx <- diag(P$values)
  Uy <- Q$vectors
  Sy <- diag(Q$values)
  ##
  nV2 <- sum(diag( (H%*%Ux) %*%Sx %*% ((t(Ux)%*%H) %*% (H%*%Uy)) %*%Sy%*% (t(Uy)%*%H) ))/n
  nV2Mean <- mean(matx)*mean(maty)
  nV2Variance <- 2*(n-4)*(n-5)/n/(n-1)/(n-2)/(n-3) * sum(diag( (H%*%Ux) %*%Sx %*% ((t(Ux)%*%H) %*% (H%*%Ux)) %*%Sx%*% (t(Ux)%*%H) )) * sum(diag( (H%*%Uy) %*%Sy %*% ((t(Uy)%*%H) %*% (H%*%Uy)) %*%Sy%*% (t(Uy)%*%H) ))/n^4 * n^2
  alpha <- (nV2Mean)^2/nV2Variance
  beta  <- nV2Variance/(nV2Mean)
  pval <- 1-pgamma(q = nV2, shape = alpha, rate = 1/beta)
  dCov <- sqrt(nV2/n)

  names(dCov) <- "dCov"
  names(nV2) <- "nV^2"
  names(nV2Mean) <-"nV^2 mean"
  names(nV2Variance) <- "nV^2 variance"
  dataname <- paste("index 1, Gamma approximation", sep = "")
  e <- list(method = paste("dCov test of independence", sep = ""),
            statistic = nV2,
            estimate = dCov,
            estimates = c(nV2,nV2Mean,nV2Variance),
            p.value = pval,
            replicates = NULL,
            data.name = dataname)
  class(e) <- "htest"
  return(e)

}

#' Formula for GAM without crossterms
#'
#' Creates a formula for \link{gam} to be used in \link{regrXonS}. For data \eqn{X=(X_1,...X_n,X_{n+1},...,X_m)}, variable to be regressed \eqn{X_{i}}, i=1...n and variables to regress on \eqn{S={X_{n+1},...,X_m}} creates formula \eqn{X_i \sim s(X_{n+1}) + ... + s(X_m)}{X_i ~ s(X_{n+1}) + ... + s(X_m)}.
#'
#' @param target.ind integer, number for the variable to be regressed
#' @param pred.inds integer(s), number(s) for the variable(s) on which we regress
#' @param var.str name of variables used to create formula, default is "x"
#' @importFrom stats formula as.formula
#' @import methods
#' @export
#' @return formula.additive.smooth() returns a formula \eqn{X_i \sim s(X_{n+1}) + ... + s(X_m)}{X_i ~ s(X_{n+1}) + ... + s(X_m)}
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk})
#' @seealso \link{regrXonS}

# frml.additive.smooth <- function(target.ind,pred.inds,var.str="x")
# {
#   as.formula(paste(c(var.str,target.ind," ~ ",paste(paste("s(",var.str,pred.inds,")",sep=""),collapse="+"),sep=""),collapse=""))
# }
frml.additive.smooth <- function(target.ind, pred.inds, var.str = "x", k ) {
  # 这里通过 sprintf 将 k 值传递到 s() 函数中
  smooth_terms <- paste(sprintf("s(%s%s, k = %d)", var.str, pred.inds, k), collapse = " + ")
  
  # 构造公式并返回
  as.formula(paste0(var.str, target.ind, "~", smooth_terms))
}

#' Formula for GAM with crossterms
#'
#' Creates a formula for \link{gam} to be used in \link{regrXonS}. For data \eqn{X=(X_1,...X_n,X_{n+1},...,X_m)}, variable to be regressed \eqn{X_{i}}, i=1...n and variables to regress on \eqn{S={X_{n+1},...,X_m}} creates formula \eqn{X_i \sim s(X_{n+1},...,X_m}{X_i ~ s(X_{n+1},...,X_m}.
#'
#' @param target.ind integer, number for the variable to be regressed
#' @param pred.inds integer(s), number(s) for the variable(s) on which we regress
#' @param var.str name of variables used to create formula, default is "x"
#' @importFrom stats formula as.formula
#' @import methods
#' @export
#' @return formula.full.smooth() returns a formula \eqn{X_i \sim s(X_{n+1},...,X_m)}{X_i ~ s(X_{n+1},...,X_m)}
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk})
#' @seealso \link{regrXonS}

# frml.full.smooth <- function(target.ind,pred.inds,var.str="x")
# {
#   as.formula(paste(c(var.str,target.ind," ~ s(",paste(paste(var.str,pred.inds,sep=""),collapse=","),")",collapse=""),collapse=""))
# }
frml.full.smooth <- function (target.ind, pred.inds, var.str = "x", k) {
    predictors <- paste(paste(var.str, pred.inds, sep = ""), collapse = ",")
    smoothing_term <- paste("s(", predictors, ", k=", k, ")", sep = "")
    formula_string <- paste0(var.str, target.ind, " ~ ", smoothing_term, sep = "")
    as.formula(formula_string)
}

#' HSIC cluster permutation conditional independence test
#'
#' Conditional independence test using HSIC and permutation with clusters.
#'
#' @details Let x and y be two samples of length n. Gram matrices K and L are defined as: \eqn{K_{i,j} = \exp\frac{(x_i-x_j)^2}{\sigma^2}}{K_{i,j} =exp((x_i-x_j)^2/sig^2)}, \eqn{L_{i,j} = \exp\frac{(y_i-y_j)^2}{\sigma^2}}{L_{i,j} =exp((y_i-y_j)^2/sig^2)} and \eqn{M_{i,j} = \exp\frac{(z_i-z_j)^2}{\sigma^2}}{M_{i,j} =exp((z_i-z_j)^2/sig^2)}. \eqn{H_{i,j} = \delta_{i,j} - \frac{1}{n}}{H_{i,j} = delta_{i,j} - 1/n}. Let \eqn{A=HKH}, \eqn{B=HLH} and \eqn{C=HMH}. \eqn{HSIC(X,Y|Z) = \frac{1}{n^2}Tr(AB-2AC(C+\epsilon I)^{-2}CB+AC(C+\epsilon I)^{-2}CBC(C+\epsilon I)^{-2}C)}{HSIC(X,Y|Z) = Tr(AB-2AC(C+\epsilon I)^{-2}CB+AC(C+\epsilon I)^{-2}CBC(C+\epsilon I)^{-2}C)/n^2}. Permutation test clusters Z and then permutes Y in the clusters of Z p times to get \eqn{Y_{(p)}} and calculates \eqn{HSIC(X,Y_{(p)}|Z)}. \eqn{pval = \frac{1(HSIC(X,Y|Z)>HSIC(Z,Y_{(p)}|Z))}{p}}{(HSIC(X,Y|Z)>HSIC(X,Y_{(p)}|Z))/p}.
#' @param x first variable
#' @param y second variable
#' @param z set of variables on which we condition
#' @param sig the with of the Gaussian kernel
#' @param p the number of permutations
#' @param numCluster number of clusters for clustering z
#' @param numCol maximum number of columns that we use for the incomplete Cholesky decomposition
#' @param eps normalization parameter for HSIC cluster test
#' @param paral number of cores used
#'
#' @importFrom kernlab inchol rbfdot
#' @importFrom parallel makeCluster clusterEvalQ parLapply stopCluster
#' @importFrom stats var kmeans pgamma
#' @import methods
#'
#' @references Tillman, R. E., Gretton, A. and Spirtes, P. (2009). Nonlinear directed acyclic structure learning with weakly additive noise model. NIPS 22, Vancouver.
#' @references K. Fukumizu et al. (2007). Kernel Measures of Conditional Dependence. NIPS 20. \url{https://papers.nips.cc/paper/3340-kernel-measures-of-conditional-dependence.pdf}
#
#' @export
#' @return hsic.clust() returns a list with class htest containing
#' \item{method}{description of test}
#' \item{statistic}{observed value of the test statistic}
#' \item{estimate}{HSIC(x,y)}
#' \item{estimates}{a vector: [HSIC(x,y), mean of HSIC(x,y), variance of HSIC(x,y)]}
#' \item{replicates}{replicates of the test statistic}
#' \item{p.value}{approximate p-value of the test}
#' \item{data.name}{desciption of data}
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
#' @seealso \link{hsic.gamma}, \link{hsic.perm}, \link{kernelCItest}
#' @examples
#' library(energy)
#' set.seed(10)
#' # x and y dependent, but independent conditionally on z
#' z <- 10*runif(300)
#' x <- sin(z) + runif(300)
#' y <- cos(z) + runif(300)
#' plot(x,y)
#' hsic.gamma(x,y)
#' hsic.perm(x,y)
#' dcov.test(x,y)
#' hsic.clust(x,y,z)


hsic.clust <- function(                x,             # the first variable for cond ind
                                       y,             # the second variable for cond ind
                                       z,             # the set of variables to wich I condition
                                       sig=1,         # the width of gaussian kernel
                                       p=100,         # number of permutation of 1-alpha level test
                                       numCluster=10, #number of clusters if we use kpc cluster permutation
                                       numCol=50,      # number of column in Incomplete Cholesky Decomposition
                                       eps=0.1,       #normalization parameter
                                       paral=1       # number of cores
                                       ){
  x = as.matrix(x)
  y = as.matrix(y)
  z = as.matrix(z)
  n = nrow(as.matrix(x))
  H = diag(n)-1/n*matrix(rep(x = 1, length.out = n*n), n, n)
  pval.perm=matrix(0,p)
  ###
  tempx = inchol(as.matrix(x), kernel = 'rbfdot', kpar = list(sigma = 1/sig), maxiter = numCol)
  iCHx = tempx@.Data
  tempy = inchol(as.matrix(y), kernel = 'rbfdot', kpar = list(sigma = 1/sig), maxiter = numCol)
  iCHy = tempy@.Data
  tempz = inchol(as.matrix(z), kernel = 'rbfdot', kpar = list(sigma = 1/sig), maxiter = numCol)
  iCHz = tempz@.Data
  a = svd(x = H%*%iCHx)
  Ux = a$u
  Sx = a$d
  ##
  b = svd(x = H%*%iCHy)
  Uy = b$u
  Sy = b$d
  ##
  c = svd(x = H%*%iCHz)
  Uz = c$u
  Sz = c$d
  ##
  if (length(Sx) != 1){ MSx <- diag(Sx^2) }  else {MSx <- as.matrix(Sx^2)}
  if (length(Sy) != 1){ MSy <- diag(Sy^2) }  else {MSy <- as.matrix(Sy^2)}
  if (length(Sz) != 1){ MSz <- diag(Sz^2/(Sz^2+eps)) }  else {MSz <- as.matrix(Sz^2/(Sz^2+eps))}
  ##
  hsic = 1/(n*n)*sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Uy%*%MSy))%*%t(Uy)  -2*Ux%*%((((MSx%*%t(Ux))%*%(Uz%*%MSz))%*%t(Uz))%*%(Uy%*%MSy))%*%t(Uy) + Ux%*%((((MSx%*%t(Ux))%*%(Uz%*%MSz))%*%t(Uz))%*%(Uy%*%((MSy%*%t(Uy))%*%(Uz%*%MSz))))%*%t(Uz)))
  if (paral == 1){
    pval.perm = lapply(X = 1:p, function(m, Ux, Sx, Uy, Sy, Uz, Sz, sig, numCol, numCluster) {
      t = kmeans(z,numCluster)
      perm = c(1:n)
      for (j in 1:numCluster){
        perm[t$cluster==j] = perm[t$cluster==j][sample(sum(t$cluster==j))]
      }
      ##
      Uyp <- Uy[perm,]
      result = 1/(n*n)*sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Uyp%*%MSy))%*%t(Uyp)  -2*Ux%*%((((MSx%*%t(Ux))%*%(Uz%*%MSz))%*%t(Uz))%*%(Uyp%*%MSy))%*%t(Uyp) + Ux%*%((((MSx%*%t(Ux))%*%(Uz%*%MSz))%*%t(Uz))%*%(Uyp%*%((MSy%*%t(Uyp))%*%(Uz%*%MSz))))%*%t(Uz)))
      return (  result ) },  Ux=Ux, Sx=Sx, Uy=Uy, Sy=Sy, Uz=Uz, Sz=Sz, sig, numCol, numCluster=numCluster)
  }
  else {
    cl <- makeCluster(mc <- getOption('cl.cores',paral))
    clusterEvalQ(cl, library(kernlab))
    clusterEvalQ(cl, library(psych))
    #clusterEvalQ(cl, source('kpc_functions.R'))?clust
    pval.perm <- parLapply(cl = cl, X = 1:p, function(m, Ux, Sx, Uy, Sy, Uz, Sz, sig, eps, numCol, numCluster) { m
      t <- kmeans(z,numCluster)
      perm <- c(1:n)
      for (j in 1:numCluster){
        perm[t$cluster==j] <- perm[t$cluster==j][sample(sum(t$cluster==j))]
      }
      ##
      Uyp <- Uy[perm,]
      result <- 1/(n*n)*sum(diag( Ux%*%((diag(Sx^2)%*%t(Ux))%*%(Uyp%*%diag(Sy^2)))%*%t(Uyp)  -2*Ux%*%((((diag(Sx^2)%*%t(Ux))%*%(Uz%*%diag(Sz^2/(Sz^2+eps))))%*%t(Uz))%*%(Uyp%*%diag(Sy^2)))%*%t(Uyp) + Ux%*%((((diag(Sx^2)%*%t(Ux))%*%(Uz%*%diag(Sz^2/(Sz^2+eps))))%*%t(Uz))%*%(Uyp%*%((diag(Sy^2)%*%t(Uyp))%*%(Uz%*%diag(Sz^2/(Sz^2+eps))))))%*%t(Uz)))
      return (  result ) },  Ux=Ux, Sx=Sx, Uy=Uy, Sy=Sy, Uz=Uz, Sz=Sz, sig=sig, eps=eps, numCol=numCol, numCluster=numCluster)
    stopCluster(cl)
  }
  pval.perm2 <- matrix(0,ncol=p)
  for (i in 1:p){pval.perm2[i]<-pval.perm[[i]]}
  pval <- mean(c(pval.perm2,hsic)>=hsic)

  HSMean <- mean(pval.perm2)
  HSVariance <- var(pval.perm2)

  names(hsic) <- "HSIC"
  names(HSMean) <-"HSIC mean"
  names(HSVariance) <- "HSIC variance"
  dataname <- paste("Cluster permutation approximation", sep = "")
  e <- list(method = paste("HSIC test of conditional independence", sep = ""),
            statistic = hsic,
            estimate = hsic,
            estimates = c(hsic,HSMean,HSVariance),
            p.value = pval,
            replicates = pval.perm2,
            data.name = dataname)
  class(e) <- "htest"
  return(e)

}

#' Hilber Schmidt Independence Criterion gamma test
#'
#' Test to check the independence between two variables x and y using HSIC.
#' The hsic.gamma() function, uses Hilbert-Schmidt independence criterion to test for independence between
#' random variables.
#'
#' @details Let x and y be two samples of length n. Gram matrices K and L are defined as: \eqn{K_{i,j} = \exp\frac{(x_i-x_j)^2}{\sigma^2}}{K_{i,j} =exp((x_i-x_j)^2/sig^2)} and \eqn{L_{i,j} = \exp\frac{(y_i-y_j)^2}{\sigma^2}}{L_{i,j} =exp((y_i-y_j)^2/sig^2)}. \eqn{H_{i,j} = \delta_{i,j} - \frac{1}{n}}{H_{i,j} = delta_{i,j} - 1/n}. Let \eqn{A=HKH} and \eqn{B=HLH}, then \eqn{HSIC(x,y)=\frac{1}{n^2}Tr(AB)}{HSIC(x,y)=Tr(AB)/n^2}. Gamma test compares HSIC(x,y) with the \eqn{\alpha}{alpha} quantile of the gamma distribution with mean and variance such as HSIC under independence hypothesis.
#' @param x data of first sample
#' @param y data of second sample
#' @param sig Gaussian kernel width for HSIC tests. Default is 1
#' @param numCol maximum number of columns that we use for the incomplete Cholesky decomposition
#'
#' @importFrom kernlab inchol rbfdot
#' @importFrom stats var kmeans pgamma
#' @import methods
#'
#' @references A. Gretton et al. (2005). Kernel Methods for Measuring Independence. JMLR 6 (2005) 2075-2129.
#'
#' @export
#' @return hsic.gamma() returns a list with class htest containing
#' \item{method}{description of test}
#' \item{statistic}{observed value of the test statistic}
#' \item{estimate}{HSIC(x,y)}
#' \item{estimates}{a vector: [HSIC(x,y), mean of HSIC(x,y), variance of HSIC(x,y)]}
#' \item{replicates}{replicates of the test statistic}
#' \item{p.value}{approximate p-value of the test}
#' \item{data.name}{desciption of data}
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
#' @seealso \link{hsic.perm}, \link{hsic.clust}, \link{kernelCItest}
#' @examples
#' library(energy)
#' set.seed(10)
#' #independence
#' x <- runif(300)
#' y <- runif(300)
#'
#' hsic.gamma(x,y)
#' hsic.perm(x,y)
#' dcov.gamma(x,y)
#' dcov.test(x,y)
#'
#' #uncorelated but not dependent
#' z <- 10*(runif(300)-0.5)
#' w <- z^2 + 10*runif(300)
#'
#' cor(z,w)
#' hsic.gamma(z,w)
#' hsic.perm(z,w)
#' dcov.gamma(z,w)
#' dcov.test(z,w)

hsic.gamma <- function(x, # first variable
                       y, # second variable
                       sig=1, #sigma for the Gaussian kernel
                       numCol=100
                       ){
  n <- length(x)
  H <- diag(n)-matrix(1/n,n,n)
  rbf <- rbfdot(sigma=1/sig)
  tempx <- inchol(as.matrix(x), kernel = rbf, maxiter = numCol)
  iCHx <- tempx@.Data
  a <- svd(x = H%*%iCHx)
  Ux <- a$u
  Sx <- a$d
  ##
  tempy <- inchol(as.matrix(y), kernel = rbf, maxiter = numCol)
  iCHy <- tempy@.Data
  b <- svd(x = H%*%iCHy)
  Uy <- b$u
  Sy <- b$d
  ##
  if (length(Sx) != 1){ MSx <- diag(Sx^2) }  else {MSx <- as.matrix(Sx^2)}
  if (length(Sy) != 1){ MSy <- diag(Sy^2) }  else {MSy <- as.matrix(Sy^2)}
  hsic <- sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Uy%*%MSy))%*%t(Uy) )) /n^2
  K <- crossprod(t(iCHx))
  L <- crossprod(t(iCHy))
  mux <- 1/(n*(n-1))*sum(K-diag(K)*diag(ncol(K)))
  muy <- 1/(n*(n-1))*sum(L-diag(L)*diag(ncol(L)))
  HSMean <- 1/n * (1 + mux*muy - mux - muy )
  HSVariance <- (2*(n-4)*(n-5)/(n*(n-1)*(n-2)*(n-3))) * sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Ux%*%MSx))%*%t(Ux) )) * sum(diag( Uy%*%((MSy%*%t(Uy))%*%(Uy%*%MSy))%*%t(Uy) ))/n^4
  alpha <- (HSMean)^2/HSVariance
  beta  <- HSVariance/(HSMean)
  pval <- 1-pgamma(q = hsic, shape = alpha, rate = 1/beta)

  names(hsic) <- "HSIC"
  names(HSMean) <-"HSIC mean"
  names(HSVariance) <- "HSIC variance"
  dataname <- paste("Gamma approximation", sep = "")
  e <- list(method = paste("HSIC test of independence", sep = ""),
            statistic = hsic,
            estimate = hsic,
            estimates = c(hsic,HSMean,HSVariance),
            p.value = pval,
            replicates = NULL,
            data.name = dataname)
  class(e) <- "htest"
  return(e)

}

#' Hilber Schmidt Independence Criterion permutation test
#'
#' Test to check the independence between two variables x and y using HSIC.
#' The hsic.perm() function, uses Hilbert-Schmidt independence criterion to test for independence between
#' random variables.
#'
#' @details Let x and y be two samples of length n. Gram matrices K and L are defined as: \eqn{K_{i,j} = \exp\frac{(x_i-x_j)^2}{\sigma^2}}{K_{i,j} =exp((x_i-x_j)^2/sig^2)} and \eqn{L_{i,j} = \exp\frac{(y_i-y_j)^2}{\sigma^2}}{L_{i,j} =exp((y_i-y_j)^2/sig^2)}. \eqn{H_{i,j} = \delta_{i,j} - \frac{1}{n}}{H_{i,j} = \delta_{i,j} - 1/n}. Let \eqn{A=HKH} and \eqn{B=HLH}, then \eqn{HSIC(x,y)=\frac{1}{n^2}Tr(AB)}{HSIC(x,y)=Tr(AB)/n^2}. Permutation test permutes y p times to get \eqn{y_{(p)}} and calculates HSIC(x,y_{(p)}). \eqn{pval = \frac{1(HSIC(x,y)>HSIC(x,y_{(p)}))}{p}}{(HSIC(x,y)>HSIC(x,y_{(p)}))/p}.
#' @param x data of first sample
#' @param y data of second sample
#' @param p Number of permutations. Default is 100
#' @param numCol maximum number of columns that we use for the incomplete Cholesky decomposition
#' @param sig Gaussian kernel width for HSIC tests. Default is 1
#'
#' @importFrom kernlab inchol rbfdot
#' @importFrom stats var kmeans pgamma
#' @import methods
#'
#' @references A. Gretton et al. (2005). Kernel Methods for Measuring Independence. JMLR 6 (2005) 2075-2129.
#'
#' @export
#' @return hsic.perm() returns a list with class htest containing
#' \item{method}{description of test}
#' \item{statistic}{observed value of the test statistic}
#' \item{estimate}{HSIC(x,y)}
#' \item{estimates}{a vector: [HSIC(x,y), mean of HSIC(x,y), variance of HSIC(x,y)]}
#' \item{replicates}{replicates of the test statistic}
#' \item{p.value}{approximate p-value of the test}
#' \item{data.name}{desciption of data}
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
#' @seealso \link{hsic.gamma}, \link{hsic.clust}, \link{kernelCItest}
#' @examples
#' library(energy)
#' set.seed(10)
#' #independence
#' x <- runif(300)
#' y <- runif(300)
#'
#' hsic.gamma(x,y)
#' hsic.perm(x,y)
#' dcov.gamma(x,y)
#' dcov.test(x,y)
#'
#' #uncorelated but not dependent
#' z <- 10*(runif(300)-0.5)
#' w <- z^2 + 10*runif(300)
#'
#' cor(z,w)
#' hsic.gamma(z,w)
#' hsic.perm(z,w)
#' dcov.gamma(z,w)
#' dcov.test(z,w)



hsic.perm <- function(   x, # first variable
                         y, # second variable
                         sig=1, #sigma for the Gaussian kernel
                         p = 100, # number of permutation in finding p-value
                         numCol=50 # number of columns for incomplete Cholesky decomposition
                         ){
  n <- length(x)
  H <- diag(n)-matrix(1/n,n,n)
  rbf <- rbfdot(sigma=1/sig)
  tempx <- inchol(as.matrix(x), kernel=rbf, maxiter = numCol)
  iCHx <- tempx@.Data
  a <- svd(x = H%*%iCHx)
  Ux <- a$u
  Sx <- a$d
  ##
  tempy <- inchol(as.matrix(y), kernel=rbf, maxiter = numCol)
  iCHy <- tempy@.Data
  b <- svd(x = H%*%iCHy)
  Uy <- b$u
  Sy <- b$d
  ##
  if (length(Sx) != 1){ MSx <- diag(Sx^2) }  else {MSx <- as.matrix(Sx^2)}
  if (length(Sy) != 1){ MSy <- diag(Sy^2) }  else {MSy <- as.matrix(Sy^2)}
  hsic <- sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Uy%*%MSy))%*%t(Uy) )) /n^2
  pval.perm <- vector(length = p)
  for( i in 1:p)
  {
    perm <- sample(n)
    pval.perm[i] <- sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Uy[perm,]%*%MSy))%*%t(Uy[perm,]) )) /n^2
  }
  pval <- mean(c(pval.perm,hsic)>=hsic)

  HSMean <- mean(pval.perm)
  HSVariance <- var(pval.perm)

  names(hsic) <- "HSIC"
  names(HSMean) <-"HSIC mean"
  names(HSVariance) <- "HSIC variance"
  dataname <- paste("Permutation approximation", sep = "")
  e <- list(method = paste("HSIC test of independence", sep = ""),
            statistic = hsic,
            estimate = hsic,
            estimates = c(hsic,HSMean,HSVariance),
            p.value = pval,
            replicates = pval.perm,
            data.name = dataname)
  class(e) <- "htest"
  return(e)

}
#' Hilber Schmidt Independence Criterion test
#'
#' Test to check the independence between two variables x and y using HSIC.
#' The hsic.test() function, uses Hilbert-Schmidt independence criterion to test for independence between two random variables.
#'
#' @details Let x and y be two samples of length n. Gram matrices K and L are defined as: \eqn{K_{i,j} = \exp\frac{(x_i-x_j)^2}{\sigma^2}}{K_{i,j} =exp((x_i-x_j)^2/sig^2)} and \eqn{L_{i,j} = \exp\frac{(y_i-y_j)^2}{\sigma^2}}{L_{i,j} =exp((y_i-y_j)^2/sig^2)}. \eqn{H_{i,j} = \delta_{i,j} - \frac{1}{n}}{H_{i,j} = delta_{i,j} - 1/n}. Let \eqn{A=HKH} and \eqn{B=HLH}, then \eqn{HSIC(x,y)=\frac{1}{n^2}Tr(AB)}{HSIC(x,y)=Tr(AB)/n^2}.
#' @param x data of first sample
#' @param y data of second sample
#' @param hsic.method method for HSIC test, either gamma test \link{hsic.gamma} or permutation test \link{hsic.perm}
#' @param p number of replicates, if 0
#' @param sig Gaussian kernel width for HSIC. Default is 1
#' @param numCol number of columns in the Incomplete Cholesky Decomposition of Gram matrices. Default is floor(length(x)/10)
#'
#' @references A. Gretton et al. (2005). Kernel Methods for Measuring Independence. JMLR 6 (2005) 2075-2129.
#'
#' @export
#' @return hsic.gamma() returns a list with class htest containing
#' \item{method}{description of test}
#' \item{statistic}{observed value of the test statistic}
#' \item{estimate}{HSIC(x,y)}
#' \item{estimates}{a vector: [HSIC(x,y), mean of HSIC(x,y), variance of HSIC(x,y)]}
#' \item{replicates}{replicates of the test statistic}
#' \item{p.value}{approximate p-value of the test}
#' \item{data.name}{desciption of data}
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
#' @seealso \link{hsic.perm}, \link{hsic.clust}, \link{kernelCItest}
#' @examples
#' library(energy)
#' set.seed(10)
#' #independence
#' x <- runif(300)
#' y <- runif(300)
#'
#' hsic.gamma(x,y)
#' hsic.perm(x,y)
#' dcov.gamma(x,y)
#' dcov.test(x,y)
#'
#' #uncorelated but not dependent
#' z <- 10*(runif(300)-0.5)
#' w <- z^2 + 10*runif(300)
#'
#' cor(z,w)
#' hsic.gamma(z,w)
#' hsic.perm(z,w)
#' dcov.gamma(z,w)
#' dcov.test(z,w)

hsic.test <-  function(x, # first variable
                       y, # second variable
                       p = 0,
                       hsic.method = c("gamma","perm"),
                       sig=1, #sigma for the Gaussian kernel
                       numCol=floor(length(x)/10)
){
  n <- length(x)
  H <- diag(n)-matrix(1/n,n,n)
  rbf <- rbfdot(sigma=1/sig)
  tempx <- inchol(as.matrix(x), kernel = rbf, maxiter = numCol)
  iCHx <- tempx@.Data
  a <- svd(x = H%*%iCHx)
  Ux <- a$u
  Sx <- a$d
  ##
  tempy <- inchol(as.matrix(y), kernel = rbf, maxiter = numCol)
  iCHy <- tempy@.Data
  b <- svd(x = H%*%iCHy)
  Uy <- b$u
  Sy <- b$d
  ##
  if (length(Sx) != 1){ MSx <- diag(Sx^2) }  else {MSx <- as.matrix(Sx^2)}
  if (length(Sy) != 1){ MSy <- diag(Sy^2) }  else {MSy <- as.matrix(Sy^2)}
  hsic <- sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Uy%*%MSy))%*%t(Uy) )) /n^2

  if (p != 0){
    if (hsic.method == "gamma"){
      K <- crossprod(t(iCHx))
      L <- crossprod(t(iCHy))
      mux <- 1/(n*(n-1))*sum(K-diag(K)*diag(ncol(K)))
      muy <- 1/(n*(n-1))*sum(L-diag(L)*diag(ncol(L)))
      HSMean <- 1/n * (1 + mux*muy - mux - muy )
      HSVariance <- (2*(n-4)*(n-5)/(n*(n-1)*(n-2)*(n-3))) * sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Ux%*%MSx))%*%t(Ux) )) * sum(diag( Uy%*%((MSy%*%t(Uy))%*%(Uy%*%MSy))%*%t(Uy) ))/n^4
      alpha <- (HSMean)^2/HSVariance
      beta  <- HSVariance/(HSMean)
      pval <- 1-pgamma(q = hsic, shape = alpha, rate = 1/beta)
    }
    else if(hsic.method == "perm"){
      pval.perm <- vector(length = p)
      for( i in 1:p)
      {
        perm <- sample(n)
        pval.perm[i] <- sum(diag( Ux%*%((MSx%*%t(Ux))%*%(Uy[perm,]%*%MSy))%*%t(Uy[perm,]) )) /n^2
      }
      pval <- mean(pval.perm>hsic)
      HSMean <- mean(pval.perm)
      HSVariance <- var(pval.perm)
    }
  }
  else {pval <- NA}

  names(hsic) <- "HSIC"
  names(HSMean) <-"HSIC mean"
  names(HSVariance) <- "HSIC variance"
  dataname <- paste("Gamma approximation", sep = "")
  e <- list(method = paste("HSIC test of independence", sep = ""),
            statistic = hsic,
            estimate = hsic,
            estimates = c(hsic,HSMean,HSVariance),
            p.value = pval,
            replicates = NULL,
            data.name = dataname)
  class(e) <- "htest"
  return(e)

}

#' Kernel Conditional Independence test
#'
#' Test to check the (conditional) dependence between two variables x and y given a set of variables S, using independence criteria. The kernelCItest() function, uses Distance Covariance or Hilbert-Schmidt Independence Criterion to test for the (conditional) independence between random variables, with an interface that can easily by used in \code{\link{skeleton}}, \code{\link{pc}} or \code{\link{kpc}}.
#' @param x,y,S It is tested, whether x and y are conditionally independent given the subset S of the remaining nodes. x, y, S all are integers, corresponding to variable or node numbers.
#' @param verbose a logical parameter, if TRUE, detailed output is provided.
#' @param suffStat a list of parameters consisting of data, ic.method, p, index, sig, numCol, numCluster, eps, paral
#' @param data numeric matrix witch collumns representing variables and rows representing samples
#' @param ic.method Method for the (conditional) independence test: Distance Covariance (permutation or gamma test), HSIC (permutation or gamma test) or HSIC cluster
#' @param p Number of permutations for Distance Covariance, HSIC permutation and HSIC cluster tests. Default is Distance Covariance
#' @param index Number in (0,2] the power of the distance in the Distance Covariance
#' @param sig Gaussian kernel width for HSIC tests. Default is 1
#' @param numCol Number of columns used in the incomplete Cholesky decomposition. Default is 50
#' @param numCluster Number of clusters for kPC clust algorithm
#' @param eps Normalization parameter for kPC clust. Default is 0.1
#' @param paral Number of cores to use for parallel calculations.
#'
#'
#' @importFrom energy dcov.test
#' @import methods
#'
#' @references G. Szekely, M. Rizzo and N. Bakirov (2007). Measuring and Testing Dependence by Correlation of Distances. The Annals of Statistics 2007, Vol. 35, No. 6, 2769-2794.
#' @references A. Gretton et al. (2005). Kernel Methods for Measuring Independence. JMLR 6 (2005) 2075-2129.
#' @references R. Tillman, A. Gretton and P. Spirtes (2009). Nonlinear directed acyclic structure learning with weakly additive noise model. NIPS 22, Vancouver.
#'
#' @export
#' @return kernelCItest() returns the p-value of the test.
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk}) and Nina Ines Bertille Desgranges
#' @examples
#' set.seed(10)
#' library(pcalg)
#' z <- 10*runif(300)
#' w <- 10*runif(300)
#' x <- sin(z) + runif(300)
#' y <- cos(z) + runif(300)
#'
#' data <- cbind(x,y,z,w)
#'
#' #conditionally independent
#' test1a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="dcc.gamma"))
#' test2a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="dcc.perm"))
#' test3a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="hsic.gamma"))
#' test4a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="hsic.perm"))
#' test5a <- kernelCItest(x=1,y=2,S=c(3),suffStat = list(data=data,ic.method="hsic.clust"))
#' test6a <- gaussCItest( x=1,y=2,S=c(3),suffStat = list(C=cor(data),n=4))
#'
#' test1a
#' test2a
#' test3a
#' test4a
#' test5a
#' test6a
#'
#' #dependent
#' test1b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="dcc.gamma"))
#' test2b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="dcc.perm"))
#' test3b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="hsic.gamma"))
#' test4b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="hsic.perm"))
#' test5b <- kernelCItest(x=1,y=2,S=c(4),suffStat = list(data=data,ic.method="hsic.clust"))
#' test6b <- gaussCItest( x=1,y=2,S=c(4),suffStat = list(C=cor(data),n=4))
#'
#' test1b
#' test2b
#' test3b
#' test4b
#' test5b
#' test6b

kernelCItest <- function (  x,
                            y,
                            S=NULL,
                            suffStat,
                            verbose=FALSE,
                            data,
                            ic.method=NULL,
                            p=NULL,
                            index=NULL,
                            sig=NULL,
                            numCol=NULL,
                            numCluster=NULL,
                            eps=NULL,
                            paral=NULL
){

  x <- as.numeric(x)
  y <- as.numeric(y)

  if(is.null(suffStat$ic.method)) {suffStat$ic.method <- ic.method}
  if(is.null(suffStat$p)) {suffStat$p <- p}
  if(is.null(suffStat$index)) {suffStat$index <- index}
  if(is.null(suffStat$sig)) {suffStat$sig <- sig}
  if(is.null(suffStat$numCol)) {suffStat$numCol <- numCol}
  if(is.null(suffStat$numCluster)) {suffStat$numCluster <- numCluster}
  if(is.null(suffStat$eps)) {suffStat$eps <- eps}
  if(is.null(suffStat$paral)) {suffStat$paral <- paral}


  if (is.null(suffStat$data)){stop("Need to provide data")}
  n <- nrow(suffStat$data)

  if (is.null(suffStat$ic.method)){suffStat$ic.method<-"dcc.perm"; if(verbose){cat("Independence criterion method was not provided, using the default setting for ic.method distance covariance \n")}}

  if (is.null(suffStat$index) & suffStat$ic.method %in% c("dcc.perm","dcc.gamma") ){suffStat$index<-1; if(verbose){cat("The index for the distance covariance test was not provided, using the default setting index=",1,"\n")}}
  if (is.null(suffStat$p) & suffStat$ic.method %in% c("dcc.perm","hsic.perm","hsic.clust") ){suffStat$p<-100; if(verbose){cat("The number of permutations for the permutation test was not provided, using the default setting p=",100,"\n")}}
  if (is.null(suffStat$sig) & suffStat$ic.method %in% c("hsic.gamma","hsic.perm","hsic.clust") ){suffStat$sig<-1; if(verbose){cat("The kernel width sigma for HSIC not provided, using the default setting sig=",1,"\n")}}
  if (is.null(suffStat$numCol) & suffStat$ic.method %in% c("hsic.gamma","hsic.perm","hsic.clust","dcc.gamma") ){suffStat$numCol<-floor(n/10); if(verbose){cat("The number of columns for the Incomplete Cholesky decomposition was not provided, using the default setting numCol=",floor(n/10),"\n")}}
  if (is.null(suffStat$eps) & suffStat$ic.method %in% c("hsic.gamma","hsic.perm","hsic.clust") ){suffStat$eps<-0.1; if(verbose){cat("The normalizing constanst eps for HSIC cluster was not provided, using the default setting eps=",0.1,"\n")}}
  if (is.null(suffStat$numCluster) & suffStat$ic.method %in% c("hsic.clust") ){suffStat$numCluster<-floor(n/10); if(verbose){cat("The number of clusters for the HSIC cluster was not provided, using the default setting numCluster=",floor(n/10),"\n")}}
  if (is.null(suffStat$paral) & suffStat$ic.method %in% c("hsic.clust") ){suffStat$paral<-1; if(verbose){cat("The number of cores to use for parallel computation was not provided, using the default setting paral=",1,"\n")}}

  if (length(S)==0L){
    if (suffStat$ic.method=='dcc.perm'){pval = dcov.test(x=suffStat$data[,x], y=suffStat$data[,y], index = suffStat$index, R = suffStat$p)$p.value}
    else if (suffStat$ic.method=='dcc.gamma'){pval = dcov.gamma(x=suffStat$data[,x], y=suffStat$data[,y], index = suffStat$index, numCol=suffStat$numCol)$p.value}
    else if (suffStat$ic.method=='hsic.perm'){pval = hsic.perm(x=suffStat$data[,x], y=suffStat$data[,y], sig = suffStat$sig, p=suffStat$p, numCol=suffStat$numCol)$p.value}
    else if ((suffStat$ic.method=='hsic.gamma') || (suffStat$ic.method=='hsic.clust')){pval = hsic.gamma(x=suffStat$data[,x], y=suffStat$data[,y], sig=suffStat$sig, numCol=suffStat$numCol)$p.value}
  }
  else{
    if (suffStat$ic.method=='hsic.clust'){pval <- hsic.clust(x=suffStat$data[,x], y=suffStat$data[,y], z=suffStat$data[,S], sig=suffStat$sig, p=suffStat$p, numCluster=suffStat$numCluster, numCol=suffStat$numCol, eps=suffStat$eps, paral=suffStat$paral)$p.value}
    else {
      residuals <- regrXonS(suffStat$data[,c(x,y)], suffStat$data[,S])
      resx <- residuals[,1]
      resy <- residuals[,2]
      if (suffStat$ic.method=='dcc.perm'){pval = dcov.test(resx, resy, index = suffStat$index, R = suffStat$p)$p.value}
      else if (suffStat$ic.method=='dcc.gamma'){pval = dcov.gamma(x=resx, y=resy, index = suffStat$index, numCol=suffStat$numCol)$p.value}
      else if (suffStat$ic.method=='hsic.perm'){pval = hsic.perm(x=resx, y=resy, sig = suffStat$sig, p=suffStat$p, numCol=suffStat$numCol)$p.value}
      else if (suffStat$ic.method=='hsic.gamma'){pval = hsic.gamma(x=resx, y=resy, sig=suffStat$sig, numCol=suffStat$numCol)$p.value}
    }
  }
  return(pval)
}

#' Estimate the WAN-PDAG using the kPC Algorithm
#'
#' Estimates the weakly additive noise partially directed acyclic graph (WAN-PDAG) from observational data, using the kPC algorithm. This is a version of \code{\link{pc}} from pcalg package, that uses HSIC (\link{hsic.gamma}, \link{hsic.perm} or \link{hsic.clust}) or distance covariance (\code{\link{dcov.test}} or \link{dcov.gamma})  independence tests and \link{udag2wanpdag} instead of \code{\link{udag2pdag}} in the last step.
#'
#' @param suffStat a \link{list} of sufficient statistics, containing all necessary elements for the conditional independence decisions in the function indepTest
#' @param indepTest A function for testing conditional independence. It is internally called as indepTest(x,y,S,suffStat), and tests conditional independence of x and y given S. Here, x and y are variables, and S is a (possibly empty) vector of variables (all variables are denoted by their column numbers in the adjacency matrix). suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence. Default is \link{kernelCItest}.
#' @param alpha significance level (number in (0,1) for the individual conditional independence tests.
#' @param labels (optional) character vector of variable (or "node") names. Typically preferred to specifying p.
#' @param p (optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
#' @param verbose If TRUE, detailed output is provided.
#' @param fixedGaps A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
#' @param fixedEdges A logical matrix of dimension p*p. If entry [i,j] or [j,i] (or both) are TRUE, the edge i-j is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
#' @param NAdelete If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
#' @param m.max	Maximal size of the conditioning sets that are considered in the conditional independence tests.
#' @param u2pd String specifying the method for dealing with conflicting information when trying to orient edges (see details below).
#' @param skel.method	Character string specifying method; the default, "stable" provides an order-independent skeleton, see skeleton.
#' @param conservative Logical indicating if the conservative PC is used. In this case, only option u2pd = "relaxed" is supported. Note that therefore the resulting object might not be extendable to a DAG. See details for more information.
#' @param maj.rule	Logical indicating that the triples shall be checked for ambiguity using a majority rule idea, which is less strict than the conservative PC algorithm. For more information, see details.
#' @param solve.confl	If TRUE, the orientation of the v-structures and the orientation rules work with lists for candidate sets and allow bi-directed edges to resolve conflicting edge orientations. In this case, only option u2pd = relaxed is supported. Note, that therefore the resulting object might not be a CPDAG because bi-directed edges might be present. See details for more information.
#' @details For more information: \code{\link{pc}}.
#'
#' @importFrom pcalg skeleton pc.cons.intern pcAlgo
#' @importClassesFrom graph graphNEL
#' @import methods
#'
#' @export
#' @return An object of class "pcAlgo" (see \code{\link{pcAlgo}}) containing an estimate of the equivalence class of the underlying DAG.
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk})
#'
#' @references Tillman, R. E., Gretton, A. and Spirtes, P. (2009). Nonlinear directed acyclic structure learning with weakly additive noise model. NIPS 22, Vancouver.
#' @examples
#' \dontrun{
#' library(pcalg)
#' set.seed(4)
#' n <- 300
#' data <- NULL
#' x1 <- 2*(runif(n)-0.5)
#' x2 <- x1 + runif(n)-0.5
#' x3 <- x1^2 + 0.6*runif(n)
#' x4 <- rnorm(n)
#' x5 <- x3 + x4^2 + 2*runif(n)
#' x6 <- 10*(runif(n)-0.5)
#' x7 <- x6^2 + 5*runif(n)
#' x8 <- 2*x7^2 + 1.5*rnorm(n)
#' x9 <- x7 + 4*runif(n)
#' data <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)
#' true <- matrix(0,9,9)
#' true[c(1),c(2,3)]<-true[c(3,4),5]<-true[c(6),c(7)]<-true[c(7),c(8)]<-true[7,9]<-1
#'
#' pc <- pc(suffStat = list(C = cor(data), n = 9),
#'          indepTest = gaussCItest,
#'          alpha = 0.9,
#'          labels = colnames(data),
#'          u2pd = "relaxed",
#'          skel.method = "stable",
#'          verbose = TRUE)
#' kpc1 <- kpc(suffStat = list(data=data, ic.method="dcc.perm"),
#'             indepTest = kernelCItest,
#'             alpha = 0.1,
#'             labels = colnames(data),
#'             u2pd = "relaxed",
#'             skel.method = "stable",
#'             verbose = TRUE)
#' kpc2 <- kpc(suffStat = list(data=data, ic.method="hsic.gamma"),
#'             indepTest = kernelCItest,
#'             alpha = 0.1,
#'             labels = colnames(data),
#'             u2pd = "relaxed",
#'             skel.method = "stable",
#'             verbose = TRUE)
#' kpc3 <- kpc(suffStat = list(data=data, ic.method="hsic.perm"),
#'             indepTest = kernelCItest,
#'             alpha = 0.1,
#'             labels = colnames(data),
#'             u2pd = "relaxed",
#'             skel.method = "stable",
#'             verbose = TRUE)
#' kpc4 <- kpc(suffStat = list(data=data, ic.method="hsic.clust"),
#'             indepTest = kernelCItest,
#'             alpha = 0.1,
#'             labels = colnames(data),
#'             u2pd = "relaxed",
#'             skel.method = "stable",
#'             verbose = TRUE)
#'
#' if (require(Rgraphviz)) {
#'  par(mfrow=c(2,3))
#'  plot(pc,main="pc")
#'  plot(kpc1,main="dpc.perm")
#'  plot(kpc2,main="kpc.gamma")
#'  plot(kpc3,main="kpc.perm")
#'  plot(kpc4,main="kpc.clust")
#'  plot(as(true,"graphNEL"),main="True DAG")
#' }
#' }

kpc <- function (suffStat, indepTest, alpha, labels, p,
                 fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
                 u2pd = c("relaxed","rand", "retry"),
                 skel.method = c("stable", "original","stable.fast"),
                 conservative = FALSE, maj.rule = FALSE,
                 solve.confl = FALSE, verbose = FALSE)
{
  cl <- match.call()
  if (!missing(p))
    stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
                1, p >= 2)
  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    }
    else if (p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }
  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if (u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }
  if (conservative && maj.rule)
    stop("Choose either conservative PC or majority rule PC!")
  skel <- skeleton(suffStat, indepTest, alpha, labels = labels,
                   method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, verbose = verbose)
  skel@call <- cl
  if (!conservative && !maj.rule) {
    wanpdag <- udag2wanpdag(gInput = skel, suffStat=suffStat, indepTest = indepTest, alpha = alpha, verbose = verbose, solve.confl = solve.confl)
  }
  else {
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha, version.unf = c(2, 1), maj.rule = maj.rule, verbose = verbose)
    print(pc.$sk)
    print(pc.$unfTripl)
    wanpdag <- udag2wanpdag(pc.$sk, suffStat=suffStat, indepTest = indepTest, alpha = alpha, verbose = verbose, unfVect = pc.$unfTripl, solve.confl = solve.confl)
  }
  wanpdag
}

#' Regress set of variables on its parents
#'
#' Uses the generalised additive model \link{gam} to non-linearly and non-parametrically regress set of variables X on a set of variables S and returns residuals of X.
#'
#' @param X numeric matrix, set of variables to be regressed. Each column represents separate variable
#' @param S numeric matrix, set of variables we will regress on. Each column represents separate variable
#' @details If the number of variables in S is \eqn{\leq 5}{<= 5} we use \link{frml.full.smooth} as formula for \link{gam} to regress X on S, otherwise we use \link{frml.additive.smooth}.
#' @importFrom mgcv gam
#' @import methods
#' @export
#' @return regrXonS() returns the residuals of X regressed on S.
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk})
#' @seealso \link{kernelCItest}
#' @examples
#' set.seed(10)
#' library(energy)
#' z <- 10*runif(300)
#' w <- 10*runif(300)
#' x <- sin(z) + runif(300)
#' y <- cos(z) + runif(300)
#' data <- cbind(x,y,z,w)
#'
#' hsic.gamma(x,y)
#' hsic.perm(x,y)
#' dcov.test(x,y)
#'
#' resid <- regrXonS(cbind(x,y),cbind(z,w))
#'
#' hsic.gamma(resid[,1],resid[,2])
#' hsic.perm(resid[,1],resid[,2])
#' dcov.test(resid[,1],resid[,2])

regrXonS <- function(    X, # set of variables to be regressed
                         S  # set of variables on which we regress
                         ){

  #library(mgcv)

  X <- as.matrix(X)
  S <- as.matrix(S)

  n1 <- ncol(X)
  n2 <- ncol(S)

  if (n2 > 2){
    formulaXS <- frml.additive.smooth
  }else{
    formulaXS <- frml.full.smooth
  }
  data <- data.frame(cbind(X, S))
  colnames(data) <- paste("x",1:(n1+n2),sep="")
  resX <- matrix(0,nrow=nrow(X),ncol=ncol(X))

  k <- max(6, min(15, round(1/4 * min(unlist(apply(data, 2, function(x) { length(unique(x)) }))))))
  for (i in 1:n1){
    formx <- formulaXS(i,(n1+1):(n1+n2),k=k)
    resX[,i] <- gam(formx,data=data)$residuals
  }
  return(resX)

}

#' Check if variable can be regressed to independence on its parents
#'
#' Uses the generalised additive model \link{gam} to non-linearly and non-parametrically regress variable V on its parents and set of variables S.
#'
#' @param G adjacency matrix, for the graph
#' @param V integer, node which we regress
#' @param S integer(s), set we regress on
#' @param suffStat sufficient statistics to perform the independence test \link{kernelCItest}
#' @param indepTest independence test to check for dependence between residuals of V and S
#' @param alpha numeric cutoff for significance level of individual partial correlation tests
#' @import methods
#' @export
#' @return regrVonPS() returns the number of p-values smaller than the cutoff, i.e 0 means residuals of V are independent of all variables in S
#' @author Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk})
#'

regrVonPS <- function(         G,           # adj matrix
                               V,           # the node I'm investigating
                               S,           # some ( not always all ) nodes to which v is undirected connected
                               suffStat,    # sufficient statistics
                               indepTest=kernelCItest,   # independence test to be used to check independence
                               alpha=0.2       # test level
                               ){

  parentsV <- which(G[,V]==1 & G[V,]==0, arr.ind = T)
  pvalvec <- matrix(0,ncol=length(S))

  test.data <- suffStat$data[,c(V,S)]
  residV <- regrXonS(suffStat$data[,V],suffStat$data[,c(S,parentsV)])

  test.suffStat <- suffStat
  test.suffStat$data <- cbind(residV,suffStat$data[,S])

  for (i in 1:length(S)){
    pvalvec[i] <- indepTest(x=1, y=(i+1), suffStat = test.suffStat)
  }
  output <- sum(pvalvec<alpha)

}

#' Last kPC Algorithm Step: Extend Object with Skeleton to Completed PDAG
#'
#' This function performs the last (generalised transitive) step in the \link{kpc} algorithm. It transforms an object of the class "pcAlgo" containing a skeleton and corresponding conditional independence information into a weakly additive noise directed acyclic graph (CPDAG). The functions first determine the v-structures in the collider step, and then performs the Generalised Transitive Step as described in Tillman et al (2009) to orient as many of the remaining edges as possible.
#'
#'
#' @param gInput "pcAlgo"-object containing skeleton and conditional indepedence information.
#' @param suffStat a list of sufficient statistics, containing all necessary elements for the conditional independence decisions in the function indepTest.
#' @param indepTest A function for testing conditional independence. It is internally called as indepTest(x,y,S,suffStat). Default is \link{kernelCItest}.
#' @param alpha significance level (number in (0,1) for the individual conditional independence tests.
#' @param verbose 0: No output; 1: Details
#' @param unfVect vector containing numbers that encode ambiguous triples (as returned by pc.cons.intern). This is needed in the conservative and majority rule PC algorithms.
#' @param solve.confl	if TRUE, the orientation of the v-structures and the orientation rules work with lists for candidate sets and allow bi-directed edges to resolve conflicting edge orientations. Note that therefore the resulting object is order-independent but might not be a PDAG because bi-directed edges can be present.
#' @param orientCollider if TRUE, collider are oriented.
#' @param rules Array of length 3 containing TRUE or FALSE for each rule. TRUE in position i means that rule i (Ri) will be applied. By default, all rules are used.gInput
#' @importFrom utils combn
#' @importFrom pcalg skeleton triple2numb
#' @importMethodsFrom graph numEdges
#' @import methods
#' @export
#' @details First we perform a collider step, that is orienting triples a-b-c as a->b<-c iff b is not in separating set of a and c. Then we orient edges a-S as a->S if b_r is independent of a set S, where b_r are the residuals of b non parametrically regressed on S and parents of b and none of the edges S_i-a can be oriented as S_i->a, that is residuals S_i_r would be independent of a.
#'
#' @return An oriented object of class "pcAlgo".
#'
#' @references Tillman, R. E., Gretton, A. and Spirtes, P. (2009). Nonlinear directed acyclic structure learning with weakly additive noise model. NIPS 22, Vancouver.
## simulate data
#' @examples
#' \dontrun{
#' library(pcalg)
#' set.seed(4)
#' n <- 300
#' data <- NULL
#' x1 <- 2*(runif(n)-0.5)
#' x2 <- x1 + runif(n)-0.5
#' x3 <- x1^2 + 0.6*runif(n)
#' x4 <- rnorm(n)
#' x5 <- x3 + x4^2 + 2*runif(n)
#' x6 <- 10*(runif(n)-0.5)
#' x7 <- x6^2 + 10*runif(n)
#' x8 <- 2*x7^2 + rnorm(n)
#' x9 <- x7 + 5*runif(n)
#' data <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)
#' true <- matrix(0,9,9)
#' true[c(1),c(2,3)]<-true[c(3,4),5]<-true[c(6),c(7)]<-true[c(7),c(8)]<-true[7,9]<-1

#' ## estimate skeleton
#' resU1 <- skeleton(suffStat = list(data=data, ic.method="dcc.perm", p=200),
#'                   indepTest = kernelCItest,
#'                   verbose = TRUE, alpha = 0.1, p=9)
#'
#' resU2 <- skeleton(suffStat = list(data=data, ic.method="hsic.gamma",
#'                              sig=1, numCol = 50),
#'                   indepTest = kernelCItest,
#'                   verbose = TRUE, alpha = 0.1, p=9)
#'
#' resU3 <- skeleton(suffStat = list(data=data, ic.method="hsic.perm",
#'                              sig=1, numCol = 50, p=200),
#'                   indepTest = kernelCItest,
#'                   verbose = TRUE, alpha = 0.1, p=9)
#'
#' resU4 <- skeleton(suffStat = list(data=data, ic.method="hsic.clust",
#'                              p=200, sig=1, numCluster=100, numCol = 50,
#'                              eps = 0.1, paral = 1),
#'                   indepTest = kernelCItest,
#'                   verbose = TRUE, alpha = 0.1, p=9)
#'
#' resU5 <- skeleton(suffStat = list(C = cor(data), n = n),
#'                   indepTest = gaussCItest,
#'                   verbose = TRUE, alpha = 0.1, p=9)
#'
#' if (require(Rgraphviz)) {
#'  par(mfrow=c(2,3))
#'  plot(resU1,main="dpc")
#'  plot(resU2,main="kpc-resid-gamma")
#'  plot(resU3,main="kpc-resid-perm")
#'  plot(resU4,main="kpc-clust")
#'  plot(resU5,main="pc")
#'  plot(as(true,"graphNEL"),main="True DAG")
#' }
#'
#' ## orient edges using three different methods
#' resD1 <- udag2wanpdag(gInput = resU1,
#'                       suffStat = list(data=data, ic.method="dcc.perm", sig=1, numCol = 50, p=200),
#'                       indepTest = kernelCItest,
#'                       verbose = TRUE, alpha = 0.1)
#' resD2 <- udag2wanpdag(gInput = resU1,
#'                       suffStat = list(data=data, ic.method="hsic.gamma", sig=1, numCol = 50),
#'                       indepTest = kernelCItest,
#'                       verbose = TRUE, alpha = 0.1)
#' resD3 <- udag2wanpdag(gInput = resU1,
#'                       suffStat = list(data=data, ic.method="hsic.perm", sig=1, numCol = 50, p=200),
#'                       indepTest = kernelCItest,
#'                       verbose = TRUE, alpha = 0.1)
#' resD4 <- udag2pdagRelaxed(gInput = resU1, verbose = T)
#' if (require(Rgraphviz)) {
#'  par(mfrow=c(2,3))
#'  plot(resD1,main="dpc")
#'  plot(resD2,main="kpc-resid-gamma")
#'  plot(resD3,main="kpc-resid-perm")
#'  plot(resD4,main="pc")
#'  plot(as(true,"graphNEL"),main="True DAG")
#' }
#' }

udag2wanpdag <- function (       gInput,
                                 suffStat,
                                 indepTest=kernelCItest,
                                 alpha = 0.2,
                                 verbose = FALSE,
                                 unfVect = NULL,
                                 solve.confl = FALSE,
                                 orientCollider = TRUE,
                                 rules = rep(TRUE, 3)
                                 ){
  rule1 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    ind <- which(pdag == 1 & t(pdag) == 0, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(pdag[b, ] == 1 & pdag[, b] == 1 & pdag[a, ] == 0 & pdag[, a] == 0)
      if (length(isC) > 0) {
        for (ii in seq_along(isC)) {
          c <- isC[ii]
          if (!solve.confl | (pdag[b, c] == 1 & pdag[c, b] == 1)) {
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) && !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pdag[b, c] <- 1
                pdag[c, b] <- 0
              }
            }
            else {
              pdag[b, c] <- 1
              pdag[c, b] <- 0
            }
            if (verbose)
              cat("\nRule 1':", a, "->", b, " and ",
                  b, "-", c, " where ", a, " and ", c,
                  " not connected and ", a, b, c, " faithful triple: ",
                  b, "->", c, "\n")
          }
          else if (pdag[b, c] == 0 & pdag[c, b] == 1) {
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) && !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pdag[b, c] <- 2
                pdag[c, b] <- 2
                if (verbose)
                  cat("\nRule 1':", a, "->", b, "<-",
                      c, " but ", b, "->", c, "also possible and",
                      a, b, c, " faithful triple: ", a,
                      "->", b, "<->", c, "\n")
              }
            }
            else {
              pdag[b, c] <- 2
              pdag[c, b] <- 2
              if (verbose)
                cat("\nRule 1':", a, "->", b, "<-", c,
                    " but ", b, "->", c, "also possible and",
                    a, b, c, " faithful triple: ", a, "->",
                    b, "<->", c, "\n")
            }
          }
        }
      }
      if (!solve.confl)
        pdag <- pdag
    }
    pdag
  }
  ################################################################################
  rule2 <- function(pdag, solve.confl = FALSE) {
    ind <- which(pdag == 1 & t(pdag) == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(pdag[a, ] == 1 & pdag[,a] == 0 & pdag[, b] == 1 & pdag[b,] == 0)
      for (ii in seq_along(isC)) {
        c <- isC[ii]
        if (!solve.confl | (pdag[a, b] == 1 & pdag[b,a] == 1)) {
          pdag[a, b] <- 1
          pdag[b, a] <- 0
          if (verbose)
            cat("\nRule 2: Chain ", a, "->", c, "->",b, ":", a, "->", b, "\n")
        }
        else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
          pdag[a, b] <- 2
          pdag[b, a] <- 2
          if (verbose)
            cat("\nRule 2: Chain ", a, "->", c, "->",b, ":", a, "<->", b, "\n")
        }
      }
      if (!solve.confl)
        pdag <- pdag
    }
    pdag
  }
  ################################################################################
  rule3 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    ind <- which(pdag == 1 & t(pdag) == 1,
                 arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      c <- which(pdag[a, ] == 1 & pdag[,a] == 1 & pdag[, b] == 1 & pdag[b,] == 0)
      if (length(c) >= 2) {
        cmb.C <- combn(c, 2)
        cC1 <- cmb.C[1, ]
        cC2 <- cmb.C[2, ]
        for (j in seq_along(cC1)) {
          c1 <- cC1[j]
          c2 <- cC2[j]
          if (pdag[c1, c2] == 0 && pdag[c2, c1] == 0) {
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, c1,a, c2), na.rm = TRUE) && !any(unfVect == triple2numb(p, c2, a, c1), na.rm = TRUE)) {
                if (!solve.confl | (pdag[a, b] == 1 &
                                    pdag[b, a] == 1)) {
                  pdag[a, b] <- 1
                  pdag[b, a] <- 0
                  if (!solve.confl)
                    pdag <- pdag
                  if (verbose)
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ",a, "->", b, "\n")
                  break
                }
                else if (pdag[a, b] == 0 & pdag[b, a] ==
                         1) {
                  pdag[a, b] <- pdag[b, a] <- 2
                  if (verbose)
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ",a, "<->", b, "\n")
                  break
                }
              }
            }
            else {
              if (!solve.confl | (pdag[a, b] == 1 & pdag[b, a] == 1)) {
                pdag[a, b] <- 1
                pdag[b, a] <- 0
                if (!solve.confl)
                  pdag <- pdag
                if (verbose)
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ",a, "->", b, "\n")
                break
              }
              else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
                pdag[a, b] <- pdag[b, a] <- 2
                if (verbose)
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ",a, "<->", b, "\n")
                break
              }
            }
          }
        }
      }
    }
    pdag
  }
  ################################################################################
  orientConflictCollider <- function(pdag, x, y, z) {
    if (pdag[x, y] == 1) {
      pdag[y, x] <- 0
    }
    else {
      pdag[x, y] <- pdag[y, x] <- 2
    }
    if (pdag[z, y] == 1) {
      pdag[y, z] <- 0
    }
    else {
      pdag[z, y] <- pdag[y, z] <- 2
    }
    pdag
  }
  ################################################################################
  checkImmor <- function(pdag,
                         V,
                         S,
                         ...){
    output <- TRUE
    udag <- pmin(pdag+t(pdag),1)
    parentsV <- which(pdag[,V]==1 & pdag[V,]==0, arr.ind = T)
    if (length(S)>1){
      test.dag1 <- udag[S,S] + diag(length(S))
      if (length(which(test.dag1==0 & t(test.dag1)==0 , arr.ind=T)) > 0) {output <- FALSE}
    }
    test.dag2 <- udag[c(parentsV),S]
    if (sum(test.dag2) < (length(S) * length(parentsV)) ) {output <- FALSE}
    output

  }
  ################################################################################

  if (numEdges(gInput@graph) == 0)
    return(gInput)
  g <- as(gInput@graph, "matrix")
  p <- nrow(g)
  pdag <- g
  if (orientCollider) {
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x)
      for (z in allZ) {
        if (g[x, z] == 0 && !((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]]))) {
          if (length(unfVect)  == 0) {
            if (!solve.confl) {
              pdag[x, y] <- pdag[z, y] <- 1
              pdag[y, x] <- pdag[y, z] <- 0
            }
            else {
              pdag <- orientConflictCollider(pdag, x, y, z)
            }
          }
          else {
            if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) && !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
              if (!solve.confl) {
                pdag[x, y] <- pdag[z, y] <- 1
                pdag[y, x] <- pdag[y, z] <- 0
              }
              else {
                pdag <- orientConflictCollider(pdag, x, y, z)
              }
            }
          }
        }
      }
    }
  }

  s <- 1
  n <- ncol(suffStat$data)

  undi.nbhd <- list()
  ind <- which(pdag==1 & t(pdag)==1, arr.ind = T)
  colnames(ind)<-NULL
  for (i in 1:n){
    undi.nbhd[[i]] <- ind[ind[,1]==i,2]
  }
  size.undi.nbhd <- sapply(undi.nbhd,length)
  nbhd.updt <- rep(0,n)

  while (max(size.undi.nbhd) >= s){
    if(verbose) cat("Checking undirected neighbourhoods of size ",s,"\n")
    for (i in 1:n){
      V <- i
      if(verbose) cat("V: ",V,". Size of neighbourhood of V: ",size.undi.nbhd[V],"\n")
      s2 <- s
      if ( (size.undi.nbhd[V]==s2) | ( size.undi.nbhd[V]<s2 & nbhd.updt[V]) ) {
        s2 <- size.undi.nbhd[V]
        while (s2 > 0){
          if(verbose) cat("Nbhd of V: ",undi.nbhd[[V]],"\n")
          if (length(undi.nbhd[[V]])==1){sub.nbhd.V <- as.matrix(undi.nbhd[[V]])}
          else {sub.nbhd.V <- as.matrix(combn(undi.nbhd[[V]],s2),nrow=s2)}
          for (j in 1:ncol(sub.nbhd.V)){
            to.update <- TRUE
            if (checkImmor(pdag,V,sub.nbhd.V[,j])){
              S <- sub.nbhd.V[,j]
              parentsV <- which(pdag[,V]==1 & pdag[V,]==0, arr.ind = T)
              pval1 <- regrVonPS(G=pdag, V=V, S=S, suffStat=suffStat, indepTest=indepTest, alpha=alpha)
              if (pval1 > 0) {
                if(verbose) cat("Reject orientation",S,"->",V,", as",V,"cannot be regressed to independence on its parents",parentsV,"and subset",S,".\n")
                to.update <- FALSE
              }
              else {
                if(verbose) cat(V,"can be regressed to independence on its parents",parentsV,"and subset",S,".\n")
                for (k in 1:length(S)){
                  W <- S[k]
                  if(verbose) cat("W: ",W,"\n")
                  s3 <- size.undi.nbhd[W]
                  if(verbose) cat("Nbhd of W: ",undi.nbhd[[W]],"\n")
                  to.continue <- TRUE
                  while (to.continue & s3 > 0){
                    if (length(undi.nbhd[[W]])==1){sub.nbhd.W <- as.matrix(undi.nbhd[[W]])}
                    else {sub.nbhd.W <- as.matrix(combn(undi.nbhd[[W]],s3),nrow=s3)}
                    for (j2 in 1:ncol(sub.nbhd.W)){
                      if (checkImmor(pdag,W,sub.nbhd.W[,j2])){
                        S2 <- sub.nbhd.W[,j2]
                        if (sum(S2==V)==1){
                          parentsV <- which(pdag[,W]==1 & pdag[W,]==0, arr.ind = T)
                          pval2 <- regrVonPS(G=pdag, V=W, S=S2, suffStat=suffStat, indepTest=indepTest, alpha=alpha)
                          if (pval2 == 0) {
                            if(verbose) cat("Reject orientation",S,"->",V,", as",W,"can be regressed to independence on its parents",parentsV,"and subset",S2,".\n")
                            to.update <- FALSE
                            to.continue <- FALSE
                          }
                          else {if(verbose) cat(W,"cannot be regressed to independence on its parents",parentsV,"and subset",S2,".\n")}
                        }
                      }
                      else {if(verbose) cat("Reject orientation",sub.nbhd.W[,j2],"->",W,", as it would create immorality.\n")}
                    }
                    s3 <- s3 - 1
                  }

                }
              }
              if (to.update){
                if(verbose) cat("Accept orientation",S,"->",V,".\n")
                pdag[V,S] <- 0
                pdag[S,V] <- 1
                if ( length(setdiff(undi.nbhd[[V]],S)) != 0 ){
                  pdag[setdiff(undi.nbhd[[V]],S),V] <- 0
                  pdag[V,setdiff(undi.nbhd[[V]],S)] <- 1
                }
                s2 <- 0
                repeat {
                  old_pdag <- pdag
                  if (rules[1]) {
                    pdag <- rule1(pdag, solve.confl = solve.confl, unfVect = unfVect)
                  }
                  if (rules[2]) {
                    pdag <- rule2(pdag, solve.confl = solve.confl)
                  }
                  if (rules[3]) {
                    pdag <- rule3(pdag, solve.confl = solve.confl, unfVect = unfVect)
                  }
                  if (all(pdag == old_pdag))
                    break
                }
                undi.nbhd <- list()
                ind <- which(pdag==1 & t(pdag)==1, arr.ind = T)
                colnames(ind)<-NULL
                for (i in 1:n){
                  undi.nbhd[[i]] <- ind[ind[,1]==i,2]
                }
                nbhd.updt[V] <- TRUE
                size.undi.nbhd <- sapply(undi.nbhd,length)
              }
            }
          }
          s2 <- s2 - 1
        }
      }
    }
    s <- s + 1
  }
  gInput@graph <- as(pdag, "graphNEL")
  gInput
}
