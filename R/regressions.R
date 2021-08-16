# center at sample mean
center_colmeans <- function(x){
  xcenter <- colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# Input check
Input_check <- function(Y, A, B, X, pi, q){
  if(!is.numeric(A)||!is.vector(A)){
    stop("Please enter a numeric vector for treatment assignments (A)!")
  }
  if(!is.numeric(B)||!is.vector(B)){
    stop("Please enter a numeric vector for stratum labels (B)!")
  }
  if(!is.numeric(Y)||!is.vector(Y)){
    stop("Please enter a numeric vector for observed outcomes (Y)!")
  }
  if(!is.numeric(X)||!is.matrix(X)){
    stop("Please enter a numeric matrix for additional covariates!")
  }
  if(is.null(X)){
    if(length(A)!=length(B)||length(A)!=length(Y)){
      stop("A, B, Y should have same lengths!")
    }
  }
  else{
    if(length(A)!=length(B)||length(A)!=length(Y)||length(A)!=nrow(X)){
      stop("Lenthgs of A, B, Y and the number of rows in X should be the same!")
    }
  }
  if(pi <= 0){
    stop("Please enter valid (positive) value for pi!")
  }
  if(q < 0){
    stop("Please enter valid (non-negative) value for q!")
  }
}

#'Difference in means
#'
#'Estimating and inferring the treatment effect based on difference in means.
#'
#'Estimating and inferring the treatment effect based on difference in means. It
#'implements the methods as described in Sections 3.1 and 4.1, Ma et al. (2020).
#'
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param B a numeric vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X an (optional) numeric design matrix containing additional covariates
#'  used in the regression.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations. Detailed information can be found in Section 2, Ma et
#'  al.(2020).
#'@param conf.level confidence level of the interval. Default is 0.95.
#'
#'@return A list of class \code{"htest"} containing the following components:
#'\item{statistic}{the value of the t-statistic.} \item{p.value}{the p-value for
#'the test.} \item{conf.int}{a confidence interval under chosen level
#'\code{conf.level} for the difference in treatment effect between treatment
#'group and control group.} \item{estimate}{estimated treatment effect
#'difference between treatment group and control group.} \item{method}{a
#'character string indicating what type of regression was performed.}
#'
#'@references Ma, W., Tu, F., & Liu, H. (2020). \emph{Regression analysis for
#'  covariate-adaptive randomization: A robust and efficient inference
#'  perspective}. arXiv preprint arXiv:2009.02287.
#'
#'@examples
#'#The code replicates the simulation setting of Model 2 in Section 6, Ma et al. (2020).
#'n <- 1000
#'pi <- 0.5
#'q <- pi*(1-pi)
#'alphavec <- c(5,8,3,12)
#'m2e<-function(x){
#'  6*exp(x)*x*(1-x)
#'}
#'X1 <- rgamma(n,2)
#'X2 <- sample(c(1,2,3),n,replace = TRUE, prob = c(0.3,0.6,0.1))
#'X3 <- rpois(n,3)
#'X4 <- rbeta(n,2,2)
#'X1_S <- rep(1,n)
#'X1_S[which(X1 >= 2.5)] <- 2
#'B <- as.numeric(interaction(X1_S,X2))
#'X <- cbind(X1,X3)
#'A <- sample(c(0,1),n,replace=TRUE,prob=c(1-pi,pi))
#'Y0 <- alphavec[2]*3.6-alphavec[1]*2-alphavec[4]*integrate(m2e,lower = 0,upper = 1)$value+
#'      alphavec[1]*X1+log(alphavec[3]*X1*log(X3+1)+1)+alphavec[4]*exp(X4)+rnorm(n,sd = 2)
#'Y1 <- alphavec[2]*X2^2+log(alphavec[3]*X1*log(X3+1)+1)+rnorm(n,sd = 1)
#'Y <- Y0*(1-A)+Y1*A
#'tau.diff(Y, A, B, X, pi, q)
#'@export
tau.diff<-function(Y, A, B, X = NULL, pi, q, conf.level = 0.95){
  Input_check(Y, A, B, X, pi, q)
  strt <- unique(B)
  strt_num <- length(strt)
  if(is.null(X)){
    l <- stats::lm(Y~A)
    estimate <- l$coefficients[2]
    sq <- rep(q, strt_num)
    stderr <- DME_var(A, B, Y, sq, strt, strt_num, pi)
  }
  else{
    p <- ncol(X)
    formtext <- "Y~A"
    for(i in 1:p){
      eval(parse(text = paste("x", i, "= X[,", i, "]", sep = "")))
      formtext <- paste(formtext, "+ x", i, sep = "")
    }
    l <- eval(parse(text = paste("stats::lm(", formtext,")",sep = "")))
    estimate <- l$coefficients[2]
    sq <- rep(q, strt_num)
    stderr <- SR_var(A,B,Y,X,sq,strt,strt_num,pi)
  }
  testmethod<-"Difference in Means"
  tstat <- estimate/stderr
  pval <- dplyr::if_else(stats::pnorm(tstat)<(1-stats::pnorm(tstat)),stats::pnorm(tstat),1-stats::pnorm(tstat))*2
  cint <- c(estimate + stderr*stats::qnorm((1-conf.level)/2),estimate - stderr*stats::qnorm((1-conf.level)/2))
  attr(cint,"conf.level") <- conf.level
  names(tstat) <- "t"
  names(estimate) <- "difference in treatment effect"
  rval<-list(statistic = tstat, p.value = pval, conf.int = cint,
             estimate = estimate, stderr = stderr,
             method = testmethod)
  class(rval) <- "htest"
  return(rval)
}

#'Regression without interaction
#'
#'Estimating and inferring the treatment effect based on regression without
#'interaction.
#'
#'Estimating and inferring the treatment effect based on regression without
#'interaction. It implements the methods as described in Sections 3.2 and 4.2,
#'Ma et al. (2020).
#'
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param B a numeric vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X an (optional) numeric design matrix containing additional covariates
#'  used in the regression.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param q a numeric value indicating the balance level of covariate-adaptive
#'  randomizations. Detailed information can be found in Section 2, Ma et
#'  al.(2020).
#'@param conf.level confidence level of the interval. Default is 0.95.
#'
#'@return A list of class \code{"htest"} containing the following components:
#'\item{statistic}{the value of the t-statistic.} \item{p.value}{the p-value for
#'the test.} \item{conf.int}{a confidence interval under chosen level
#'\code{conf.level} for the difference in treatment effect between treatment
#'group and control group.} \item{estimate}{estimated treatment effect
#'difference between treatment group and control group.} \item{method}{a
#'character string indicating what type of regression was performed.}
#'
#'@references Ma, W., Tu, F., & Liu, H. (2020). \emph{Regression analysis for
#'  covariate-adaptive randomization: A robust and efficient inference
#'  perspective}. arXiv preprint arXiv:2009.02287.
#'
#'@examples
#'#The code replicates the simulation setting of Model 2 in Section 6, Ma et al. (2020).
#'n <- 1000
#'pi <- 0.5
#'q <- pi*(1-pi)
#'alphavec <- c(5,8,3,12)
#'m2e<-function(x){
#'  6*exp(x)*x*(1-x)
#'}
#'X1 <- rgamma(n,2)
#'X2 <- sample(c(1,2,3),n,replace = TRUE, prob = c(0.3,0.6,0.1))
#'X3 <- rpois(n,3)
#'X4 <- rbeta(n,2,2)
#'X1_S <- rep(1,n)
#'X1_S[which(X1 >= 2.5)] <- 2
#'B <- as.numeric(interaction(X1_S,X2))
#'X <- cbind(X1,X3)
#'A <- sample(c(0,1),n,replace=TRUE,prob=c(1-pi,pi))
#'Y0 <- alphavec[2]*3.6-alphavec[1]*2-alphavec[4]*integrate(m2e,lower = 0,upper = 1)$value +
#'      alphavec[1]*X1+log(alphavec[3]*X1*log(X3+1)+1)+alphavec[4]*exp(X4)+rnorm(n,sd = 2)
#'Y1 <- alphavec[2]*X2^2+log(alphavec[3]*X1*log(X3+1)+1)+rnorm(n,sd = 1)
#'Y <- Y0*(1-A)+Y1*A
#'tau.adj(Y, A, B, X, pi, q)
#'@export
tau.adj<-function(Y, A, B, X = NULL, pi, q, conf.level = 0.95){
  Input_check(Y, A, B, X, pi, q)
  strt <- unique(B)
  strt_num <- length(strt)
  if(is.null(X)){
    l <- stats::lm(Y~A+factor(B))
    estimate <- l$coefficients[2]
    sq <- rep(q, strt_num)
    stderr <- FEE_var(A, B, Y, sq, strt, strt_num, pi)
  }
  else{
    p <- ncol(X)
    formtext <- "Y~A+factor(B)"
    for(i in 1:p){
      eval(parse(text = paste("x", i, "= X[,", i, "]", sep = "")))
      formtext <- paste(formtext, "+ x", i, sep = "")
    }
    l <- eval(parse(text = paste("stats::lm(", formtext,")", sep = "")))
    estimate <- l$coefficients[2]
    sq <- rep(q, strt_num)
    stderr <- AN_var(A, B, Y, X, sq, strt, strt_num, pi)
  }
  testmethod<-"Regression without Interaction"
  tstat <- estimate/stderr
  pval <- dplyr::if_else(stats::pnorm(tstat)<(1-stats::pnorm(tstat)),stats::pnorm(tstat),1-stats::pnorm(tstat))*2
  cint <- c(estimate + stderr*stats::qnorm((1-conf.level)/2),estimate - stderr*stats::qnorm((1-conf.level)/2))
  attr(cint,"conf.level") <- conf.level
  names(tstat) <- "t"
  names(estimate) <- "difference in treatment effect"
  rval<-list(statistic = tstat, p.value = pval, conf.int = cint,
             estimate = estimate, stderr = stderr,
             method = testmethod)
  class(rval) <- "htest"
  return(rval)
}

#'Regression with interaction
#'
#'Estimating and inferring the treatment effect based on regression with
#'interaction.
#'
#'Estimating and inferring the treatment effect based on regression with interaction. It
#'implements the methods as described in Sections 3.3 and 4.3, Ma et al. (2020).
#'
#'@param Y a numeric vector of observed outcomes. Its length should be the same
#'  as the number of subjects.
#'@param A a numeric vector of treatment assignments. Its length should be the
#'  same as the number of subjects.
#'@param B a numeric vector of stratum labels. Its length should be the same as
#'  the number of subjects.
#'@param X an (optional) numeric design matrix containing additional covariates
#'  used in the regression.
#'@param pi a numeric value for the target treatment proportion in each stratum.
#'@param conf.level confidence level of the interval. Default is 0.95.
#'
#'@return 
#'  A list of class \code{"htest"} containing the following components: 
#'  \item{statistic}{the value of the t-statistic.}
#'  \item{p.value}{the p-value for the test.}
#'  \item{conf.int}{a confidence interval under chosen level \code{conf.level} for the difference
#'  in treatment effect between treatment group and control group.}
#'  \item{estimate}{estimated treatment effect difference between treatment
#'  group and control group.}
#'  \item{method}{a character string indicating what type of regression was performed.} 
#'
#'@references Ma, W., Tu, F., & Liu, H. (2020). \emph{Regression analysis for
#'  covariate-adaptive randomization: A robust and efficient inference
#'  perspective}. arXiv preprint arXiv:2009.02287.
#'
#'@examples
#'#The code replicates the simulation setting of Model 2 in Section 6, Ma et al. (2020).
#'n <- 1000
#'pi <- 0.5
#'alphavec <- c(5,8,3,12)
#'m2e<-function(x){
#'  6*exp(x)*x*(1-x)
#'}
#'X1 <- rgamma(n,2)
#'X2 <- sample(c(1,2,3),n,replace = TRUE, prob = c(0.3,0.6,0.1))
#'X3 <- rpois(n,3)
#'X4 <- rbeta(n,2,2)
#'X1_S <- rep(1,n)
#'X1_S[which(X1 >= 2.5)] <- 2
#'B <- as.numeric(interaction(X1_S,X2))
#'X <- cbind(X1,X3)
#'A <- sample(c(0,1),n,replace=TRUE,prob=c(1-pi,pi))
#'Y0 <- alphavec[2]*3.6-alphavec[1]*2-alphavec[4]*integrate(m2e,lower = 0,upper = 1)$value+
#'      alphavec[1]*X1+log(alphavec[3]*X1*log(X3+1)+1)+alphavec[4]*exp(X4)+rnorm(n,sd = 2)
#'Y1 <- alphavec[2]*X2^2+log(alphavec[3]*X1*log(X3+1)+1)+rnorm(n,sd = 1)
#'Y <- Y0*(1-A)+Y1*A
#'tau.interact(Y, A, B, X, pi)
#'@export
tau.interact<-function(Y, A, B, X = NULL, pi, conf.level = 0.95){
  Input_check(Y, A, B, X, pi, 0)
  strt <- unique(B)
  strt_num <- length(strt)
  if(is.null(X)){
    dummy<-stats::model.matrix(~factor(B))
    dummy<-dummy[,-1]
    dummy_cent <- center_colmeans(dummy)
    l <- stats::lm(Y~A+dummy+A:dummy_cent)
    estimate <- l$coefficients[2]
    stderr <- SE_var(A, B, Y, strt, strt_num, pi)
  }
  else{
    p <- ncol(X)
    dummy<-stats::model.matrix(~factor(B))
    dummy<-dummy[,-1]
    dummy_cent <- center_colmeans(dummy)
    X_cent <- center_colmeans(X)
    formtext <- "Y~A+dummy+A:dummy_cent"
    for(i in 1:p){
      eval(parse(text = paste("x", i, "= X[,", i, "]", sep = "")))
      formtext <- paste(formtext, "+ x", i, sep = "")
      eval(parse(text = paste("x_cent", i, "= X_cent[,", i, "]", sep = "")))
      formtext <- paste(formtext, "+ A:x_cent", i, sep = "")
    }
    l <- eval(parse(text = paste("stats::lm(", formtext, ")",sep = "")))
    estimate <- l$coefficients[2]
    stderr <- IR_var(A, B, Y, X, strt, strt_num, pi)
  }
  testmethod<-"Regression with interaction"
  tstat <- estimate/stderr
  pval <- dplyr::if_else(stats::pnorm(tstat)<(1-stats::pnorm(tstat)),stats::pnorm(tstat),1-stats::pnorm(tstat))*2
  cint <- c(estimate + stderr*stats::qnorm((1-conf.level)/2),estimate - stderr*stats::qnorm((1-conf.level)/2))
  attr(cint,"conf.level") <- conf.level
  names(tstat) <- "t"
  names(estimate) <- "difference in treatment effect"
  rval<-list(statistic = tstat, p.value = pval, conf.int = cint,
             estimate = estimate, stderr = stderr,
             method = testmethod)
  class(rval) <- "htest"
  return(rval)
}