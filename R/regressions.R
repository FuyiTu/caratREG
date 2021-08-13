# center at sample mean
center_colmeans <- function(x){
  xcenter <- colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

# Input check
Input_check <- function(A, B, Y, X, pi, q){
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
#'Estimating and inferring the treatment effect based on difference in means and
#'(optional) additional covatiates. It follows the idea of Ma et al. (2020)
#'<arXiv:2009.02287>, Section 3.1 and Section 4.1.
#'
#'@param A numeric vector, containing subjects' treatment assignments. Its
#'  length should be the same as the number of subjects.
#'@param B numeric vector, containing subjects' stratum labels. Its length
#'  should be the same as the number of subjects.
#'@param Y numeric vector, containing subjects' observed outcomes. Its length
#'  should be the same as the number of subjects.
#'@param X design matrix, containing additional covariates used in the
#'  regression (optional). Each column represents a covariate.
#'@param pi numeric, the target treatment proportion in each stratum.
#'@param q numeric, indicating the balance level of covariate-adaptive
#'  randomizations. Detailed information can be found in Ma et al.(2020),
#'  section 2.
#'@param conf.level confidence level of the interval. Default is 0.95.
#'
#'@return It returns an object of class \code{"htest"}.
#'
#'  The function \code{print} is used to obtain results. The generic accessor
#'  functions \code{statistic}, \code{p.value}, \code{conf.int} and others
#'  extract various useful features of the value returned by \code{dme.test}.
#'
#'  An object of class \code{"htest"} is a list containing at least the
#'  following components: \describe{ \item{statistic}{the value of the
#'  t-statistic.} \item{p.value}{the p-value of the test,the null hypothesis is
#'  rejected if p-value is less than the significance level.} \item{conf.int}{a
#'  confidence interval under chosen level \code{conf.level} for the difference in
#'  treatment effect between treatment group and control group.}
#'  \item{estimate}{estimated treatment effect difference between treatment
#'  group and control group.} \item{method}{a character string indicating what
#'  type of test was performed.} }
#'
#'@references Ma, W., Tu, F., & Liu, H. (2020). \emph{Regression analysis for
#'  covariate-adaptive randomization: A robust and efficient inference
#'  perspective}. arXiv preprint arXiv:2009.02287.
#'
#'@examples
#'n <- 1000
#'pi <- 0.5
#'q <- pi*(1-pi)
#'X1 <- rbeta(n, 2, 2)
#'X2 <- runif(n, -5, 3)
#'B <- sample(1:4, n, replace = TRUE)
#'A <- sample(c(0, 1), n, replace = TRUE, prob = c(1-pi, pi))
#'Y0 <- 3*log(X1+3)*X2 + rnorm(n,sd = 1)
#'Y1 <- 1 + 2*X1 + rnorm(n,sd = 2)
#'Y <- Y0*(1-A) + Y1*A
#'X <- cbind(X1, X2)
#'tau.diff(A, B, Y, X, pi, q)
#'@export
tau.diff<-function(A, B, Y, X = NULL, pi, q, conf.level = 0.95){
  Input_check(A, B, Y, X, pi, q)
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
#'interaction and (optional) additional covatiates. It follows the idea of Ma et
#'al. (2020) <arXiv:2009.02287>, Section 3.2 and Section 4.2.
#'
#'@param A numeric vector, containing subjects' treatment assignments. Its
#'  length should be the same as the number of subjects.
#'@param B numeric vector, containing subjects' stratum labels. Its length
#'  should be the same as the number of subjects.
#'@param Y numeric vector, containing subjects' observed outcomes. Its length
#'  should be the same as the number of subjects.
#'@param X design matrix, containing additional covariates used in the
#'  regression (optional). Each column represents a covariate.
#'@param pi numeric, the target treatment proportion in each stratum.
#'@param q numeric, indicating the balance level of covariate-adaptive
#'  randomizations. Detailed information can be found in Ma et al.(2020),
#'  section 2.
#'@param conf.level confidence level of the interval. Default is 0.95.
#'
#'@return It returns an object of class \code{"htest"}.
#'
#'  The function \code{print} is used to obtain results. The generic accessor
#'  functions \code{statistic}, \code{p.value}, \code{conf.int} and others
#'  extract various useful features of the value returned by \code{adj.test}.
#'
#'  An object of class \code{"htest"} is a list containing at least the
#'  following components: \describe{ \item{statistic}{the value of the
#'  t-statistic.} \item{p.value}{the p-value of the test,the null hypothesis is
#'  rejected if p-value is less than the significance level.} \item{conf.int}{a
#'  confidence interval under chosen level \code{conf.level} for the difference in
#'  treatment effect between treatment group and control group.}
#'  \item{estimate}{estimated treatment effect difference between treatment
#'  group and control group.} \item{method}{a character string indicating what
#'  type of test was performed.} }
#'
#'@references Ma, W., Tu, F., & Liu, H. (2020). \emph{Regression analysis for
#'  covariate-adaptive randomization: A robust and efficient inference
#'  perspective}. arXiv preprint arXiv:2009.02287.
#'
#'@examples
#'n <- 1000
#'pi <- 0.5
#'q <- pi*(1-pi)
#'X1 <- rbeta(n, 2, 2)
#'X2 <- runif(n, -5, 3)
#'B <- sample(1:4, n, replace = TRUE)
#'A <- sample(c(0, 1), n, replace = TRUE, prob = c(1-pi, pi))
#'Y0 <- 3*log(X1+3)*X2 + rnorm(n,sd = 1)
#'Y1 <- 1 + 2*X1 + rnorm(n,sd = 2)
#'Y <- Y0*(1-A) + Y1*A
#'X <- cbind(X1, X2)
#'tau.adj(A, B, Y, X, pi, q)
#'@export
tau.adj<-function(A, B, Y, X = NULL, pi, q, conf.level = 0.95){
  Input_check(A, B, Y, X, pi, q)
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
#'Estimating and inferring the treatment effect based on regression with
#'interaction and (optional) additional covatiates. It follows the idea of Ma et
#'al. (2020) <arXiv:2009.02287>, Section 3.3 and Section 4.3.
#'
#'@param A numeric vector, containing subjects' treatment assignments. Its
#'  length should be the same as the number of subjects.
#'@param B numeric vector, containing subjects' stratum labels. Its length
#'  should be the same as the number of subjects.
#'@param Y numeric vector, containing subjects' observed outcomes. Its length
#'  should be the same as the number of subjects.
#'@param X design matrix, containing additional covariates used in the
#'  regression (optional). Each column represents a covariate.
#'@param pi numeric, the target treatment proportion in each stratum.
#'@param q numeric, indicating the balance level of covariate-adaptive
#'  randomizations. Detailed information can be found in Ma et al.(2020),
#'  section 2.
#'@param conf.level confidence level of the interval. Default is 0.95.
#'
#'@return It returns an object of class \code{"htest"}.
#'
#'  The function \code{print} is used to obtain results. The generic accessor
#'  functions \code{statistic}, \code{p.value}, \code{conf.int} and others
#'  extract various useful features of the value returned by \code{inter.test}.
#'
#'  An object of class \code{"htest"} is a list containing at least the
#'  following components: \describe{ \item{statistic}{the value of the
#'  t-statistic.} \item{p.value}{the p-value of the test,the null hypothesis is
#'  rejected if p-value is less than the significance level.} \item{conf.int}{a
#'  confidence interval under chosen level \code{conf} for the difference in
#'  treatment effect between treatment group and control group.}
#'  \item{estimate}{estimated treatment effect difference between treatment
#'  group and control group.} \item{method}{a character string indicating what
#'  type of test was performed.} }
#'
#'@references Ma, W., Tu, F., & Liu, H. (2020). \emph{Regression analysis for
#'  covariate-adaptive randomization: A robust and efficient inference
#'  perspective}. arXiv preprint arXiv:2009.02287.
#'
#'@examples
#'n <- 1000
#'pi <- 0.5
#'q <- pi*(1-pi)
#'X1 <- rbeta(n, 2, 2)
#'X2 <- runif(n, -5, 3)
#'B <- sample(1:4, n, replace = TRUE)
#'A <- sample(c(0, 1), n, replace = TRUE, prob = c(1-pi, pi))
#'Y0 <- 3*log(X1+3)*X2 + rnorm(n,sd = 1)
#'Y1 <- 1 + 2*X1 + rnorm(n,sd = 2)
#'Y <- Y0*(1-A) + Y1*A
#'X <- cbind(X1, X2)
#'tau.interact(A, B, Y, X, pi, q)
#'@export
tau.interact<-function(A, B, Y, X = NULL, pi, q, conf.level = 0.95){
  Input_check(A, B, Y, X, pi, q)
  strt <- unique(B)
  strt_num <- length(strt)
  if(is.null(X)){
    dummy<-stats::model.matrix(~factor(B))
    dummy<-dummy[,-1]
    dummy_cent <- center_colmeans(dummy)
    l <- stats::lm(Y~A+dummy+A:dummy_cent)
    estimate <- l$coefficients[2]
    sq <- rep(q, strt_num)
    stderr <- SE_var(A, B, Y, sq, strt, strt_num, pi)
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
    sq <- rep(q, strt_num)
    stderr <- IR_var(A, B, Y, X, sq, strt, strt_num, pi)
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