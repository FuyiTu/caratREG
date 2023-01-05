n <- 1000
pi <- 0.5
q <- pi*(1-pi)
X1 <- rgamma(n,2)
X2 <- sample(c(1,2,3),n,replace = TRUE, prob = c(0.3,0.6,0.1))
X3 <- rpois(n,3)
X4 <- rbeta(n,2,2)
X1_S <- ifelse(X1 < 2.5, 1, 2)
B <- as.numeric(interaction(X1_S,X2))
X <- cbind(X1,X3)
A <- sample(c(0,1),n,replace=TRUE,prob=c(1-pi,pi))
alphavec <- c(5,8,3,12)
Y0 <- alphavec[1]*X1+log(alphavec[3]*X1*log(X3+1)+1)+alphavec[4]*exp(X4)+rnorm(n,sd = 2)
Y1 <- 1.483708+alphavec[2]*X2^2+log(alphavec[3]*X1*log(X3+1)+1)+rnorm(n,sd = 1)
Y <- Y0*(1-A)+Y1*A
tau.diff(Y, A, B, X=NULL, pi, q)

test_that("NULL X works", {
  expect_true(is.numeric(tau.diff(Y, A, B, X=NULL, pi, q)$stat))
  expect_true(is.numeric(tau.adj(Y, A, B, X=NULL, pi, q)$stat))
  expect_true(is.numeric(tau.interact(Y, A, B, X=NULL, pi, q)$stat))
})
