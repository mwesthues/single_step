# Efficient LOOCV With Prediction Errors
pacman::p_load("BGLR", "tidyverse")
data("mice")
X <- scale(mice.X)
h2 <- 0.5
QTL <- seq(from = 50, to = ncol(X), length = 20)
nQTL <- length(QTL)
n <- nrow(X)
b <- rep(1, nQTL) * sqrt(h2 / nQTL)
signal <- X[, QTL] %*% b
error <- rnorm(n, sd = sqrt(1 - h2))
y <- signal + error

G <- X %>%
  list() %>%
  map(~tcrossprod(scale(.)) / ncol(.)) %>%
  map(chol) %>%
  map(t) %>%
  .[[1]]

## Using a flat prior
system.time(fm1 <- BGLR(y = y, 
                        ETA = list(list(X = G, model = "BRR")), 
                        nIter = 30000, 
                        burnIn = 15000,
                        verbose = FALSE,
                        saveAt = "./tmp/"))
X_prime <- cbind(matrix(1, nrow = nrow(X), ncol = 1),
                 X)
b_hat <- fm1$ETA[[1]]$b
mu_hat <- fm1$mu
b_hat_prime <- matrix(c(mu_hat, b_hat), ncol = 1)
lambda <- c(fm1$varE / fm1$ETA[[1]]$varB, 1)
Diag <- cbind(diag(nrow = ncol(X), ncol = ncol(X)),
              matrix(1, nrow = ncol(X), ncol = 1))
Diag <- rbind(Diag,
              matrix(1, nrow = 1, ncol = ncol(Diag)))
H <- X_prime %*% solve(t(X_prime) %*% X_prime + Diag * lambda) %*% t(X_prime)
e_hat <- vapply(seq_len(nrow(X_prime)), FUN = function(j) {
t_x_prime <- X_prime[j, , drop = FALSE]
numerator <- as.numeric(y[j] - t_x_prime %*% b_hat_prime)
denominator <- 1 - H[j, j]
numerator / denominator
}, FUN.VALUE = numeric(1))
press <- sum(e_hat ^ 2)
rmse <- sqrt(press / length(e_hat))
y_hat <- y - e_hat
# Leave-one-out prediction accuracy
r <- cor(y, y_hat)




# Claas Heuer, October 2016
#
# Fast implementation of leave-one-out cross-validation
# without reestimating variance components in a ridge
# regression framework.
# Based on Rohan's paper
library(pacman)
p_load(rrBLUP, Matrix)
LOOCVFast <- function(y, M) {
  if(!is.vector(y) | !is.numeric(y)) stop("y has to be a numeric vector")
  if(anyNA(y)) stop("no NAs allowed in y")
  if(nrow(M) != length(y)) stop("dimensions for y and M dont match")
  if(anyNA(M)) stop("no NAs allowed in M")
  # check dimensions for most efficient scheme
  if(nrow(M) < ncol(M)) {
    G <- tcrossprod(M)
    L <- t(chol(G))
  } else {
    L <- M
  }
  mod <- mixed.solve(y = y, Z = L)
  bs <- c(mod$beta, mod$u)
  vA <- mod$Vu
  vE <- mod$Ve
  lambda <- vE / vA
  D <- c(0, rep(lambda, ncol(L)))
  X <- cbind(1, L)
  Htemp <- crossprod(X)
  diag(Htemp) <- diag(Htemp) + D
  H <- X %*% solve(Htemp) %*% t(X)
  E <- (y - X %*% bs) / (1 - diag(H))
  yhat <- y - as.numeric(E)
  return(list(y = y, yhat = yhat, e = as.numeric(E), cor = cor(y, yhat)))
}
out <- LOOCVFast(c(y), X)
