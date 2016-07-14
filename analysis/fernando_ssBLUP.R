id1 <- c(3, 5, 6)
id2 <- c(1, 2, 4)
Ainv <- matrix(c(1.5, 0, -1, 0.5, 0, 0,
                 0, 2, 0, -1, -1, 0,
                 -1, 0, 2, -1, 0, 0, 
                 0.5, -1, -1, 2.5, 1, -1,
                 0, -1, 0, 1, 2, -1,
                 0, 0, 0, -1, -1, 2),
               nrow = 6, ncol = 6,
               dimnames = list(c(id1, id2), c(id1, id2)))
A_up11 <- Ainv[id1, id1]
A_up22 <- Ainv[id2, id2]
A_up12 <- Ainv[id1, id2]
M2 <- matrix(c(1, 2, 1,
               2, 1, 1,
               1, 1, 0, 
               1, 1, 1, 
               0, 2, 1,
               0, 0, 0,
               1, 1, 2,
               2, 1, 1, 
               1, 1, 2,
               0, 1, 1),
             nrow = 3, ncol = 10,
             dimnames = list(id2,
                             paste0("m", seq_len(10))))
M2 <- scale(M2, center = TRUE, scale = FALSE)

Pheno <- data.frame(Individual = seq_len(6),
                    y = c(0, 1.25, -0.34, 1.3, 1.27, 0.46),
                    stringsAsFactors = FALSE)

# Equation 21
M1 <- solve(A_up11, -A_up12 %*% M2)

# Equation 22
J2 <- matrix(rep(-1, times = ncol(A_up12)), nrow = ncol(A_up12), ncol = 1)
J1 <- solve(A_up11, -A_up12 %*% J2)


# --- Equation 20
# y
y1 <- Pheno[Pheno$Individual %in% id1, "y"]
y2 <- Pheno[Pheno$Individual %in% id2, "y"]
y <- c(y1, y2)
X1 <- matrix(1, nrow = length(id1), ncol = 1)
X1_prime <- cbind(X1, J1)
X2 <- matrix(1, nrow = length(id2), ncol = 1)
X2_prime <- cbind(X2, J2)
X_prime <- rbind(X1_prime, X2_prime)
Z1 <- diag(1, nrow = length(id1), ncol = length(id1))
Z2 <- diag(1, nrow = length(id2), ncol = length(id2))
W1 <- Z1 %*% M1
W2 <- Z2 %*% M2
W <- rbind(W1, W2)
U <- rbind(Z1, matrix(0, nrow = length(id2), ncol = ncol(Z1)))


# --- Equation 23
sigma2_g <- 0.2
sigma2_a <- sigma2_g / 10
sigma2_e <- 9 * sigma2_g
# Left hand side
lhs1 <- cbind(crossprod(X_prime), 
              crossprod(X_prime, W),
              crossprod(X1_prime, Z1))
lhs2 <- cbind(crossprod(W, X_prime),
              crossprod(W) + diag(1, nrow = ncol(W), 
                                  ncol = ncol(W)) * sigma2_e / sigma2_a,
              crossprod(W1, Z1))
lhs3 <- cbind(crossprod(Z1, X1_prime),
              crossprod(Z1, W1),
              crossprod(Z1) + A_up11 * sigma2_e / sigma2_g)
lhs <- rbind(lhs1, lhs2, lhs3)

# Right hand side
rhs <- rbind(crossprod(X_prime, y),
             crossprod(W, y),
             crossprod(Z1, y1))
c(solve(lhs, rhs))


# -- Equation 32
P <- diag(1, nrow = nrow(Z1), ncol = ncol(Z1)) - 
  Z1 %*% solve(crossprod(Z1) + A_up11 * sigma2_e / sigma2_g) %*% t(Z1)
lhs11 <- crossprod(X1_prime, P) %*% X1_prime + crossprod(X2_prime)
lhs12 <- crossprod(X1_prime, P) %*% W1 + crossprod(X2_prime, W2)
lhs21 <- crossprod(W1, P) %*% X1_prime + crossprod(W2, X2_prime)
lhs22 <- crossprod(W1, P) %*% W1 + crossprod(W2) + 
  diag(1, nrow = ncol(W1), ncol = ncol(W1)) * sigma2_e / sigma2_a

rhs1 <- crossprod(X1_prime, P) %*% y1 + crossprod(X2_prime, y2)
rhs2 <- crossprod(W1, P) %*% y1 + crossprod(W2, y2)
solve(rbind(cbind(lhs11, lhs12), cbind(lhs21, lhs22)),
      rbind(rhs1, rhs2))


# ---------------------------------------------------------------------------
pacman::p_load("cpgen")
id <- seq_len(6)
sire <- c(rep(NA, times = 3), rep(1, times = 3))
dam <- c(rep(NA, times = 3), 2, 2, 3)
y <- c(NA_real_, 1.25, -0.34, 1.3, 1.27, 0.46)
dat <- data.frame(id = id, sire = sire, dam = dam, y = y)
M <- matrix(c(1, 2, 1,
              2, 1, 1,
              1, 1, 0, 
              1, 1, 1, 
              0, 2, 1,
              0, 0, 0,
              1, 1, 2,
              2, 1, 1, 
              1, 1, 2,
              0, 1, 1),
            nrow = 3, ncol = 10)
M.id <- seq_len(3)
var_y <- var(y, na.rm = TRUE)
var_e <- (10 * var_y / 21)
var_a <- var_e
var_m <- var_e / 10
dfree <- 500
par_random <- list(list(method = "ridge", scale = var_m, df = dfree),
                   list(method = "ridge", scale = var_a, df = dfree))
set_num_threads(1)
mod <- cSSBR(data = dat,
             M = M,
             M.id = M.id,
             par_random = par_random,
             scale_e = var_e, 
             df_e = dfree,
             niter = 5e+04,
             burnin = 3e+04)
print(round(mod$Effect_1$posterior$estimates_mean, digits = 2))
print(round(mod$SSBR$Breeding_Values, digits = 2))
