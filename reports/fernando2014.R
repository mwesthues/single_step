if (!require("pacman")) install.packages("pacman")
pacman::p_load("clusterGeneration", "ggplot2", "MASS")

# Number of observations
n_obs <- 50
# Number of observations without genotypes.
n_1 <- seq(from = 1, to = n_obs / 2, by = 1)
# Number of observations with genotypes.
n_2 <- seq(from = max(n_1) + 1, to = max(n_obs), by = 1)

# Generate a positive definite numerator relationship matrix.
A <- genPositiveDefMat(dim = n_obs, covMethod = "eigen",
                       rangeVar = c(0, 2))$Sigma
# Generate a positive definite genomic relationship matrix.
G <- genPositiveDefMat(dim = n_obs, covMethod = "eigen", 
                       rangeVar = c(0, 2))$Sigma

# Block matrices
## Numerator relationship matrix
A11 <- A[n_1, n_1]
A12 <- A[n_1, n_2]
A21 <- A[n_2, n_1]
A22 <- A[n_2, n_2]
## Genomic relationship matrix
G11 <- G[n_1, n_1]
G12 <- G[n_1, n_2]
G21 <- G[n_2, n_1]
G22 <- G[n_2, n_2]


# Legarra et al. (2009), Equation (4) 
H11 <- A11 + A12 %*% solve(A22) %*% (G22 - A22) %*% solve(A22) %*% A21
H12 <- A12 %*% solve(A22) %*% solve(G22)
H21 <- G22 %*% solve(A22) %*% A21
H22 <- G22
H <- cbind(rbind(H11, H21), rbind(H12, H22))
Hinv <- solve(H)


# Legarra et al. (2009), Equation (6) 
Diag <- diag(1, nrow = n_obs / 2, ncol = n_obs / 2)
Zero <- matrix(0, nrow = n_obs / 2, ncol = n_obs / 2)
block1 <- cbind(rbind(A12 %*% solve(A22),
                      Zero),
                rbind(Zero,
                      Diag))
block2 <- rbind(Diag, Diag)
block3 <- G22 - A22
block4 <- cbind(Diag, Diag)
block5 <- cbind(rbind(solve(A22) %*% A21,
                      Zero),
                rbind(Zero, Diag))
H_a <- A + block1 %*% block2 %*% block3 %*% block4 %*% block5
Hinv_a <- solve(H_a)


# Christensen and Lund (2010)
Hinv_b <- solve(A) + cbind(rbind(Zero, Zero),
                           rbind(Zero, solve(G22) - solve(A22)))


# Vitezica et al. (2011)
alpha_ <- 1 / n_obs ^ 2 * (sum(A22) - sum(G22))
H11_c <- A11
H12_c <- A12
H21_c <- A21
H22_c <- A22 + solve(G22 + 
                     matrix(1, nrow = length(n_2), length(n_2)) * alpha_) -
         solve(A22)
Hinv_c <- cbind(rbind(H11_c, H21_c), rbind(H12_c, H22_c))



# Compare different methods
tri_hinv <- Hinv[upper.tri(Hinv, diag = FALSE)]
tri_hinv_a <- Hinv_a[upper.tri(Hinv_a, diag = FALSE)]
tri_hinv_b <- Hinv_b[upper.tri(Hinv_b, diag = FALSE)]
tri_hinv_c <- Hinv_c[upper.tri(Hinv_c, diag = FALSE)]
hinv_df <- data.frame(Values = c(tri_hinv, tri_hinv_a, tri_hinv_b, tri_hinv_c),
                      Type = c(rep("Legarra_Eq4", times = length(tri_hinv)),
                               rep("Legarra_Eq6", times = length(tri_hinv_a)),
                               rep("Christ_Eq8", times = length(tri_hinv_b)),
                               rep("Vitezica_Eq2", times = length(tri_hinv_c))))
hinv_df$Type <- as.factor(hinv_df$Type)
ggplot(hinv_df, aes(x = Values, fill = Type)) +
  geom_density(alpha = 0.4) +
  theme(legend.position = "top")


### --- Simulate a response, with random variables that are correlated to it as
### specified by a genomic relationship matrix "G" defined above.
# http://www.quantumforest.com/2011/10/simulating-data-following-a-given-
# covariance-structure/
M <- genPositiveDefMat(dim = n_obs + 1, covMethod = "eigen",
                       rangeVar = c(0, 2))$Sigma
# Cholesky decomposition             
L <- chol(M)
# R chol function produces an upper triangular version of L
# so we have to transpose it.
# Just to be sure we can have a look at t(L) and the
# product of the Cholesky decomposition by itself
# t(L)
# t(L) %*% L

# Random variables that follow an M correlation matrix
r <- t(L) %*% matrix(rnorm((n_obs + 1) * n_obs), nrow = n_obs + 1,
                     ncol = n_obs)
r <- t(r)
rdata <- as.data.frame(r)
names(rdata) <- c('resp', paste0("pred", seq_len(n_obs)))

