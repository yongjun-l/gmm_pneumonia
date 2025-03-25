# Kline and Michael's code
pneumonia_data <- function(
    n = 1000,
    k = 3,
    p_L = 0.05, # The proportion of n that truly has pneumonia
    p_Z = c(-2.0, 0.0, 2.0), # Intercepts representing the baseline effect if L = 0
    slope_Z = c(-1.0, 0.5, 0.05), # Slopes that weight the influence of L on Z
    p_A = 0.01, # Baseline ignoring influence of Z on A
    p_AZ = 0.5, # Additional Effect of Z = z on A
    p_L_star = -4.0, # Baseline L_star = 1 if A and L are 0
    p_LL_star = 5.0, # Influence of L on L_star
    p_AL_star = 0.01, # Influence of A on L_star
    p_Y = -3, # Baseline Y = 1 if L_star, A, and L are 0
    p_YL_star = 0.01, # Influence of L_star on Y
    p_YA = -2, # Influence of A on Y
    p_YL = 0.05 # Influence of L on Y
){

  # Defunct
  # # Simulate a latent factor U to induce correlation between L and Z
  # # This is to make sure L and Z are correlated since we don't have a directional
  # # arrow for that part of the dag.
  # U <- rnorm(n)

  # Get L for all the subjects
  n_L <- ceiling(n*p_L)
  prob_L <- c(rep(0, times = n-n_L), rep(1, times = n_L))
  L <- sample(prob_L, n)

  # Get Z for all subjects
  linear_preds <- sapply(1:k, function(x) p_Z[x] + slope_Z[x] * L)
  exp_preds <- exp(linear_preds)
  prob_Z <- t(apply(exp_preds, 1, function(x) x / sum(x))) # Individual probability that subject i is assigned to Z = z
  Z <- apply(prob_Z, 1, function(probs) sample(1:k, 1, prob = probs))

  # Get A conditional on Z
  prob_A <- plogis(p_A + p_AZ * Z)
  A <- rbinom(n, 1, prob_A)

  # Get L_star
  prob_L_star <- plogis(p_L_star + p_AL_star * A + p_LL_star * L)
  L_star <- rbinom(n, 1, prob_L_star)

  # Get Y
  prob_Y <- plogis(p_Y + p_YL_star * L_star + p_YA * A + p_YL * L)
  Y <- rbinom(n, 1, prob_Y)

  # Return the data
  dat <- data.frame(L, Z, A, L_star, Y)

  param <- list(
    # P(L)
    p_L = p_L,
    # P(Z|L)
    p_Z1_L0 = c(exp(p_Z) / sum(exp(p_Z)))[1],
    p_Z2_L0 = c(exp(p_Z) / sum(exp(p_Z)))[2],
    p_Z1_L1 = c(exp(p_Z + slope_Z) / sum(exp(p_Z + slope_Z)))[1],
    p_Z2_L1 = c(exp(p_Z + slope_Z) / sum(exp(p_Z + slope_Z)))[2],
    # P(A|Z)
    p_A_Z1 = plogis(p_A + p_AZ * 1),
    p_A_Z2 = plogis(p_A + p_AZ * 2),
    p_A_Z3 = plogis(p_A + p_AZ * 3),
    # P(L*|A, L)
    p_Ls_A0L0 = plogis(p_L_star + p_AL_star * 0 + p_LL_star * 0),
    p_Ls_A0L1 = plogis(p_L_star + p_AL_star * 0 + p_LL_star * 1),
    p_Ls_A1L0 = plogis(p_L_star + p_AL_star * 1 + p_LL_star * 0),
    p_Ls_A1L1 = plogis(p_L_star + p_AL_star * 1 + p_LL_star * 1),
    # P(Y|L*, A, L)
    p_Y_Ls0A0L0 = plogis(p_Y + p_YL_star * 0 + p_YA * 0 + p_YL * 0),
    p_Y_Ls0A0L1 = plogis(p_Y + p_YL_star * 0 + p_YA * 0 + p_YL * 1),
    p_Y_Ls0A1L0 = plogis(p_Y + p_YL_star * 0 + p_YA * 1 + p_YL * 0),
    p_Y_Ls0A1L1 = plogis(p_Y + p_YL_star * 0 + p_YA * 1 + p_YL * 1),
    p_Y_Ls1A0L0 = plogis(p_Y + p_YL_star * 1 + p_YA * 0 + p_YL * 0),
    p_Y_Ls1A0L1 = plogis(p_Y + p_YL_star * 1 + p_YA * 0 + p_YL * 1),
    p_Y_Ls1A1L0 = plogis(p_Y + p_YL_star * 1 + p_YA * 1 + p_YL * 0),
    p_Y_Ls1A1L1 = plogis(p_Y + p_YL_star * 1 + p_YA * 1 + p_YL * 1)
  )

  return(list(dat=dat, param=param))
}

# Yongjun's simulation version
pneumonia_data_yj <- function(
    n=1000,
    k=3
){

  # Defunct
  # # Simulate a latent factor U to induce correlation between L and Z
  # # This is to make sure L and Z are correlated since we don't have a directional
  # # arrow for that part of the dag.
  # U <- rnorm(n)

  param <- list(
    # P(L)
    p_L = 0.8,
    # P(Z|L)
    p_Z1_L0 = 0.2,
    p_Z2_L0 = 0.3,
    p_Z1_L1 = 0.1,
    p_Z2_L1 = 0.5,
    # P(A|Z)
    p_A_Z1 = 0.45,
    p_A_Z2 = 0.55,
    p_A_Z3 = 0.65,
    # P(L*|A, L)
    p_Ls_A0L0 = 0.1,
    p_Ls_A0L1 = 0.7,
    p_Ls_A1L0 = 0.05,
    p_Ls_A1L1 = 0.8,
    # P(Y|L*, A, L)
    p_Y_Ls0A0L0 = 0.05,
    p_Y_Ls0A0L1 = 0.15,
    p_Y_Ls0A1L0 = 0.01,
    p_Y_Ls0A1L1 = 0.1,
    p_Y_Ls1A0L0 = 0.03,
    p_Y_Ls1A0L1 = 0.12,
    p_Y_Ls1A1L0 = 0.05,
    p_Y_Ls1A1L1 = 0.15
  )
  
  param <- list(
    # P(L)
    p_L = 0.2,
    # P(Z|L)
    p_Z1_L0 = 0.2,
    p_Z2_L0 = 0.2,
    p_Z1_L1 = 0.2,
    p_Z2_L1 = 0.2,
    # P(A|Z)
    p_A_Z1 = 0.2,
    p_A_Z2 = 0.2,
    p_A_Z3 = 0.2,
    # P(L*|A, L)
    p_Ls_A0L0 = 0.2,
    p_Ls_A0L1 = 0.2,
    p_Ls_A1L0 = 0.2,
    p_Ls_A1L1 = 0.2,
    # P(Y|L*, A, L)
    p_Y_Ls0A0L0 = 0.2,
    p_Y_Ls0A0L1 = 0.2,
    p_Y_Ls0A1L0 = 0.2,
    p_Y_Ls0A1L1 = 0.2,
    p_Y_Ls1A0L0 = 0.2,
    p_Y_Ls1A0L1 = 0.2,
    p_Y_Ls1A1L0 = 0.2,
    p_Y_Ls1A1L1 = 0.2
  )

  # Get L for all the subjects
  L <- rbinom(n, 1, param$p_L)

  # Get Z for all subjects
  prob_Z <- matrix(0, nrow = n, ncol = k)

  # Corrected Code
  prob_Z[L == 0, ] <- matrix(c(param$p_Z1_L0, param$p_Z2_L0, 1-param$p_Z1_L0-param$p_Z2_L0), nrow = sum(L == 0), ncol = k, byrow = TRUE)
  prob_Z[L == 1, ] <- matrix(c(param$p_Z1_L1, param$p_Z2_L1, 1-param$p_Z1_L1-param$p_Z2_L1), nrow = sum(L == 1), ncol = k, byrow = TRUE)
  Z <- apply(prob_Z, 1, function(probs) sample(1:k, 1, prob = probs))

  # Get A conditional on Z
  prob_A <- ifelse(Z == 1, param$p_A_Z1, ifelse(Z == 2, param$p_A_Z2, param$p_A_Z3))
  A <- rbinom(n, 1, prob_A)

  # Get L_star
  prob_L_star <- ifelse(A == 0 & L == 0, param$p_Ls_A0L0,
                        ifelse(A == 0 & L == 1, param$p_Ls_A0L1,
                               ifelse(A == 1 & L == 0, param$p_Ls_A1L0, param$p_Ls_A1L1)))
  L_star <- rbinom(n, 1, prob_L_star)

  # Get Y
  prob_Y <- ifelse(L_star == 0 & A == 0 & L == 0, param$p_Y_Ls0A0L0,
                   ifelse(L_star == 0 & A == 0 & L == 1, param$p_Y_Ls0A0L1,
                          ifelse(L_star == 0 & A == 1 & L == 0, param$p_Y_Ls0A1L0,
                                 ifelse(L_star == 0 & A == 1 & L == 1, param$p_Y_Ls0A1L1,
                                        ifelse(L_star == 1 & A == 0 & L == 0, param$p_Y_Ls1A0L0,
                                               ifelse(L_star == 1 & A == 0 & L == 1, param$p_Y_Ls1A0L1,
                                                      ifelse(L_star == 1 & A == 1 & L == 0, param$p_Y_Ls1A1L0, param$p_Y_Ls1A1L1)))))))
  Y <- rbinom(n, 1, prob_Y)

  # Return the data
  dat <- data.frame(L, Z, A, L_star, Y)

  return(list(dat=dat, param=param))
}


# Implement GMM
gmm_functions <- function(theta, dat) {

  p_L <- theta[1]

  p_Z1_L0 <- theta[2]
  p_Z2_L0 <- theta[3]
  p_Z1_L1 <- theta[4]
  p_Z2_L1 <- theta[5]

  p_A_Z1 <- theta[6]
  p_A_Z2 <- theta[7]
  p_A_Z3 <- theta[8]

  p_Ls_A0L0 <- theta[9]
  p_Ls_A0L1 <- theta[10]
  p_Ls_A1L0 <- theta[11]
  p_Ls_A1L1 <- theta[12]

  p_Y_Ls0A0L0 <- theta[13]
  p_Y_Ls0A0L1 <- theta[14]
  p_Y_Ls0A1L0 <- theta[15]
  p_Y_Ls0A1L1 <- theta[16]
  p_Y_Ls1A0L0 <- theta[17]
  p_Y_Ls1A0L1 <- theta[18]
  p_Y_Ls1A1L0 <- theta[19]
  p_Y_Ls1A1L1 <- theta[20]


  counts <- table(dat$Z, dat$A, dat$L_star, dat$Y) / nrow(dat)


  # Estimating equations
  # Z A Ls Y
  m1  <- counts[1,1,1,1] - ((1-p_L) * p_Z1_L0             * (1-p_A_Z1) * (1-p_Ls_A0L0) * (1-p_Y_Ls0A0L0) + p_L * p_Z1_L1             * (1-p_A_Z1) * (1-p_Ls_A0L1) * (1-p_Y_Ls0A0L1))
  m2  <- counts[1,1,1,2] - ((1-p_L) * p_Z1_L0             * (1-p_A_Z1) * (1-p_Ls_A0L0) * p_Y_Ls0A0L0     + p_L * p_Z1_L1             * (1-p_A_Z1) * (1-p_Ls_A0L1) * p_Y_Ls0A0L1)
  m3  <- counts[1,1,2,1] - ((1-p_L) * p_Z1_L0             * (1-p_A_Z1) * p_Ls_A0L0     * (1-p_Y_Ls1A0L0) + p_L * p_Z1_L1             * (1-p_A_Z1) * p_Ls_A0L1     * (1-p_Y_Ls1A0L1))
  m4  <- counts[1,1,2,2] - ((1-p_L) * p_Z1_L0             * (1-p_A_Z1) * p_Ls_A0L0     * p_Y_Ls1A0L0     + p_L * p_Z1_L1             * (1-p_A_Z1) * p_Ls_A0L1     * p_Y_Ls1A0L1)
  m5  <- counts[1,2,1,1] - ((1-p_L) * p_Z1_L0             * p_A_Z1     * (1-p_Ls_A1L0) * (1-p_Y_Ls0A1L0) + p_L * p_Z1_L1             * p_A_Z1     * (1-p_Ls_A1L1) * (1-p_Y_Ls0A1L1))
  m6  <- counts[1,2,1,2] - ((1-p_L) * p_Z1_L0             * p_A_Z1     * (1-p_Ls_A1L0) * p_Y_Ls0A1L0     + p_L * p_Z1_L1             * p_A_Z1     * (1-p_Ls_A1L1) * p_Y_Ls0A1L1)
  m7  <- counts[1,2,2,1] - ((1-p_L) * p_Z1_L0             * p_A_Z1     * p_Ls_A1L0     * (1-p_Y_Ls1A1L0) + p_L * p_Z1_L1             * p_A_Z1     * p_Ls_A1L1     * (1-p_Y_Ls1A1L1))
  m8  <- counts[1,2,2,2] - ((1-p_L) * p_Z1_L0             * p_A_Z1     * p_Ls_A1L0     * p_Y_Ls1A1L0     + p_L * p_Z1_L1             * p_A_Z1     * p_Ls_A1L1     * p_Y_Ls1A1L1)
  m9  <- counts[2,1,1,1] - ((1-p_L) * p_Z2_L0             * (1-p_A_Z2) * (1-p_Ls_A0L0) * (1-p_Y_Ls0A0L0) + p_L * p_Z2_L1             * (1-p_A_Z2) * (1-p_Ls_A0L1) * (1-p_Y_Ls0A0L1))
  m10 <- counts[2,1,1,2] - ((1-p_L) * p_Z2_L0             * (1-p_A_Z2) * (1-p_Ls_A0L0) * p_Y_Ls0A0L0     + p_L * p_Z2_L1             * (1-p_A_Z2) * (1-p_Ls_A0L1) * p_Y_Ls0A0L1)
  m11 <- counts[2,1,2,1] - ((1-p_L) * p_Z2_L0             * (1-p_A_Z2) * p_Ls_A0L0     * (1-p_Y_Ls1A0L0) + p_L * p_Z2_L1             * (1-p_A_Z2) * p_Ls_A0L1     * (1-p_Y_Ls1A0L1))
  m12 <- counts[2,1,2,1] - ((1-p_L) * p_Z2_L0             * (1-p_A_Z2) * p_Ls_A0L0     * p_Y_Ls1A0L0     + p_L * p_Z2_L1             * (1-p_A_Z2) * p_Ls_A0L1     * p_Y_Ls1A0L1)
  m13 <- counts[2,2,1,1] - ((1-p_L) * p_Z2_L0             * p_A_Z2     * (1-p_Ls_A1L0) * (1-p_Y_Ls0A1L0) + p_L * p_Z2_L1             * p_A_Z2     * (1-p_Ls_A1L1) * (1-p_Y_Ls0A1L1))
  m14 <- counts[2,2,1,2] - ((1-p_L) * p_Z2_L0             * p_A_Z2     * (1-p_Ls_A1L0) * p_Y_Ls0A1L0     + p_L * p_Z2_L1             * p_A_Z2     * (1-p_Ls_A1L1) * p_Y_Ls0A1L1)
  m15 <- counts[2,2,2,1] - ((1-p_L) * p_Z2_L0             * p_A_Z2     * p_Ls_A1L0     * (1-p_Y_Ls1A1L0) + p_L * p_Z2_L1             * p_A_Z2     * p_Ls_A1L1     * (1-p_Y_Ls1A1L1))
  m16 <- counts[2,2,2,2] - ((1-p_L) * p_Z2_L0             * p_A_Z2     * p_Ls_A1L0     * p_Y_Ls1A1L0     + p_L * (1-p_Z1_L1-p_Z2_L1) * p_A_Z2     * p_Ls_A1L1     * p_Y_Ls1A1L1)
  m17 <- counts[3,1,1,1] - ((1-p_L) * (1-p_Z1_L0-p_Z2_L0) * (1-p_A_Z3) * (1-p_Ls_A0L0) * (1-p_Y_Ls0A0L0) + p_L * (1-p_Z1_L1-p_Z2_L1) * (1-p_A_Z3) * (1-p_Ls_A0L1) * (1-p_Y_Ls0A0L1))
  m18 <- counts[3,1,1,2] - ((1-p_L) * (1-p_Z1_L0-p_Z2_L0) * (1-p_A_Z3) * (1-p_Ls_A0L0) * p_Y_Ls0A0L0     + p_L * (1-p_Z1_L1-p_Z2_L1) * (1-p_A_Z3) * (1-p_Ls_A0L1) * p_Y_Ls0A0L1)
  m19 <- counts[3,1,2,1] - ((1-p_L) * (1-p_Z1_L0-p_Z2_L0) * (1-p_A_Z3) * p_Ls_A0L0     * (1-p_Y_Ls1A0L0) + p_L * (1-p_Z1_L1-p_Z2_L1) * (1-p_A_Z3) * p_Ls_A0L1     * (1-p_Y_Ls1A0L1))
  m20 <- counts[3,1,2,2] - ((1-p_L) * (1-p_Z1_L0-p_Z2_L0) * (1-p_A_Z3) * p_Ls_A0L0     * p_Y_Ls1A0L0     + p_L * (1-p_Z1_L1-p_Z2_L1) * (1-p_A_Z3) * p_Ls_A0L1     * p_Y_Ls1A0L1)
  m21 <- counts[3,2,1,1] - ((1-p_L) * (1-p_Z1_L0-p_Z2_L0) * p_A_Z3     * (1-p_Ls_A1L0) * (1-p_Y_Ls0A1L0) + p_L * (1-p_Z1_L1-p_Z2_L1) * p_A_Z3     * (1-p_Ls_A1L1) * (1-p_Y_Ls0A1L1))
  m22 <- counts[3,2,1,2] - ((1-p_L) * (1-p_Z1_L0-p_Z2_L0) * p_A_Z3     * (1-p_Ls_A1L0) * p_Y_Ls0A1L0     + p_L * (1-p_Z1_L1-p_Z2_L1) * p_A_Z3     * (1-p_Ls_A1L1) * p_Y_Ls0A1L1)
  m23 <- counts[3,2,2,1] - ((1-p_L) * (1-p_Z1_L0-p_Z2_L0) * p_A_Z3     * p_Ls_A1L0     * (1-p_Y_Ls1A1L0) + p_L * (1-p_Z1_L1-p_Z2_L1) * p_A_Z3     * p_Ls_A1L1     * (1-p_Y_Ls1A1L1))

  # Quadratic form
  vec <- c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23)

  return(sum(vec^2))
}


# debugonce(pneumonia_data)
set.seed(362025)
test_dat <- pneumonia_data_yj(n=500000)
param <- unlist(test_dat$param)
dat <- test_dat$dat

# Sanity check
table(dat$L) / nrow(dat)
table(dat$L_star) / nrow(dat)
table(dat$L, dat$Z) / rowSums(table(dat$L, dat$Z))
table(dat$Z, dat$A) / rowSums(table(dat$Z, dat$A))

# Step 3: Optimize using optim()
theta_start <- rep(0.1,20)  # Initial guess
gmm_solution <- optim(theta_start, gmm_functions, dat = dat, method = "L-BFGS-B", lower = rep(0,20), upper = rep(1,20))
cbind(param,gmm_solution$par, abs(param-gmm_solution$par))

# Simulations
m=1
n=500000
probs_iter <- matrix(NA, nrow = m, ncol = 20)
for (i in 1:m) {
  print(i)
  test_dat <- pneumonia_data_yj(n=n)
  param <- unlist(test_dat$param)
  dat <- test_dat$dat

  # Step 3: Optimize using optim()
  theta_start <- rep(0.1,20)  # Initial guess
  gmm_solution <- optim(theta_start, gmm_functions, dat = dat, method = "L-BFGS-B", lower = rep(0,20), upper = rep(1,20))
  probs_iter[i,] <- gmm_solution$par
}

# BB package



rslts <- cbind(param,colMeans(probs_iter),apply(probs_iter, 2, sd),abs(param-colMeans(probs_iter)))
colnames(rslts) <- c("True", "Estimate", "SE", "Bias")
rslts <- round(rslts, 4)



library(flextable)
library(officer)

# Convert matrix to data frame and preserve row names
df <- as.data.frame(rslts)
df <- cbind(RowName = rownames(df), df)

# Create a Word document with a table
doc <- read_docx() |>
  body_add_flextable(flextable(df))

# Save the Word file
print(doc, target = "m500_n1000_yj_param2.docx")






































