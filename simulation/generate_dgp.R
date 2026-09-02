Xfun <- function(n, d, indep_covariate, rho = 0) {
if (indep_covariate == TRUE) {
    X <- matrix(
      rnorm(n * d), #changing to unif gives moderate overlap setting
      nrow = n,
      ncol = d
    )
  } else if (indep_covariate == FALSE) {
      X <- matrix(
      rnorm(n * d),
      nrow = n,
      ncol = d
    )
    fac <- rnorm(n)
    X <- X * sqrt(1 - rho) +
    fac * sqrt(rho)
    X <- pnorm(X)
  } else {
    stop("Neither independent nor correlated covariates chosen.")
  }
  return(X)
}

psfun <- function(X) {
  (1 + pbeta(X[, 1], 2, 4)) / 4
}
psfun_overlap <- function(X) {
  1 / (1 + exp(-X))
}


# S1: Complex linear model
gen_data_s1 <- function(
    sim = 10,
    n = 1000,
    ntest = 1000,
    d = 20,
    indep_covariate = TRUE,
    rho = 0
) {
  # Simulation IDs
  sim_id <- rep(1:sim, each = n)
   X <- Xfun(
    n * sim,
    d,
    indep_covariate,
    rho
  )
  colnames(X) <- paste0("X", 1:d)
  beta_1 <- runif(d, min = 0, max = 10)
  beta_0 <- runif(d, min = 0, max = 10)
  mu_0 <- X %*% beta_0
  mu_1 <- X %*% beta_1
  CATE <- mu_1 - mu_0
  Y1 <- mu_1 + rnorm(n * sim)
  Y0 <- mu_0 + rnorm(n * sim)
  
  ps <- psfun(X)
  #  pstest <-  psfun_overlap(1+4*Xtest[, 1])

  T_ind <- as.numeric(
    runif(n * sim) < ps
  )
  
  Yobs <- Y0
  Yobs[T_ind == 1] <- Y1[T_ind == 1]
  Xtest <- Xfun(
    ntest * sim,
    d,
    indep_covariate,
    rho
  )
  
  colnames(Xtest) <- paste0("Xtest", 1:d)
  
  mu_0_test <- Xtest %*% beta_0
  mu_1_test <- Xtest %*% beta_1
  
  CATE_test <- mu_1_test - mu_0_test
  
  Y1test <- mu_1_test + rnorm(ntest * sim)
  Y0test <- mu_0_test + rnorm(ntest * sim)
  
  pstest <- psfun(Xtest)
 #pstest <-  psfun_overlap(1+4*Xtest[, 1]) for overlap violation
  
  T_ind_test <- as.numeric(
    runif(ntest * sim) < pstest
  )
  # combine data

  dat <- tibble(
    sim_id,
    Y1 = Y1[, 1],
    Y0 = Y0[, 1],
    CATE = CATE[, 1],
    ps,
    T_ind,
    Yobs = Yobs[, 1],
    as_tibble(X),
    Y1test = Y1test[, 1],
    Y0test = Y0test[, 1],
    CATE_test = CATE_test[, 1],
    pstest,
    T_ind_test,
    as_tibble(Xtest)
  )
  
  return(
    dat %>%
      arrange(sim_id)
  )
}

# S2: Nonlinear CATE 

gen_data_s2 <- function(
    sim = 1,
    n = 1000,
    ntest = 1000,
    d = 10,
    indep_covariate,
    homosc_err,
    rho
) {
  
  sim_id <- rep(1:sim, each = n)
  
  errdist <- rnorm
  
  X <- Xfun(n * sim, d, indep_covariate, rho)
  colnames(X) <- paste0("X", seq(d))
  
  Y0 <- rep(0, n * sim)
  
  tau <- taufun(X)
  std <- sdfun(X, homosc_err)
  
  Y1 <- tau + std * errdist(n * sim)
  
  ps <- psfun(X)
  
  T_ind <- as.numeric(
    runif(n * sim) < ps
  )
  
  Yobs <- Y0
  Yobs[T_ind == 1] <- Y1[T_ind == 1]
  
  # Test data
  Xtest <- Xfun(ntest * sim, d, indep_covariate, rho)
  colnames(Xtest) <- paste0("Xtest", seq(d))
  
  Y0test <- rep(0, ntest * sim)
  
  tautest <- taufun(Xtest)
  stdtest <- sdfun(Xtest, homosc_err)
  
  Y1test <- tautest + stdtest * errdist(ntest * sim)
  
  pstest <- psfun(Xtest)
  
  T_ind_test <- as.numeric(
    runif(ntest * sim) < pstest
  )
  
  dat <- tibble(
    sim_id,
    Y1,
    Y0,
    tau,
    ps,
    T_ind,
    Yobs,
    as_tibble(X),
    Y1test,
    Y0test,
    tautest,
    pstest,
    T_ind_test,
    as_tibble(Xtest)
  )
  
  dat <- arrange(dat, sim_id)
  
  return(dat)
}
# S3: highly nonlinear with interaction terms 
gen_data_s3 <- function(
    sim = 1000,
    n = 1000,
    ntest = 1000,
    d = 10,
    indep_covariate = TRUE,
    rho = 0
) {
  
  sim_id <- rep(1:sim, each = n)
 beta_0 <- function(X) {
    100 +
      4 * X[, 1] +
      X[, 2] -
      3 * X[, 3]
  }
  
  CATE_fun <- function(X) {
    6 * sin(2 * X[, 1]) +
      3 * (X[, 2] + 3) * X[, 3] +
      9 * tanh(0.5 * X[, 4]) +
      3 * X[, 5] * (2 * as.numeric(X[, 4] > 0) - 1) +
      3 * X[, 6] +
      2 * X[, 7] +
      X[, 8] -
      2 * X[, 9] -
      4 * X[, 10]
  }
 X <- Xfun(
    n * sim,
    d,
    indep_covariate,
    rho
  )
  colnames(X) <- paste0("X", 1:d)
  ps <- psfun_overlap(
    1 + 4 * X[, 1]
  )
  tau <- CATE_fun(X)
  
  mu0 <- beta_0(X) - 0.5 * tau
  mu1 <- beta_0(X) + 0.5 * tau
  
  Y0 <- mu0 + rnorm(n * sim)
  Y1 <- mu1 + rnorm(n * sim)
  T_ind <- rbinom(n * sim, 1, ps)
  
  Yobs <- ifelse(
    T_ind == 1,
    Y1,
    Y0
  )
  
   Xtest <- Xfun(
    ntest * sim,
    d,
    indep_covariate,
    rho
  )
  
  colnames(Xtest) <- paste0("Xtest", 1:d)
  pstest <- psfun_overlap(
    1 + 4 * Xtest[, 1]
  )
  tau_test <- CATE_fun(Xtest)
  mu0_test <- beta_0(Xtest) - 0.5 * tau_test
  mu1_test <- beta_0(Xtest) + 0.5 * tau_test
  Y0test <- mu0_test + rnorm(ntest * sim)
  Y1test <- mu1_test + rnorm(ntest * sim)
  T_ind_test <- rbinom(ntest * sim,1,pstest)
  dat <- tibble(
    sim_id,
    Y1,
    Y0,
    Yobs,
    T_ind,
    CATE = tau,
    ps,
    as_tibble(X),
    Y1test,
    Y0test,
    T_ind_test,
    CATE_test = tau_test,
    pstest,
    as_tibble(Xtest)
  )
  
  return(
    dat %>%
      arrange(sim_id)
  )
}

# Generate datasets

set.seed(123)
sim_s1 <- gen_data_s1(
  sim = 1000,
  n = 1000,
  ntest = 1000,
  d = 10,
  indep_covariate = TRUE,
  rho = 0
)

saveRDS(
  sim_s1,
  here::here(
    "simulation",
    "data",
    "sim_s1.rds"
  )
)


sim_s3 <- gen_data_s3(
  sim = 1000,
  n = 1000,
  ntest = 1000,
  d = 10,
  indep_covariate = TRUE,
  rho = 0
)

saveRDS(
  sim_s3,
  here::here(
    "simulation",
    "data",
    "sim_s3_moderate_overlap_2.rds"
  )
)
