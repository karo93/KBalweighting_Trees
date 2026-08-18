gen_data_linear <- function(sim = 10,   # Number of simulations
                            n = 1000,  # Training set size
                            ntest = 1000,  # Test set size
                            d = 10,    # Number of features
                            sigma = 1,  # Noise level
                            indep_covariate, 
                            homosc_err, 
                            rho){  # Independent or correlated X
  
  # Simulated IDs
  sim_id <- rep(1:sim, each = n)
  
  # Generate Covariates
  X <- Xfun(n * sim, d, indep_covariate, rho)
  colnames(X) <- paste("X", seq(d), sep="")
  
  # True Treatment Effect (Linear)
  beta_tau <- runif(d, 0.5, 2)  # Random coefficients for tau
  tau <- X %*% beta_tau  # Linear function of X
  
  # True Outcome Function (Linear)
  beta_Y <- runif(d, -2, 2)  # Random coefficients for control outcome
  Y0 <- X %*% beta_Y  # Control outcome
  
  # Add treatment effect
  Y1 <- Y0 + tau + rnorm(n * sim, mean = 0, sd = sigma)
  
  # Treatment Assignment (Propensity Score)
  #ps <-  psfun_overlap(1+4*X[, 1])
  ps <- psfun(X)
  
  T_ind <- rbinom(n * sim, 1, ps)  # Treatment assigned probabilistically
  
  # Observed Outcome (Treatment or Control)
  obs <- rbinom(n*sim, size=1, prob=0.5)
  Yobs <- ifelse(T_ind == 1, Y1, Y0) 
  
  # Generate Test Data
  Xtest <- Xfun(ntest * sim, d, indep_covariate, rho)
  colnames(Xtest) <- paste("Xtest", seq(d), sep="")
  tau_test <- Xtest %*% beta_tau
  Y0test <- Xtest %*% beta_Y + rnorm(ntest * sim, mean = 0, sd = sigma)
  Y1test <- Y0test + tau_test
  #pstest <-  psfun_overlap(1+4*Xtest[, 1])
  pstest <- psfun(Xtest)
  T_ind_test <- rbinom(ntest * sim, 1, pstest)
  Yobs_test <- ifelse(T_ind_test == 1, Y1test, Y0test)
  
  # Final Dataset
  dat <- tibble(
    sim_id,
    Y1=Y1[,1], Y0=Y0[,1], tau=tau[,1], ps, T_ind, Yobs, as_tibble(X),
    Y1test=Y1test[,1], Y0test=Y0test[,1], tau_test=tau_test[,1], pstest, T_ind_test, as_tibble(Xtest)
  )
  
  return(arrange(dat, sim_id))
}

gen_data_linear_no_treatment<- function(sim = 100,   # Number of simulations
                            n = 1000,  # Training set size
                            ntest = 1000,  # Test set size
                            d = 10,    # Number of features
                            sigma = 1,  # Noise level
                            indep_covariate, 
                            homosc_err, 
                            rho){  # Independent or correlated X
  
  # Simulated IDs
  sim_id <- rep(1:sim, each = n)
  
  # Generate Covariates
  X <- Xfun(n * sim, d, indep_covariate, rho)
  colnames(X) <- paste("X", seq(d), sep="")
  
  # True Treatment Effect (Linear)
  beta_0 <- runif(d, 1, 30)  # Random coefficients for tau
  mu_0<- X %*% beta_0
  # True Outcome Function (Linear)
  Y0 <- mu_0 + + rnorm(n * sim, mean = 0, sd = sigma)
  mu_1<-mu_0
  tau<- mu_1-mu_0
  # Add treatment effect
  Y1 <- mu_1 + rnorm(n * sim, mean = 0, sd = sigma)
  
  # Treatment Assignment (Propensity Score)
  ps <-  psfun_overlap(1+4*X[, 1])
  #ps <- psfun(X)
  
  T_ind <- as.numeric(runif(n*sim) < ps)
  
  # Observed Outcome (Treatment or Control)
  Yobs <- Y0
  Yobs[T_ind == 1] <- Y1[T_ind == 1]
  
  #obs <- rbinom(n*sim, size=1, prob=0.5)
  #Yobs <- ifelse(T_ind == 1, Y1, Y0) 
  
  # Generate Test Data
  Xtest <- Xfun(ntest * sim, d, indep_covariate, rho)
  colnames(Xtest) <- paste("Xtest", seq(d), sep="")
  mu_0_test<- Xtest %*% beta_0
  mu_1_test<-mu_0_test
  Y0test <- mu_0_test + rnorm(ntest * sim, mean = 0, sd = sigma)
  Y1test <- mu_1_test + rnorm(ntest * sim, mean = 0, sd = sigma)
  pstest <-  psfun_overlap(1+4*Xtest[, 1])
  #pstest <- psfun(Xtest)
  T_ind_test <- as.numeric(runif(n*sim) < pstest)
  
  tau_test<- mu_1_test-mu_0_test
  
  # Final Dataset
  dat <- tibble(
    sim_id,
    Y1=Y1[,1], Y0=Y0[,1], tau=tau[,1], ps, T_ind, Yobs, as_tibble(X),
    Y1test=Y1test[,1], Y0test=Y0test[,1], tau_test=tau_test[,1], pstest, T_ind_test, as_tibble(Xtest)
  )
  
  return(arrange(dat, sim_id))
}
gen_data_linear_no_treatment <- gen_data_linear_no_treatment(sim=1000,n=1000, d=20, ntest = 1000, indep_covariate = T, homosc_err=T, rho=0)
saveRDS(object=gen_data_linear_no_treatment, file=here::here("00_Data/Kuenzel/sim_S4_low_unbalanced.rds"))
