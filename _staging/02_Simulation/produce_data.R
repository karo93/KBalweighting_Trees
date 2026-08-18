source(file=here::here("02_Simulation/packages.R"))
source(file=here::here("02_Simulation/helper.R"))


set.seed(1234)


#### params ####

sim_all <- 1000
d_low <- 10
d_high <- 100

#### low-dimensional ####


## Expr 1 (equation (27) of Wager and Athey (2018)): homoscedastic errors + independent covariates 
gen_data(sim=sim_all, n=1000, ntest=1000, d=d_low,
         indep_covariate = T, 
         homosc_err = T, 
         rho=0.9) %>%
  saveRDS(file=here::here("00_Data/sim_expr1_low.rds"))

## Expr 2 (equation (27) of Wager and Athey (2018)) heteroscedastic errors + independent covariates
gen_data(sim=sim_all, n=1000, ntest=1000, d=d_low,
         indep_covariate = T, 
         homosc_err = F, 
         rho=0.9) %>%
  saveRDS(file=here::here("00_Data/sim_expr2_low.rds"))

#overlap
gen_data_overlap(sim=sim_all, n=1000, ntest=1000, d=d_low,
         indep_covariate = T, 
         homosc_err = T, 
         rho=0.9) %>%
saveRDS(file=here::here("00_Data/sim_expr1_low_overlap.rds"))

gen_data_overlap(sim=sim_all, n=1000, ntest=1000, d=d_low,
         indep_covariate = T, 
         homosc_err = F, 
         rho=0.9) %>%
  saveRDS(file=here::here("00_Data/sim_expr2_low_overlap.rds"))

#IV

simulated_data <- simulate_data_iv(sim=1000, # number of simulations
                 n=1000, # number of observations
                 ntest=1000, # number of unseen test points
                 d=10,
                 model_type=5,
                 error_dist="exp") 
file_name <- paste0("sim_expr_", model_type, "_", error_dist, ".rds")
saveRDS(simulated_data, file = here::here(paste0("00_Data/IV/", file_name)))

# IV&Overlap
simulated_data <- simulate_data_iv_overlap(sim=1000, # number of simulations
                                   n=1000, # number of observations
                                   ntest=1000, # number of unseen test points
                                   d=10,
                                   model_type=5, 
                                   error_dist="uniform") 
file_name <- paste0("sim_expr_no_overlap_", model_type, "_", error_dist, ".rds")
saveRDS(simulated_data, file = here::here(paste0("00_Data/IV/no_overlap/", file_name)))

gen_data_iv(sim=sim_all, n=1000, ntest=1000, d=d_low,
         indep_covariate = T, 
         homosc_err = F, 
         rho=0.9) %>%
  saveRDS(file=here::here("00_Data/Final/IV/sim_expr1_low_iv_hete.rds"))

gen_data_overlap_iv(sim=sim_all, n=1000, ntest=1000, d=d_low,
            indep_covariate = T, 
            homosc_err = F, 
            rho=0.9) %>%
  saveRDS(file=here::here("00_Data/Final/IV/sim_expr_low_iv_rnorm_hete.rds"))
# 
# 
# ## Expr 3 (equation (27) of Wager and Athey (2018)), homoscedastic errors + correlated covariates 
# gen_data(sim=sim_all, n=1000, ntest=1000, d=d_low,
#          indep_covariate = F, 
#          homosc_err = T, 
#          rho=0.9) %>%
#   saveRDS(file=here::here("00_data/01_simulations/sim_expr3_low.rds"))
# 
# 
# ## Expr 4 (equation (27) of Wager and Athey (2018)) heteroscedastic errors + correlated covariates 
# gen_data(sim=sim_all, n=1000, ntest=1000, d=d_low,
#          indep_covariate = F, 
#          homosc_err = F, 
#          rho=0.9) %>%
#   saveRDS(file=here::here("00_data/01_simulations/sim_expr4_low.rds"))
# 

#### high-dimensional ####


## Expr 1 (equation (27) of Wager and Athey (2018)): homoscedastic errors + independent covariates 
gen_data(sim=sim_all, n=1000, ntest=1000, d=d_high,
         indep_covariate = T, 
         homosc_err = T, 
         rho=0.9) %>%
  saveRDS(file=here::here("00_data/01_simulations/sim_expr1_high.rds"))

## Expr 2 (equation (27) of Wager and Athey (2018)) heteroscedastic errors + independent covariates
gen_data(sim=sim_all, n=1000, ntest=1000, d=d_high,
         indep_covariate = T, 
         homosc_err = F, 
         rho=0.9) %>%
  saveRDS(file=here::here("00_data/01_simulations/sim_expr2_high.rds"))
################

gen_data_1 = function(n=1000, # number of observations
                      d=20, # number of independent variables
                      d_mu=10, # number of independent variables, which effect mu, default=8
                      prop.score=0.5, 
                      alph = 0# if 0, X's are independent, if alph>0, X's are correlated  
){
  
  
  X <- Xfunc_2(d=d, alph=alph, n=n)
  beta <- runif(d, min=-5, max=5)
  mu_0 <- X[,1:d_mu] %*% beta[1:d_mu] + as.vector(ifelse(X[, 1] > 0.5, 5, 0))
  CATE_true <- as.vector(ifelse(X[, 2] > 0.1, 8, 0))
  mu_1 <- mu_0 + CATE_true
  
  
  Y1 <- mu_1 + rnorm(n)
  Y0 <- mu_0 + rnorm(n)
  T_ind <- rbinom(n, size=1, prob=prop.score)
  Yobs <- ifelse(T_ind == 1, Y1, Y0)
  
  
  dat = tibble(
    Y1=Y1[,1],
    Y0=Y0[,1],
    Yobs,
    T_ind, 
    CATE_true,
    ps=prop.score, 
    as_tibble(X)
  )
  
  return(dat)
  
}


# uses the C-vine method for
#'   simulating correlation matrices.
#'   Lewandowski D, Kurowicka D, Joe H (2009) Generating random correlation matrices based on vines and extended onion method.
#'   https://github.com/soerenkuenzel/causalToolbox/blob/master/R/ExmpleSetups_causaleffects.R
simulate_correlation_matrix <- function(dim, alpha) {
  
  betaparam <- 1 / alpha
  
  P <- matrix(nrow = dim, ncol = dim)           # storing partial correlations
  S <- diag(dim)
  
  for (k in 1:(dim - 1)) {
    for (i in (k + 1):dim) {
      P[k, i] <- rbeta(1, betaparam, betaparam) # sampling from beta
      P[k, i] <- (P[k, i] - 0.5) * 2     # linearly shifting to [-1, 1]
      p <- P[k, i]
      if (k > 1) {
        for (l in (k - 1):1) {
          # converting partial correlation to raw correlation
          p <- p * sqrt((1 - P[l, i] ^ 2) * (1 - P[l, k] ^ 2)) + P[l, i] *
            P[l, k]
          # p  = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
        }
      }
      S[k, i] <- p
      S[i, k] <- p
    }
  }
  
  # permuting the variables to make the distribution permutation-invariant
  permutation <- sample(1:dim)
  S <- S[permutation, permutation]
  return(S)
}

Xfunc_2 <- function(d, alph, n){
  corr.mat <- simulate_correlation_matrix(dim=d, alpha=alph)
  x <- mvtnorm::rmvnorm(n, sigma=corr.mat)
  colnames(x) <- paste("X", seq(d), sep="")
  return(x)
}

sim=1000
n=1000

set.seed(1234)
train <- replicate(sim, gen_data_1(n=n, d=10, d_mu=8, prop.score = 0.1, alph = 0), simplify=FALSE)
sim_id <- tibble(sim_id=rep(1:sim, each=n))
train <- bind_rows(train)
train <- bind_cols(sim_id, train)
saveRDS(object=train, file=here::here("00_Data/S1/Unbalanced/sim_low_s1_train.rds"))

test <- replicate(sim, gen_data_1(n=n, d=10, d_mu=8, prop.score = 0.1, alph = 0), simplify=FALSE)
sim_id <- tibble(sim_id=rep(1:sim, each=n))
test <- bind_rows(test)
test <- bind_cols(sim_id, test)
saveRDS(object=test, file=here::here("00_Data/S1/Unbalanced/sim_low_s1_test.rds"))

