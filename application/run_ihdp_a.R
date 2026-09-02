library(tidyverse)
library(here)

method <- "xlearner_hard" 

#Setting A
Test_data <- readRDS(here::here("application", "data", "IHDP_A_test.rds"))
Train_data <- readRDS(here::here("application", "data", "IHDP_A_train.rds"))
set.seed(123)
process_data <- function(i) {
  IHDPA_test_i <- Test_data[[i]]
  IHDPA_train_i <- Train_data[[i]]
  
  X_test <- IHDPA_test_i %>% select(6:30)
  X_train <- IHDPA_train_i %>% select(6:30)
  Y_test <- IHDPA_test_i$Yobs
  Y_train <- IHDPA_train_i$Yobs
  t_train <- IHDPA_train_i$treatment
  t_test <- IHDPA_test_i$treatment
  
  x_preds <- X_learner(X = X_train %>% as.data.frame(),
                        Y = Y_train,
                        T_ind = t_train,
                        Xtest = X_test %>% as.data.frame())
  
  results <- x_preds %>% 
    as_tibble() %>% 
    mutate(true_ATE = mean(IHDPA_test_i$mu1-IHDPA_test_i$mu0), 
           true_ITE = IHDPA_test_i$mu1-IHDPA_test_i$mu0,  
           id = rep(i, nrow(IHDPA_test_i)))
  
  return(results)
}
tictoc::tic()
all_results <- map_dfr(1:1000, process_data)
tictoc::toc()
output_dir <- here::here("application", "raw_results",  "IHDP_A", method)

output_path <- file.path(output_dir, paste0(method, "_results.rds"))
saveRDS(all_results, file = output_path)


# Note that IHDP application contains both continuous and binary covariates. We have to specify the kbal weight differently, since we have continous and binary covariates
# e.g for CF KBAl
names_bin<- c("X4",paste0("X", 7:25))

CF_Cf_CI_W <- function(X, Y, T_ind, Xtest){
  kbalout_treated <- kbal::kbal(allx = X, mixed_data = TRUE, cat_columns = names_bin,
                                sampled = T_ind, sampledinpop = TRUE, ebal.tol = 1e-4,
                                b = length(X), printprogress = T, linkernel = FALSE,
                                scale_data = FALSE,cont_scale = 1,maxnumdims = 140 )
  kbalout_control <- kbal::kbal(allx = X, sampled = 1 - T_ind, sampledinpop = TRUE,
                                mixed_data = TRUE, cat_columns = names_bin, ebal.tol = 1e-4,
                                b = length(X), printprogress = T, linkernel = FALSE,
                                scale_data = FALSE, cont_scale = 1,maxnumdims = 140)
  weights_kbal <- ifelse(kbalout_treated$w == 1, kbalout_control$w, kbalout_treated$w)
  
  fit <- grf::causal_forest(X, Y, T_ind, sample.weights =weights_kbal, num.trees = 4000) 
  pred <- predict(fit, Xtest, estimate.variance = TRUE)
  
  CI <- data.frame(low = pred[, 1] - 1.96 * sqrt(pred[, 2]),
                   high = pred[, 1] + 1.96 * sqrt(pred[, 2]))
  
  return(data.frame(lower = CI$low, upper = CI$high, pred = pred[, 1]))
}
