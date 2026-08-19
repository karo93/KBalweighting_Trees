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



