library(tidyverse)
library(here)

method <- "xlearner_hard" 
set.seed(123)

cl<-detectCores() - 2
doParallel::registerDoParallel(cl)

Test_data <- readRDS("~/00_Data/05_application/Test_data.rds")
Train_data <- readRDS("~00_Data/05_application/Train_data.rds")
names_bin<- paste0("X", 7:25)
names_bin<-c("X4",names_bin)

set.seed(123)
process_data <- function(i) {
  Test_data_small <- Test_data[[i]]
  Train_data_small <- Train_data[[i]]
  
  X_test <- Test_data_small %>% select(6:30)
  X_train <- Train_data_small %>% select(6:30)
  Y_test <- Test_data_small$Yobs
  Y_train <- Train_data_small$Yobs
  t_train <- Train_data_small$treatment
  t_test <- Test_data_small$treatment
  
  x_preds <- X_learner(X = X_train %>% as.data.frame(),
                        Y = Y_train,
                        T_ind = t_train,
                        Xtest = X_test %>% as.data.frame())
  
  results <- x_preds %>% 
    as_tibble() %>% 
    mutate(true_ATE = mean(Test_data_small$mu1-Test_data_small$mu0), 
           true_ITE = Test_data_small$mu1-Test_data_small$mu0,  
           id = rep(i, nrow(Test_data_small)))
  
  return(results)
}
tictoc::tic()
all_results <- map_dfr(1:1000, process_data)
tictoc::toc()

output_dir <- here::here(paste0("03_results/application/", method))
output_path <- file.path(output_dir, paste0(method, "_results.rds"))
saveRDS(all_results, file = output_path)



