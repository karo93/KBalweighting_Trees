#Data generating
library(bartcs)
library(synthpop)
library(here)
library(tidyverse)
splitRule <- "X_learner_fkbal" 

Test_data <-  readRDS("~/Documents/GitHub/Clone Github/Opti_ML/00_Data/Application/IHDPA/new_test_data.rds")
Train_data <- readRDS("~/Documents/GitHub/Clone Github/Opti_ML/00_Data/Application/IHDPA/new_train_data.rds")
names_bin<- c("X4",paste0("X", 7:25))
results_data <- lapply(1:1000, function(i) {
  print(i)
  Test_data_small <- Test_data[[i]]
  Train_data_small <- Train_data[[i]]
  
  X_test <- Test_data_small[,6:30] 
  X_train <- Train_data_small[,6:30] 
  Y_test <- Test_data_small$Yobs
  Y_train <- Train_data_small$Yobs
  t_train <- Train_data_small$treatment
  t_test <- Test_data_small$treatment
  
  x_preds <- X_learner_kbal_scale_proxy(X = X_train %>% as.data.frame(),
                       Y = Y_train,
                       T_ind = t_train,
                       T_ind_test=t_test,
                       Xtest = X_test %>% as.data.frame())
  
  results <- x_preds %>% 
    as_tibble() %>% 
    mutate(true_ATE = mean(Test_data_small$mu1-Test_data_small$mu0), 
           true_ITE = Test_data_small$mu1-Test_data_small$mu0,  
           id = rep(i, nrow(Test_data_small)))
  
  return(results)
})

# Specify the path
output_path <- "/Users/karo/Documents/GitHub/Clone Github/Opti_ML/03_results/application/X_learner_fkbal/X_learner_fkbal_train.rds"

output_dir <- here::here(paste0("03_results/application/", splitRule))
output_path <- file.path(output_dir, paste0(splitRule, "_results.rds"))
saveRDS(results_data, file = output_path)

# #IPW: 
# IPW.forest <- grf::regression_forest(X_train, t_train)
# prop.hat <- predict(IPW.forest)$predictions
# IPW<-ifelse(t_train == 1, 1 / prop.hat, 1 / (1 - prop.hat))
# range(prop.hat) #0.1241357 0.2431056
# 
# data_rf<- cbind(X_train,t_train)
# rf <- randomForest::randomForest(as.factor(t_train) ~ ., data = X_train, ntree = 2000)
# 
# # Propensity scores from regression trees
# Propensity_Tree <- predict(rf, type = "prob")[, 2]
# range(Propensity_Tree) #0.0111576 0.6097240
# 
# log_reg <- glm(t_train ~ ., data = data_rf, family = binomial)
# # Calculate propensity scores
# Propensity_Logistic = predict(log_reg, type = "response")
# IPW<-ifelse(t_train == 1, 1 / Propensity_Logistic, 1 / (1 - Propensity_Logistic))
# 
# range(Propensity_Logistic) #0.001932245 0.879072624
# # Calculate RMSE for logistic regression
# rmse_logistic <- sqrt(mean((t_train - Propensity_Logistic)^2))
# 
# # Calculate RMSE for regression trees (random forest)
# rmse_tree <- sqrt(mean((t_train - Propensity_Tree)^2))
# rmse_grf <- sqrt(mean((t_train - prop.hat)^2))
# 
# # Print RMSE values
# cat("RMSE for Logistic Regression:", rmse_logistic, "\n")
# cat("RMSE for Regression Trees:", rmse_tree, "\n")


#################### Functions 
X_learner_kbal_scale_proxy<-function(Y, X, Xtest, T_ind){
  m0 <- Rforestry::forestry(x=as.matrix(X[T_ind == 0,]), y= Y[T_ind == 0],
                            ntree = 1000,
                            replace = TRUE,
                            sample.fraction = 0.8,
                            mtry = round(ncol(X) * 13 / 20),
                            nodesizeSpl = 2,
                            nodesizeAvg = 1,
                            nodesizeStrictSpl = 2,
                            nodesizeStrictAvg = 1,
                            splitratio = 1,
                            middleSplit = TRUE,
                            OOBhonest = TRUE)
  m1 <- Rforestry::forestry(x=as.matrix(X[T_ind == 1,]), y= Y[T_ind == 1],
                            ntree = 1000,
                            replace = TRUE,
                            sample.fraction = 0.8,
                            mtry = round(ncol(X) * 13 / 20),
                            nodesizeSpl = 2,
                            nodesizeAvg = 1,
                            nodesizeStrictSpl = 2,
                            nodesizeStrictAvg = 1,
                            splitratio = 1,
                            middleSplit = TRUE,
                            OOBhonest = TRUE)
  r_0 <- predict(m1,as.matrix(X[T_ind == 0,])) - Y[T_ind == 0]
  r_1 <- Y[T_ind == 1]- predict(m0,as.matrix(X[T_ind == 1,]))
  
  mx0 <- Rforestry::forestry(x=as.matrix(X[T_ind == 0,]), y= r_0,
                             ntree = 1000,
                             replace = TRUE,
                             sample.fraction = 0.7,
                             mtry = round(ncol(X) * 17 / 20),
                             nodesizeSpl = 5,
                             nodesizeAvg = 6,
                             nodesizeStrictSpl = 3,
                             nodesizeStrictAvg = 1,
                             splitratio = 1,
                             middleSplit = TRUE,
                             OOBhonest = TRUE)
  
  mx1 <- Rforestry::forestry(x=as.matrix(X[T_ind == 1,]), y= r_1,
                             ntree = 1000,
                             replace = TRUE,
                             sample.fraction = 0.7,
                             mtry = round(ncol(X) * 17 / 20),
                             nodesizeSpl = 5,
                             nodesizeAvg = 6,
                             nodesizeStrictSpl = 3,
                             nodesizeStrictAvg = 1,
                             splitratio = 1,
                             middleSplit = TRUE,
                              OOBhonest = TRUE)
   kbalout_treated <- kbal::kbal(allx = Xtest, mixed_data = TRUE, cat_columns = names_bin,
                          sampled = T_ind_test, sampledinpop = TRUE, ebal.tol = 1e-6,
                          b = length(X), printprogress = T, linkernel = FALSE,
                          scale_data = FALSE,maxnumdims = 150)
  kbalout_control <- kbal::kbal(allx = Xtest, sampled = 1 - T_ind_test, sampledinpop = TRUE,
                          mixed_data = TRUE, cat_columns = names_bin, ebal.tol = 1e-6,
                          b = length(Xtest), printprogress = T, linkernel = FALSE,
                          scale_data = FALSE,maxnumdims = 150)

#weights_kbal <- ifelse(kbalout_treated$w == 1, kbalout_control$w, kbalout_treated$w)
  
  f_kbal<-kbalout_control$w/(kbalout_control$w+kbalout_treated$w)
  #x_cate_test <- (kbalout_control$w)/2 * predict(mx0, as.matrix(Xtest)) + 
   # (kbalout_treated$w)/2 * predict(mx1, as.matrix(Xtest))
 x_cate_test <- f_kbal * predict(mx0, as.matrix(Xtest)) + 
    (1-f_kbal) * predict(mx1, as.matrix(Xtest))
  
  return(x_cate_test)
}
CF_Cf_CI_W <- function(X, Y, T_ind, Xtest){
  kbalout_treated <- kbal::kbal(allx = X, mixed_data = TRUE, cat_columns = names_bin,
                                sampled = T_ind, sampledinpop = TRUE, ebal.tol = 1e-4,
                                b = length(X), printprogress = T, linkernel = FALSE,
                                scale_data = FALSE,cont_scale = 1,maxnumdims = 140 )
  kbalout_control <- kbal::kbal(allx = X, sampled = 1 - T_ind, sampledinpop = TRUE,
                                mixed_data = TRUE, cat_columns = names_bin, ebal.tol = 1e-4,
                                b = length(X), printprogress = T, linkernel = FALSE,
                                scale_data = FALSE,cont_scale = 1,maxnumdims = 140)
  weights_kbal <- ifelse(kbalout_treated$w == 1, kbalout_control$w, kbalout_treated$w)
  
  fit <- grf::causal_forest(X, Y, T_ind, sample.weights =weights_kbal, num.trees = 4000) 
  pred <- predict(fit, Xtest, estimate.variance = TRUE)
  
  CI <- data.frame(low = pred[, 1] - 1.96 * sqrt(pred[, 2]),
                   high = pred[, 1] + 1.96 * sqrt(pred[, 2]))
  
  return(data.frame(lower = CI$low, upper = CI$high, pred = pred[, 1]))
}
true_ATE = mean(Test_data_small$mu1-Test_data_small$mu0)
true_ITE = Test_data_small$mu1-Test_data_small$mu0
RMSE_w<-caret::RMSE(results$value, results$true_ITE)
abs_ATE_w<-abs(mean(x_preds$pred)- true_ATE)
