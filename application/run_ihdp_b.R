library(tidyverse)
library(here)

method <- "xlearner_hard" 

#Setting B
dat1 <- readRDS(
  here::here("application", "data", "IHDP_B.rds"))

set.seed(123)

results_data <- lapply(1:1000, function(i) {
  data_i <- dat1 %>%
    dplyr::filter(iteration == i)
  # Random train/test split
  test_id <- sample(
    seq_len(nrow(data_i)),
    floor(nrow(data_i) / 2)
  )
  Train <- data_i[-test_id, ]
  Test  <- data_i[test_id, ]
  
  # Covariates and outcomes
  X_train <- Train[, 2:26]
  X_test  <- Test[, 2:26]
  Y_train <- Train$y
  Y_test  <- Test$y
  t_train <- Train$z
  t_test  <- Test$z
  
  # Estimate CATE
  x_preds <- CF_Cf_CI_IPW(
    X = X_train %>% as.data.frame(),
    Y = Y_train,
    T_ind = t_train,
    Xtest = X_test %>% as.data.frame()
  )
  
  # Store results
  x_preds %>% as_tibble() %>% mutate(
      true_ATE = mean(Test$mu_1 - Test$mu_0),
      true_ITE = Test$mu_1 - Test$mu_0,
      iteration = i
    )
})

# Combine replications
all_results <- dplyr::bind_rows(results_data)

# Save raw results
output_dir <- here::here(
  "application",
  "raw_results",
  "IHDP_B",
  method
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  all_results,
  file = file.path(
    output_dir,
    paste0(method, "_results.rds")
  )
)

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
