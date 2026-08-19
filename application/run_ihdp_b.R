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