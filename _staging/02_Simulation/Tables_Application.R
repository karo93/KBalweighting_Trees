library(tidyverse)

# Calculate the mean of each column using tidyverse and the pipe operator %>%
# mean_values_grf_ipw <- grf_ipw_log_results %>%
#   summarise(across(everything(), mean))
# 
# mean_values_grf_ipw_trimmed <- grf_ipw_trimmed_results %>%
#   summarise(across(everything(), mean))
# 
mean_values_grf <- grf_results %>%
  summarise(across(everything(), mean))
# 
# mean_values_grf_w_1_ <- grf_w_1_results %>%
#   summarise(across(everything(), mean))
# 
# mean_values_x <- xlearner_hard %>%
#   summarise(across(everything(), mean))
# 
# mean_values_x_bal <- x_learner_kbal_scale_proxy_results %>%
#   summarise(across(everything(), mean))
# 
# mean_values_x_fkbal <- X_learner_fkbal %>%
#   summarise(across(everything(), mean))

mean_values_grf <- IHDP_B_no_overlap_grf %>%
  summarise(across(everything(), mean))

mean_values_grf_ipw <- IHDP_B_no_overlap_grf_ipw %>%
  summarise(across(everything(), mean))

mean_values_grf_w <- IHDP_B_no_overlap_grf_w %>%
  summarise(across(everything(), mean))

mean_values_grf_w_x4 <- IHDP_B_no_overlap_grf_w_x4 %>%
  summarise(across(everything(), mean))

mean_values_x <- IHDP_B_no_overlap_x %>%
  summarise(across(everything(), mean))

mean_values_xfkbal <- IHDP_B_no_overlap_xfkbal %>%
  summarise(across(everything(), mean))


# Print the result


print(mean_values)
combined_results_B <- bind_rows(
  "GRF" = mean_values_grf,
  "GRF IPW" = mean_values_grf_ipw,
  #"GRF IPW trimmed" = mean_values_grf_ipw_trimmed,
  "GRF KBal" = mean_values_grf_w,
  "GRF KBal x4" = mean_values_grf_w_x4,
  "XRF" = mean_values_x,
  "XRF FKBal" = mean_values_xfkbal,
 # "XRF KBal" = mean_values_x_bal,
  .id = "Method"
)
colnames(combined_results_B) <- c("Method", "ID", "Coverage", "RMSE", "MAE","BIAS", "Abs.Error ATE")
combined_results <- combined_results %>% select(-ID)

saveRDS(combined_results_B, "combined_results_B_fkbal.rds")
combined_results_B <- xtable(combined_results_B, digits = 4)
print(combined_results_B)
