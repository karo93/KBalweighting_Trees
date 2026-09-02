# Define file paths
path_in  <- here::here("simulation", "statistics_results")
path_out <- here::here("simulation", "aggregated_results")

dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

aggregate_method_results <- function(method_folder) {
  
  file_paths <- list.files(
    here::here(path_in, method_folder),
    pattern = "\\.rds$",
    full.names = TRUE
  )
  
  data_list <- list()
  
  for (file in file_paths) {
    dataset_name <- sub("\\.rds$", "", basename(file))
    data_list[[dataset_name]] <- readRDS(file)
  }
  
  results_df <- data.frame(
    Scenario = character(),
    Coverage = numeric(),
    RMSE = numeric(),
    MAE = numeric(),
    Bias = numeric(),
    AbsBias = numeric(),
    Width = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (name in names(data_list)) {
    
    df <- data_list[[name]]
    
    results_df <- rbind(
      results_df,
      data.frame(
        Scenario = name,
        Coverage = mean(df$coverage_ITE, na.rm = TRUE),
        RMSE = mean(df$precision_RMSE, na.rm = TRUE),
        MAE = mean(df$precision_MAE, na.rm = TRUE),
        Bias = mean(df$precision_BIAS, na.rm = TRUE),
        AbsBias = mean(abs(df$precision_BIAS), na.rm = TRUE),
        Width = mean(df$width_ITE, na.rm = TRUE)
      )
    )
  }
  
  results_df[order(results_df$Scenario), ]
}

# the different methods
xrf <- aggregate_method_results("x_learner_hard")
xrf_kbal <- aggregate_method_results("x_learner_fkbal")
grf_kbal <- aggregate_method_results("grf_w")
grf <- aggregate_method_results("grf")
grf_ipw <- aggregate_method_results("grf_w_trimmed")
methods_list <- list(
  "XRF" = xrf,
  "XRF KBAL" = xrf_kbal,
  "GRF KBAL" = grf_kbal,
  "GRF" = grf,
  "GRF IPW" = grf_ipw
)
combined_results <- dplyr::bind_rows(
  methods_list,
  .id = "Method"
)

combined_results <- combined_results %>%
  dplyr::arrange(Scenario, Method)

saveRDS(
  combined_results,
  here::here(path_out, "simulation_results_summary.rds")
)

# example for S3

s3_results <- combined_results %>%
  dplyr::filter(grepl("^sim_s3", Scenario)) %>%
  dplyr::mutate(
    Overlap = dplyr::case_when(
      Scenario == "sim_s3" ~ "good",
      Scenario == "sim_s3_moderate_overlap" ~ "moderate",
      Scenario == "sim_s3_no_overlap" ~ "poor",
      TRUE ~ "unknown"
    )
  ) %>%
  dplyr::select(
    Overlap,
    Method,
    Coverage,
    RMSE,
    MAE,
    Bias,
    AbsBias,
    Width
  ) %>%
  dplyr::arrange(
    factor(Overlap, levels = c("good", "moderate", "poor")),
    Method
  )
