mean_RMSE <- function(dimensions="low"){
  scen_names <- c("sim_expr1", "sim_expr2", "sim_expr3")
  choose_dataset<- c(paste(scen_names[1], "_", dimensions, ".rds", sep=""),
                     paste(scen_names[2], "_", dimensions,"_","overlap_unbalance", ".rds", sep=""), 
                     paste(scen_names[3], "_", dimensions,"_","overlap_balance", ".rds", sep=""))
                     
  #methods_names <- c("grf", "grf_w","grf_ipw", "xlearner_hard","xlearner_kbal_new")
  methods_names <- c("grf", "grf_ipw", "grf_w", "xlearner_hard","x_learner_kbal_scale_proxy")
  values <- data.frame(Method=NA)
  for(j in 1:length(choose_dataset)){
    for(i in 1:length(methods_names)){
      dat_temp <- readRDS(file=here::here(paste("04_results_statistics/", methods_names[i], "/", choose_dataset[j], sep="")))
      values[i,j+1] <- dat_temp %>% summarise(m_RMSE = mean(precision_RMSE)) %>% pull(m_RMSE) 
    }
  }
  values[,1] <- methods_names
  colnames(values) <- c("method", "overlap", "weak overlap", "no overlap")
  #values$method <- c("CF", "CF kbal","CF IPW", "XRF", "XRF weighted")
  values$method <- c("CF","CF IPW trimmed", "CF KBal", "XRF", "XRF KBal")
  saveRDS(values, file=here::here(paste("04_results_statistics/00_tables/", dimensions, ".rds", sep="")))
  return(values)
}
mean_RMSE_low<- mean_RMSE(dimensions = "low")
library(xtable)
# mean_RMSE_high_table <- xtable(mean_RMSE_high, digits = 4)
#  print(mean_RMSE_high_table)
# ############
#  
#  mean_bias_high<- mean_bias(dimensions = "high")
#  ?xtable
#  mean_bias_high <- as.data.frame(mean_bias_high)
#  mean_bias_high_table <- xtable(mean_bias_high, digits = 4)
 print(mean_RMSE_low)
 ##################
mean_bias <- function(dimensions="low"){
  scen_names <- c("sim_expr1", "sim_expr2", "sim_expr3")
  choose_dataset<- c(paste(scen_names[1], "_", dimensions, ".rds", sep=""),
                     paste(scen_names[2], "_", dimensions,"_","overlap_unbalance", ".rds", sep=""), 
                     paste(scen_names[3], "_", dimensions,"_","overlap_balance", ".rds", sep=""))
  
  #methods_names <- c("grf", "grf_w","grf_ipw", "xlearner_hard","xlearner_kbal_new")
  methods_names <- c("grf", "grf_ipw_trimmed", "grf_w", "xlearner_hard","x_learner_kbal_scale_proxy")
  values <- data.frame(Method=NA)
  for(j in 1:length(choose_dataset)){
    for(i in 1:length(methods_names)){
      dat_temp <- readRDS(file=here::here(paste("04_results_statistics/", methods_names[i], "/", choose_dataset[j], sep="")))
      values[i,j+1] <- dat_temp %>% summarise(m_bias = mean(precision_BIAS)) %>% pull(m_bias) 
    }
  }
  values[,1] <- methods_names
  colnames(values) <- c("method", "overlap", "weak overlap", "no overlap")
  #values$method <- c("CF", "CF weighted", "CF IPW trimmed", "XRF", "XRF weighted")
  values$method <- c("CF","CF IPW trimmed", "CF KBal", "XRF", "XRF KBal")
  saveRDS(values, file=here::here(paste("04_results_statistics/00_tables/", dimensions, ".rds", sep="")))
  return(values)
}
mean_bias_low<-mean_bias(dimensions = "low")
mean_coverage <- function(dimensions="low"){
  scen_names <- c("sim_expr1", "sim_expr2", "sim_expr3")
  choose_dataset<- c(paste(scen_names[1], "_", dimensions, ".rds", sep=""),
                     paste(scen_names[2], "_", dimensions,"_","overlap_unbalance", ".rds", sep=""), 
                     paste(scen_names[3], "_", dimensions,"_","overlap_balance", ".rds", sep=""))
  
  #methods_names <- c("grf", "grf_w","grf_ipw", "xlearner_hard","xlearner_kbal_new")
  methods_names <- c("grf", "grf_ipw_trimmed", "grf_w")
  values <- data.frame(Method=NA)
  for(j in 1:length(choose_dataset)){
    for(i in 1:length(methods_names)){
      dat_temp <- readRDS(file=here::here(paste("04_results_statistics/", methods_names[i], "/", choose_dataset[j], sep="")))
      values[i,j+1] <- dat_temp %>% summarise(m_cov = mean(coverage)) %>% pull(m_cov) 
    }
  }
  values[,1] <- methods_names
  colnames(values) <- c("method", "overlap", "weak overlap", "no overlap")
  #values$method <- c("CF", "CF kbal","CF IPW", "XRF", "XRF weighted")
  values$method <- c("CF","CF IPW trimmed", "CF KBal")
  saveRDS(values, file=here::here(paste("04_results_statistics/00_tables/", dimensions, ".rds", sep="")))
  return(values)
}
mean_width <- function(dimensions="low"){
  scen_names <- c("sim_expr1", "sim_expr2", "sim_expr3")
  choose_dataset<- c(paste(scen_names[1], "_", dimensions, ".rds", sep=""),
                     paste(scen_names[2], "_", dimensions,"_","overlap_unbalance", ".rds", sep=""), 
                     paste(scen_names[3], "_", dimensions,"_","overlap_balance", ".rds", sep=""))
  
  #methods_names <- c("grf", "grf_w","grf_ipw", "xlearner_hard","xlearner_kbal_new")
  methods_names <- c("grf", "grf_ipw_trimmed", "grf_w")
  values <- data.frame(Method=NA)
  for(j in 1:length(choose_dataset)){
    for(i in 1:length(methods_names)){
      dat_temp <- readRDS(file=here::here(paste("04_results_statistics/", methods_names[i], "/", choose_dataset[j], sep="")))
      values[i,j+1] <- dat_temp %>% summarise(m_width = mean(width)) %>% pull(m_width) 
    }
  }
  values[,1] <- methods_names
  colnames(values) <- c("method", "overlap", "weak overlap", "no overlap")
  #values$method <- c("CF", "CF kbal","CF IPW", "XRF", "XRF weighted")
  values$method <- c("CF","CF IPW trimmed", "CF KBal")
  saveRDS(values, file=here::here(paste("04_results_statistics/00_tables/", dimensions, ".rds", sep="")))
  return(values)
}

mean_coverage_low<-mean_coverage(dimensions = "low")
mean_width_low<-mean_width(dimensions = "low")

# Save an object to a file
save(mean_bias_low, file = "mean_bias_low.RData")
save(mean_RMSE_low, file = "mean_RMSE_low.RData")
save(mean_coverage_low, file = "mean_coverage_low.RData")
save(mean_RMSE_low, file = "mean_RMSE_low.RData")

# Restore the object
load(file = "my_data.RData")

