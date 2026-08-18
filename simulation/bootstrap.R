method <- "x_learner_fkbal" 
dataset_names <- list.files(here::here("simulation", "data" ))

# get bootstrap preds 

boot_pred <- function(b, treat_idx, contr_idx, n_1, n_0, dataTrain, dataTest, subsize_factor=1){
  
  samp_treat <- sample(x = treat_idx, size = n_1*subsize_factor, replace = T)
  samp_contr <- sample(x = contr_idx, size = n_0*subsize_factor, replace = T)
  samp_full <- c(samp_treat, samp_contr)
  samp_dataTrain <- dataTrain[samp_full,]
  samp_preds <- if (method == "xlearner_hard") {
      samp_preds<- X_learner(Y=samp_dataTrain$Yobs,
                             X=samp_dataTrain %>% select("X1":last_col()) %>% as.data.frame(),
                             Xtest = dataTest %>% select("Xtest1":last_col()) %>% as.data.frame(),
                             T_ind=samp_dataTrain$T_ind)
    }
  else if (method == "x_learner_fkbal") {
    samp_preds<- X_learner_fkbal(Y=samp_dataTrain$Yobs,
                                X=samp_dataTrain %>% select("X1":last_col()) %>% as.data.frame(),
                                Xtest = dataTest %>% select("Xtest1":last_col()) %>% as.data.frame(),
                                T_ind=samp_dataTrain$T_ind,
                                T_ind_test=dataTest$T_ind_test)
  }
  return(samp_preds)
}

bootstrap_ci <- function(dataset_names, dataset_num, sim_id_num, B_boot, subsize_fac){
  
  # load dataset
  dat1<-readRDS(here::here("simulation","data", dataset_names[dataset_num] ))
  results_dat<-readRDS(here::here("simulation","raw_results", method , dataset_names[dataset_num]))
  # prepare data
  Xcol_select <- floor(((ncol(dat1)+3)/2))
  Xtestcol_select <- Xcol_select +1
  
  datTrain <- dat1 %>% filter(sim_id==sim_id_num) %>%
    select(2:Xcol_select)
  
  datTest <- dat1 %>% filter(sim_id==sim_id_num) %>%
    select(Xtestcol_select:ncol(dat1))
  
  treat_idx <- which(datTrain$T_ind %in% 1)
  contr_idx <- which(datTrain$T_ind  %in% 0)
  n_1 <- length(treat_idx)
  n_0 <- length(contr_idx)

  # get full sample point predictions
  preds_full_sample <- results_dat[[sim_id_num]]$value
  true_CATE <- results_dat[[sim_id_num]]$true_CATE
  
  preds_bootstrap <- future.apply::future_sapply(1:B_boot, boot_pred,
                                                 treat_idx=treat_idx,
                                                 contr_idx=contr_idx,
                                                 n_1=n_1,
                                                 n_0=n_0, 
                                                 dataTrain=datTrain,
                                                 dataTest=datTest,
                                                 subsize_factor=subsize_fac,
                                                 future.seed=123,
                                                 future.packages = c("rpart")
  )
  # get bootstrapped sds
  sd_preds <- apply(preds_bootstrap, 1, sd) 
  boot_mean_preds <- apply(preds_bootstrap, 1, mean) 
  CIalpha <- 0.05
  z_score <- qnorm(1-CIalpha/2)
  lower <-  preds_full_sample-z_score*sd_preds
  upper <- preds_full_sample+z_score*sd_preds

  results <- sd_preds %>%
    as_tibble() %>% rename(sd_preds = value) %>%
    mutate(
      lower = lower,
      upper = upper,
      pred = preds_full_sample,
      pred_boot = boot_mean_preds,
      true_CATE = true_CATE,
      covered = (lower <= true_CATE & true_CATE <= upper)
    )
  return(results)
}

# to iterate over multiple sim_ids
#future::plan(multisession, workers = 7)
#options(future.globals.maxSize=5e+08)
#future::plan(sequential)

multiple_mapply <- function(data_id, sims_n, boots_n, subsize_fac){ 
  mapply(bootstrap_ci,
         sim_id_num = sims_n,
         MoreArgs = list(dataset_names=dataset_names,
                         B_boot = boots_n, 
                         dataset_num = data_id,
                         subsize_fac = subsize_fac))
}

tictoc::tic() # to track time

Data_Boot <- lapply(seq_along(dataset_names), # ~ 20 datasets currently
                    multiple_mapply, # see function above
                    sims_n=1:1000, # sim_ids, could iterate over all by 1:1000 
                    boots_n=50, # bootstrap reps
                    subsize_fac=1 # subsize factor resampling
)
tictoc::toc()
names(Data_Boot) <- dataset_names
out_dir <- here::here("simulation","bootstrap_results", method )
dir.create( out_dir, recursive = TRUE, showWarnings = FALSE )
saveRDS( Data_Boot, file = file.path(out_dir, "list_boot_50.rds"))

as_sim_list <- function(x) {
  stopifnot(is.list(x), !is.null(dim(x)))
  dn <- dimnames(x)
  
  # Require rownames
  if (is.null(dn) || is.null(dn[[1]])) {
    stop("x must have dimnames with row names (e.g., 'pred', 'lower', ...).")
  }
  
  rows <- dn[[1]]
  nsim <- dim(x)[2]
  
  # Helper: safely extract (row, col) and return the underlying vector
  get_cell <- function(r, j) {
    if (!r %in% rows) return(NULL)
    v <- x[r, j][[1]]
    v
  }
  
  # Decide expected length per sim from pred (fallback: first available row)
  get_n <- function(j) {
    v <- get_cell("pred", j)
    if (!is.null(v)) return(length(v))
    # fallback: find first non-null cell in that column
    for (r in rows) {
      v2 <- get_cell(r, j)
      if (!is.null(v2)) return(length(v2))
    }
    stop("Could not determine vector length for sim ", j)
  }
  
  lapply(seq_len(nsim), function(j) {
    n <- get_n(j)
    
    # Extract and standardize lengths (replicate scalars; error on weird lengths)
    pull <- function(r) {
      v <- get_cell(r, j)
      if (is.null(v)) return(rep(NA, n))
      if (length(v) == 1L && n > 1L) return(rep(v, n))
      if (length(v) != n) stop("Length mismatch for row '", r, "' in sim ", j,
                               ": got ", length(v), " expected ", n)
      v
    }
    
    data.frame(
      pred      = pull("pred"),
      lower     = pull("lower"),
      upper     = pull("upper"),
      sd_preds  = pull("sd_preds"),
      pred_boot = pull("pred_boot"),
      true_CATE = pull("true_CATE"),
      covered   = pull("covered"),
      sim_id    = rep.int(j, n),
      check.names = FALSE
    )
  })
}
Data_Boot_simlists <- lapply(Data_Boot, as_sim_list)
idx <- 1:3
names(Data_Boot_simlists) <- dataset_names[idx]  # nice to have
for (j in seq_along(idx)) {
  i <- idx[j]
  saveRDS(
    Data_Boot_simlists[[j]],
    file = file.path(out_dir, paste0("Bootstrap_50_", dataset_names[i]))
  )
}
