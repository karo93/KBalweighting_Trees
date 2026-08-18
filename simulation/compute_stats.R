# load dataset by choosing the datasetname and splitrule

method <- "xlearner_hard"
dataset_names <- list.files(here::here("simulation", "raw_results", method ))

lapply(dataset_names, 
       function(k){
         choose_dataset <- k
         dataset<-readRDS(here::here("simulation","raw_results", method , choose_dataset ))
         dplyr::bind_rows(dataset) %>%
          mutate(id=rep(1:length(dataset), each=nrow(dataset[[1]]))) %>%
          group_by(id) %>% 
          summarise(#coverage_ITE = mean(true_CATE > lower & true_CATE < upper), 
             precision_RMSE = caret::RMSE(value, true_CATE), # sometime value, with CF its "pred"
             precision_MAE =mean(abs(value - true_CATE)),
             precision_BIAS =mean((value - true_CATE))
          #   precision_abs_ATE = abs(mean(value)-mean(true_ATE))) # commented out if Metalearner
        ) %>%
           saveRDS(file=here::here("simulation","statistics", method, choose_dataset))
       })
