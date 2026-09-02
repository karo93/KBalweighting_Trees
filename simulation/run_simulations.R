# load packages and functions
source(file = here::here("simulation", "packages.R"))
source(file = here::here("simulation", "helper.R"))

method <- "grf_ipw" # change method here
dataset_names <- list.files(here::here("simulation", "data"))
set.seed(123)

for(j in 1:length(dataset_names)){
  dat1 <- readRDS(here::here("simulation", "data", dataset_names[j]))
  lapply(1:1000,
         function(i){
           # prepare data
           Xcol_select <- floor(((ncol(dat1)+3)/2))
           Xtestcol_select<-Xcol_select+1
           datTrain <- dat1 %>% filter(sim_id==i) %>% select(2:Xcol_select)
           datTest <- dat1 %>% filter(sim_id==i) %>% select(Xtestcol_select:last_col())
           x_preds<- CF_Cf_CI_IPW_trimmed(X = datTrain %>% select(X1:last_col()) %>% as.data.frame(),
                                Y = datTrain$Yobs,
                                T_ind = datTrain$T_ind,
                                Xtest = datTest %>% select(Xtest1:last_col())  %>% as.data.frame())
                                #T_ind_test= datTest$T_ind)
           
           
           # change x_preds to the method you are using
           
           results <- x_preds %>% 
             as_tibble() %>% 
             mutate(true_CATE=datTest$CATE_test, id=rep(i, length(datTest$CATE_test)))
           return(results)
         }
  ) %>% 
    saveRDS(file = here::here("simulation", "raw_results", method, dataset_names[j]))
}
