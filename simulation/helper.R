##Helper File##


#### Data generation ####


gen_data = function(sim=1, # number of simulations
                    n=1000, # number of observations
                    ntest=1000, # number of unseen test points
                    d=10, # number of independent variables
                    indep_covariate, 
                    homosc_err, 
                    rho
){
  
  sim_id <- rep(1:sim, each=n)
  
  errdist <- rnorm
  
  X <- Xfun(n*sim, d, indep_covariate, rho)
  colnames(X) <- paste("X", seq(d), sep="")
  Y0 <- rep(0, n*sim)
  tau <- taufun(X) 
  std <- sdfun(X, homosc_err)
  Y1 <- tau + std * errdist(n*sim)
  ps <- psfun(X)
  T_ind <- as.numeric(runif(n*sim) < ps)
  Yobs <- Y0
  Yobs[T_ind == 1] <- Y1[T_ind == 1]
  
  Xtest <- Xfun(ntest*sim, d, indep_covariate, rho)
  colnames(Xtest) <- paste("Xtest", seq(d), sep="")
  Y0test <- rep(0, ntest*sim)
  tautest <- taufun(Xtest) # true CATE
  stdtest <- sdfun(Xtest, homosc_err)
  Y1test <- tautest + stdtest * errdist(ntest*sim) # true ITE
  pstest <- psfun(Xtest)
  T_ind_test <- as.numeric(runif(ntest) < pstest)
 
  dat = tibble(
    sim_id, 
    Y1,
    Y0,
    tau,
    ps, 
    T_ind, 
    Yobs, 
    as.tibble(X),
    Y1test,
    Y0test,
    tautest,
    pstest, 
    T_ind_test,
    as.tibble(Xtest),
  )
  
  arrange(dat, sim_id)
  
  return(dat)
  
}

Xfun <- function(n, d, indep_covariate, rho){
  if(indep_covariate==T){
    # independent covariates
    Xfunmat <-   matrix(runif(n * d), nrow = n, ncol = d) ######### runif 
    return(Xfunmat)
  }else if(indep_covariate==F){
    # correlated covariates
    X <- matrix(rnorm(n * d), nrow = n, ncol = d)
    fac <- rnorm(n)
    X <- X * sqrt(1 - rho) + fac * sqrt(rho)
    Xfunmat <- pnorm(X)
    return(Xfunmat)
  }else{
    stop("Neither indep covariates nor corr covariates chosen.")
  }
}
 taufun <- function(X){
   2 / (1 + exp(-12 * (X[, 1] - 0.5))) * 2 / (1 + exp(-12 * (X[, 2] - 0.5)))
 }


sdfun <- function(X, homosc_err){
  if(homosc_err==T){
    return(rep(1, nrow(X)))
  }else if(homosc_err==F){
    return(-log(X[, 1] + 1e-9))
  }else{
    stop("Neither homosc. nor heterosc. errors chosen.")
  }
  
}
psfun <- function(X){
  (1 + pbeta(X[, 1], 2, 4)) / 4
}


gen_data_overlap= function(sim=1, # number of simulations
                           n=1000, # number of observations
                           ntest=1000, # number of unseen test points
                           d=10, # number of independent variables
                           indep_covariate, 
                           homosc_err, 
                           rho
){
  
  sim_id <- rep(1:sim, each=n)
  
  errdist <- rnorm
  
  X <- Xfun(n*sim, d, indep_covariate, rho)
  colnames(X) <- paste("X", seq(d), sep="")
  Y0 <- rep(0, n*sim)
  tau <- taufun(X) 
  std <- sdfun(X, homosc_err)
  Y1 <- tau + std * errdist(n*sim)
  ps <- psfun_overlap(1+4*X[, 1])
  #T_ind <- as.numeric(runif(n*sim) < ps)
  T_ind <- rbinom(n*sim,size=1, prob=ps) #sum(T_id)
  Yobs <- Y0
  Yobs[T_ind == 1] <- Y1[T_ind == 1]
  
  Xtest <- Xfun(ntest*sim, d, indep_covariate, rho)
  colnames(Xtest) <- paste("Xtest", seq(d), sep="")
  Y0test <- rep(0, ntest*sim)
  tautest <- taufun(Xtest) # true CATE
  stdtest <- sdfun(Xtest, homosc_err)
  Y1test <- tautest + stdtest * errdist(ntest*sim) # true ITE
  #pstest <- psfun(Xtest)
  pstest <-  psfun_overlap(1+4*Xtest[, 1])
  #T_ind_test <- as.numeric(runif(ntest) < pstest)
  T_ind_test <- rbinom(n*sim,size=1, prob=pstest) #sum(T_id)
  
  Yobs_test <- Y0test
  Yobs_test[T_ind_test == 1] <- Y1test[T_ind_test == 1]
  true_ITE_test <-  Y1test - Y0test
  
  dat = tibble(
    sim_id, 
    Y1,
    Y0,
    tau,
    ps, 
    T_ind, 
    Yobs, 
    as.tibble(X),
    Y1test,
    Y0test,
    tautest,
    pstest,
    T_ind_test,
    as_tibble(Xtest),
  )
  
  arrange(dat, sim_id)
  
  return(dat)
  
}

psfun_overlap <- function(X){
  1/ (1 + exp(-X))
}
X_learner<-function(Y, X, Xtest, T_ind){
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
  # Propensity score model
  m_prop <- Rforestry::forestry(x=as.matrix(X), y=T_ind,
                                ntree = 500,
                                replace = TRUE,
                                sample.fraction =  0.5,
                                mtry = ncol(X),
                                nodesizeSpl = 11,
                                nodesizeAvg = 33,
                                nodesizeStrictSpl = 2,
                                nodesizeStrictAvg = 1,
                                splitratio = 1,
                                middleSplit = FALSE,
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
  colnames(Xtest)<-colnames(X)
  prop_scores <- predict(m_prop, Xtest)
  x_cate_test <- (prop_scores * predict(mx0, as.matrix(Xtest)) +
                    (1-prop_scores) * predict(mx1, as.matrix(Xtest)))
  return(x_cate_test)
}
X_learner_fkbal<-function(Y, X, Xtest, T_ind,T_ind_test){
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
  
  kbalout_treated= kbal::kbal(allx=Xtest,sampled=T_ind_test,sampledinpop=T,
                              ebal.tol=1e-6, b=length(Xtest), printprogress =F, linkernel = T)
  kbalout_control= kbal::kbal(allx=Xtest,sampled=1-T_ind_test,sampledinpop=T,
                              ebal.tol=1e-6, b=length(Xtest), printprogress =F, linkernel = T)
  
  f_kbal<-kbalout_control$w/(kbalout_control$w+kbalout_treated$w)

  x_cate_test <- f_kbal * predict(mx0, as.matrix(Xtest)) + 
    (1-f_kbal) * predict(mx1, as.matrix(Xtest))
  
  return(x_cate_test)
}

# #### kbal functions
kbal_function<-function(X, T_ind){
  X <- as.data.frame(X)
  data <- as.data.frame(data)
  kbalout_treated= kbal::kbal(allx=X,sampled=T_ind,sampledinpop=T,
                              ebal.tol=1e-6, b=ncol(X), printprogress =F, linkernel = TRUE)
  kbalout_control= kbal::kbal(allx=X,sampled=1-T_ind,sampledinpop=T,
                              ebal.tol=1e-6, b=ncol(X), printprogress =F, linkernel = TRUE)
  weights_kbal=numeric(nrow(X))
  
  for(k in 1:nrow(X)){
    if(kbalout_treated$w[k]==1){
      weights_kbal[k]<-kbalout_control$w[k]
    }
    else{
      weights_kbal[k]<-kbalout_treated$w[k]
    }
  }
  weights_output<-weights_kbal
  return(weights_output)
}

CF_Cf_CI <- function(X, Y, T_ind, Xtest){
  fit <- grf::causal_forest(X, Y, T_ind, num.trees = 4000)
  ATE<- grf::average_treatment_effect(fit,target.sample = "all")
  preds <- predict(fit, Xtest, estimate.variance = TRUE)
  #return(pred)
  CI <- data.frame(low = preds[, 1] - 1.96 * sqrt(preds[, 2]),
                   high = preds[, 1] + 1.96 * sqrt(preds[, 2]))
return(data.frame(lower = CI$low, upper = CI$high, pred = preds[, 1],sd_pred=preds[, 2], ATE=rep(ATE[1], length(CI$low))))
}

CF_Cf_CI_W <- function(X, Y, T_ind, Xtest){
   kbalout_treated= kbal::kbal(allx=X,sampled=T_ind,sampledinpop=T,
                               ebal.tol=1e-6, b=length(X), printprogress =F, linkernel = TRUE)
   kbalout_control= kbal::kbal(allx=X,sampled=1-T_ind,sampledinpop=T,
                               ebal.tol=1e-6, b=length(X), printprogress =F, linkernel = TRUE)
   weights_kbal <- ifelse(kbalout_treated$w == 1, kbalout_control$w, kbalout_treated$w)

  fit <- grf::causal_forest(X, Y, T_ind, sample.weights =weights_kbal, num.trees = 4000) 
  ATE<- grf::average_treatment_effect(fit,target.sample = "all")
  preds <- predict(fit, Xtest, estimate.variance = TRUE)
  
  CI <- data.frame(low = preds[, 1] - 1.96 * sqrt(preds[, 2]),
                   high = preds[, 1] + 1.96 * sqrt(preds[, 2]))
  return(data.frame(lower = CI$low, upper = CI$high, pred = preds[, 1],sd_pred=preds[, 2], ATE=rep(ATE[1], length(CI$low))))
}
CF_Cf_CI_IPW <- function(X, Y, T_ind, Xtest){
  log_reg <- glm(T_ind ~ ., data = X, family = binomial)
  # Calculate propensity scores
  Propensity_Logistic = predict(log_reg, type = "response")
  min_threshold <- 0.1
  max_threshold <- 0.9
  # Apply the trimming
  Propensity_trimmed <- pmin(pmax(Propensity_Logistic, min_threshold), max_threshold)
  
  IPW<-ifelse(T_ind == 1, 1 / Propensity_Logistic, 1 / (1 - Propensity_Logistic))
  fit <- grf::causal_forest(X, Y, W=T_ind, sample.weights =IPW, num.trees = 4000) 
  ATE<- grf::average_treatment_effect(fit,target.sample = "all")
  
  preds <- predict(fit, Xtest, estimate.variance = TRUE)
  
  CI <- data.frame(low = preds[, 1] - 1.96 * sqrt(preds[, 2]),
                   high = preds[, 1] + 1.96 * sqrt(preds[, 2]))
  
  return(data.frame(lower = CI$low, upper = CI$high, pred = preds[, 1],sd_pred=preds[, 2], ATE=rep(ATE[1], length(CI$low))))
}
CF_Cf_CI_IPW_trimmed <- function(X, Y, T_ind, Xtest){
  log_reg <- glm(T_ind ~ ., data = X, family = binomial)
  # Calculate propensity scores
  Propensity_Logistic = predict(log_reg, type = "response")
  min_threshold <- 0.1
  max_threshold <- 0.9
  # Apply the trimming
  Propensity_trimmed <- pmin(pmax(Propensity_Logistic, min_threshold), max_threshold)
 IPW_trimmed<-ifelse(T_ind == 1, 1 / Propensity_trimmed, 1 / (1 - Propensity_trimmed))
  
  fit <- grf::causal_forest(X, Y, W=T_ind, sample.weights =IPW_trimmed, num.trees = 4000) 
  ATE<- grf::average_treatment_effect(fit,target.sample = "all")
  
  preds <- predict(fit, Xtest, estimate.variance = TRUE)
  
  CI <- data.frame(low = preds[, 1] - 1.96 * sqrt(preds[, 2]),
                   high = preds[, 1] + 1.96 * sqrt(preds[, 2]))
  
  return(data.frame(lower = CI$low, upper = CI$high, pred = preds[, 1],sd_pred=preds[, 2], ATE=rep(ATE[1], length(CI$low))))
}

