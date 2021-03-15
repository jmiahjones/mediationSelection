#
# Author: Jeremiah Jones
# rm(list=objects())

main <- function(
  n,
  num_simulations,
  abbv_scn,
  coef_setting,
  weight_gam,
  use_sl,
  suffix_arg,
  cores
) {
  
  assertthat::assert_that(is.logical(use_sl))
  
  ######################## Load Procedures ########################
  # need this in the function environment
  source("./R/simulation-helpers.R", local=T)
  source("./R/mediation-funs-postsel.R", local=T)
  
  num_noise = 7
  
  num_bootstraps=100L
  small.setting <- substr(coef_setting, 1,1)=="s"
  is.randomized <- substr(abbv_scn, 1,1)=="r"
  
  if(small.setting){
    # small
    # chose alphas*betas=16/sqrt(n) for the first 3 mediators
    alphas = c(4, 4*n^(1/4), 4*n^(1/4), rep(0, num_noise)) # D-M coef
    betas = c(4, 4*n^(-1/4), 4*n^(-1/4), rep(0,num_noise)) # M-Y coef
    alphas = alphas/(n^(1/4))
    betas = betas/(n^(1/4))
  } else {
    # nonsmall
    # chose alphas*betas=.8 + O(1/n) for the first 3 mediators
    alpha0 <- c(1, 2, 2, rep(0, num_noise))
    beta0 <- c(.8, .4, .4, rep(0, num_noise))
    # alpha.h = c(rep(1, 3), rep(1, num_noise)) # D-M coef
    # beta.h = c(rep(1, 3), rep(0, num_noise)) # M-Y coef
    alphas = alpha0 #+ alpha.h/(n^(1/2))
    betas = beta0 #+ beta.h/(n^(1/2))
  }
  
  
  
  noise = c(0,0,0,rep(0,num_noise))
  variances <- c(1,1,1, rep(1, num_noise-1), 1)
  NDE = 2
  NIE = sum(alphas * betas)
  rho = 0.0 # covariance between errors
  rho2 = 0.0 # noise covariance
  corr_meds = c(1,9)
  p = length(alphas)
  covariates_size = 3
  V = 10 # number of folds for crossvalidation
  
  
  
  ################# Initialization ########################
  
  start = Sys.time()
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(SuperLearner)
  library(glmnet)
  # library(gam)
  library(mgcv)
  library(earth)
  library(randomForest)
  # library(e1071)
  # library(xgboost)
  library(parallel)
  library(doParallel)
  library(foreach)
  # print(sessionInfo())
  
  
  candidate.mediators = paste0("m.", 1:p)
  true_M <- which(alphas*betas != 0)
  oracle.mediators <- candidate.mediators[true_M]
  m.cols="m"
  # x.cols = paste0("x.", 1:covariates_size)
  x.cols="x"
  expit <- plogis
  
  # Setup a lambda grid
  lambda.len = 301
  lambdas = n^(-0.75)*(2^seq(from=-2, to=10, length.out=lambda.len))
  
  set.seed(841664)
  folds <- caret::createFolds(y=1:n, k=V)
  
  ################# Files ########################
  
  savefile_suffix <- paste(n, num_simulations, abbv_scn, coef_setting,
                           weight_gam_str, use_sl, suffix_arg, sep="-")
  # print(savefile_suffix)
  # slfile_y_suffix <- paste(n, num_simulations, abbv_scn, coef_setting,
  #                          suffix_arg, sep="-")
  # slfile_md_suffix <- paste(n, num_simulations, substr(abbv_scn, 1,2), 
  #                           coef_setting,
  #                           suffix_arg, sep="-")
  
  logfile <- paste0("./logs/cluster-", savefile_suffix, ".out")
  # savefile <- paste0("./cache/", savefile_suffix, "-var-selection.RData")
  # slfile_y <- paste0("./cache/sl/", slfile_y_suffix, "-train-y.RData")
  # slfile_md <- paste0("./cache/sl/", slfile_md_suffix, "-train-md.RData")
  
  if(!dir.exists("./logs")){
    dir.create("./logs")
  }
  if(file.exists(logfile)){
    file.remove(logfile)
  }
  
  
  
  #################### Confounding Mechanisms ########################
  
  
  if(substr(abbv_scn, 1,1)=="n"){
    # nonlinear function
    propensity <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      # expit((0*x1 + 1*x2 + 0*x3 + 0*x4 + 0*x5+0*x1^2*x2+ 1*x1*x3 - 0)/1)
      expit((x2 + x1*x2)/1.25)
    }
  } else if(substr(abbv_scn, 1,1)=="l"){
    # linear function
    propensity <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      expit((1*x1 + 1*x2 )/1.25)
    }
  } else if(substr(abbv_scn, 1,1)=="r"){
    # linear function
    propensity <- function(x) {
      0.5
    }
  }
  
  
  if(substr(abbv_scn, 2,2)=="n"){
    # nonlinear function
    psi_m <- function(x) { 
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      1*x1^2 + x2 - x3 # 0*x3 - 0*x4 - 0*x5 + 0*x2*(x1-0.5)^2
    }
  } else if(substr(abbv_scn, 2,2)=="l"){
    # linear function
    psi_m <- function(x) { 
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      1*x1 + x2 - x3
    }
  }
  
  
  
  if(substr(abbv_scn, 3,3)=="n"){
    # nonlinear function
    psi_y <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      1*(2*(x1-0.5)^2 + x2 + 2*x3)
    }
  } else if(substr(abbv_scn, 3,3)=="l"){
    # linear function
    psi_y <- function(x) {
      x1 <- x[1]
      x2 <- x[2]
      x3 <- x[3]
      x4 <- x[4]
      x5 <- x[5]
      (2*(x1-0.5) + x2 + 2*x3)
    }
  }
  
  
  ######################## Simulations ########################
  
  if(cores > 1){
    cl <- parallel::makeCluster(cores, "FORK", outfile=logfile)
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl, 2018)
  } else {
    foreach::registerDoSEQ()
  }
  
  
  
  # tryCatch({
  simulations <- foreach(sim.idx=1:num_simulations) %dopar% {
    
    #### Simulation mechanism ####
    
    # set seed to ensure similar xdm each time -- allows for saving sl
    set.seed(65411 + sim.idx)
    # confounders
    # x <- replicate(covariates_size, runif(n, min=0, max=1))
    x <- replicate(covariates_size, rnorm(n, sd=0.5))
    
    # treatment confounding
    d.prob <- apply(x, 1, function(row){
      propensity <- propensity(row)
    })
    d <- rbinom(n=n, size=1, prob=d.prob)
    
    Sigma <- 
      matrix(c(rep(rho, p**2)), ncol=p)
    
    Sigma[corr_meds, p] <- rho2
    Sigma[p, corr_meds] <- rho2
    Sigma[corr_meds, corr_meds] <- rho2
    diag(Sigma) <- rep(1, p)
    
    Sigma <- diag(sqrt(variances)) %*% Sigma %*%
      diag(sqrt(variances))
    # if(p>5){
    #   Sigma[,as.logical(noise)] <- rho2
    #   Sigma[as.logical(noise),] <- rho2
    # }
    
    epsilon_m <- MASS::mvrnorm(n, mu=rep(0, p), Sigma=Sigma)
    
    # mediator confounding
    confounding_m <- apply(x, 1, function(row){
      psi_m(row)
    })
    # browser()
    # mediators
    m <- foreach(i=1:n, .combine=rbind) %do%{
      d[i] * alphas + confounding_m[i] + epsilon_m[i,]
    }
    # m[,p] <- m[,-p] %>% rowMeans
    # m[,p] <- m[,1]^2 + m[,2]^2 + rowMeans(m[,4:(p-1)])
    # m <- m + epsilon_m
    
    # outcome confounding
    confounding_y <- apply(x, 1, function(row){
      psi_y(row)
    })
    
    set.seed(7243957 + sim.idx)
    # outcome equation
    epsilon_y <- rnorm(n, mean=0, sd=1)
    y <- (NDE * d) + (m %*% betas) + (confounding_y) + epsilon_y
    
    # for troubleshooting only
    true_em <- foreach(i=1:n, .combine=rbind) %do%{
      (d.prob[i] * alphas) + confounding_m[i]
    }
    true_ey <- as.numeric(NDE*d.prob + true_em %*% as.matrix(betas) + confounding_y)
    
    stopifnot(max(abs(
      true_ey - confounding_y - NDE*d.prob - true_em %*% betas
    )) < 1e-10)
    
    # if(sim.idx %% 10 == 0)
    #   print(paste("Completed simulation", sim.idx))
    
    # save the true centering variables for future use
    # em = replicate(p,confounding_m) + replicate(p,propensity)*alphas
    # ey = NDE*propensity + (em %*% betas) + confounding_y
    
    # return(
    tibble(
      x=x, d=d, m=m, y=y,
      propensity=d.prob,
      true_em = true_em,
      true_ey = true_ey,
      psi_mx = confounding_m,
      psi_yx = confounding_y,
      epsilon = epsilon_y,
      epsilon_m = epsilon_m
    )
    # )
  }
  # return(simulations)
  print("Created data.")
  # parallel::clusterExport(cl, "simulations")
  
  
  #### Estimate mu (if necessary) ####
  if(use_sl){
    source("./R/estimate-nuisance.R", local=T)
    # returns mu.hats
    # mu.hats <- estimate_nuisance()
    print("Estimated mu-hats with SuperLearner.")
  } else {
    print("Robinsonizing with True Mu functions.")
  }
  
  
  opt <- tibble(
    # n = n, # already calculated
    num_simulations = num_simulations,
    scenario = abbv_scn,
    coef_setting = coef_setting,
    use_sl = use_sl,
    suffix = suffix_arg,
    cores = cores
  )
  
  #### Use the methods ####
  results <- foreach(sim.idx=1:num_simulations, sim=simulations,
                     .combine=rbind
  ) %dopar% {
    
    if(use_sl){
      mus <- mu.hats[[sim.idx]]
      dc <- sim$d - mus$mu.dx
      m_0 <- sim$m - mus$mu.mxis
      y_0 <- sim$y - mus$mu.yx
      
      mu_err <- cbind(
        muderr=(mus$mu.dx-sim$propensity)^2,
        mumerr=(mus$mu.mxis-sim$true_em)^2,
        muyerr=(mus$mu.yx-sim$true_ey)^2
      ) %>% colMeans
      
    } else {
      dc <- sim$d - sim$propensity
      m_0 <- sim$m - sim$true_em
      y_0 <- sim$y - sim$true_ey
      
      mu_err <- c(
        muderr=0,
        mumerr=rep(0, ncol(sim$m)),
        muyerr=0
      )
    }
    
    names(mu_err)[-c(1, length(mu_err))] <- paste0("mumerr", 1:p)
    this_opt <- tibble(opt, mu_err %>% as.list %>% as_tibble)
    
    print("Beginning prd")
    prd_results <- cv_sl_estimates_w_sel(
      y_0, m_0, dc, weight.version="product", weight.gam=weight_gam,
      lambdas=lambdas, folds=folds, opt=this_opt, 
      num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2349871
    ) %>% mutate(sim=sim.idx)
    
    print("Beginning mix")
    mix_results <- cv_sl_estimates_w_sel(
      y_0, m_0, dc, weight.version="mixture", weight.gam=weight_gam,
      lambdas=lambdas, folds=folds, opt=this_opt,
      num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2349871
    ) %>% mutate(sim=sim.idx)
    
    print("Beginning adp")
    adp_results <- cv_sl_estimates_w_sel(
      y_0, m_0, dc, weight.version="adaptive", weight.gam=weight_gam,
      lambdas=lambdas, folds=folds, opt=this_opt, 
      num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2349871
    ) %>% mutate(sim=sim.idx)
    
    print("Beginning full")
    full_results <- cv_sl_estimates_no_sel(
      y_0, m_0, dc, model_name="full", opt=this_opt, 
      num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2349871
    ) %>% mutate(sim=sim.idx)
    
    print("Beginning oracle")
    oracle_results <- cv_sl_estimates_no_sel(
      y_0, m_0[,true_M], dc, model_name="oracle", opt=this_opt, 
      num_bootstraps=num_bootstraps, do.boot=TRUE, boot.seed=2349871
    ) %>% mutate(sim=sim.idx)
    
    return(
      rbind(
        prd_results, mix_results, adp_results,
        full_results, oracle_results
      ) %>% mutate(
        coverage_NDE = (lower_NDE <= NDE) & (NDE <= upper_NDE),
        coverage_NIE = (lower_NIE <= NIE) & (NIE <= upper_NIE),
        err_NDE = NDE_hat - NDE,
        err_NIE = NIE_hat - NIE
        # coverage=if_else(
        #   target == "NDE",
        #   lower <= NDE & NDE <= upper,
        #   lower <= NIE & NIE <= upper
        # )
      )
    )
  }
  
  
  print("Scenario Complete!")
  stop = Sys.time()
  print(stop - start)
  
  
  return(results)
  
}
