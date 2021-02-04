#
# Author: Jeremiah Jones
# rm(list=objects())
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  stop("No commandline arguments found!")
}
# 
# # abbv_scn <- "nnl"
# 
stopifnot(length(args)==7)
n <- as.numeric(args[1])
num_simulations <- as.numeric(args[2])
abbv_scn <- args[3]
coef_setting <- args[4]
weight_gam <- as.numeric(args[5])
use_sl <- args[6]
suffix_arg <- args[7]


if(is.na(weight_gam)){
  # weight_gam was not numeric, so assume we want to cv it
  weight_gam_str <- gsub(" ", "", args[5])
  weight_gam <- c(1,2,3)
} else {
  weight_gam_str <- weight_gam
}

# stopifnot(basename(getwd()) == "pseudo-mse-with-ml")


# n <- 500
# num_simulations <- 5
# abbv_scn <- "lnn"
# suffix <- "ml"

####################
## parameters
# n = 2000
# num_simulations = 200
num_noise = 7

num_bootstraps=1000L
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
covariates_size = 5
V = 10 # number of folds for crossvalidation
splits = ifelse(n<2000, 4L, 2L) # number of sample splits for cross-fitting
####################

####################
#### Files ####
####################
savefile_suffix <- paste(n, num_simulations, abbv_scn, coef_setting,
                         weight_gam_str, use_sl, suffix_arg, sep="-")
print(savefile_suffix)
slfile_y_suffix <- paste(n, num_simulations, abbv_scn, coef_setting,
                         suffix_arg, sep="-")
slfile_md_suffix <- paste(n, num_simulations, substr(abbv_scn, 1,2), 
                          coef_setting,
                          suffix_arg, sep="-")

logfile <- paste0("./logs/", savefile_suffix, "-var-selection.stdout")
savefile <- paste0("./cache/", savefile_suffix, "-var-selection.RData")
slfile_y <- paste0("./cache/sl/", slfile_y_suffix, "-train-y.RData")
slfile_md <- paste0("./cache/sl/", slfile_md_suffix, "-train-md.RData")

if(!dir.exists("./logs")){
  dir.create("./logs")
}
if(!dir.exists("./cache")){
  dir.create("./cache")
}
if(file.exists(logfile)){
  file.remove(logfile)
}

start = Sys.time()

library(dplyr)
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


candidate.mediators = paste0("m.", 1:p)
true_M <- which(alphas*betas != 0)
oracle.mediators <- candidate.mediators[true_M]
x.cols = paste0("x.", 1:covariates_size)

expit <- plogis

########################
## Confounding Mechanisms
########################

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
    expit((0*x1 + 1*x2 + 0*x3 + 0*x4 + 0*x5+0*x1^2*x2+0*x1*x3 - 0)/1)
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
    1*x1^2 + x2 # 0*x3 - 0*x4 - 0*x5 + 0*x2*(x1-0.5)^2
  }
} else if(substr(abbv_scn, 2,2)=="l"){
  # linear function
  psi_m <- function(x) { 
    x1 <- x[1]
    x2 <- x[2]
    x3 <- x[3]
    x4 <- x[4]
    x5 <- x[5]
    1*x1 + x2
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
    1*(2*(x1-0.5)^2 + x2 + 2*x3 - x4 - x5)
  }
} else if(substr(abbv_scn, 3,3)=="l"){
  # linear function
  psi_y <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    x3 <- x[3]
    x4 <- x[4]
    x5 <- x[5]
    (2*(x1-0.5) + x2 + 2*x3 - x4 - x5)
  }
}









set.seed(2018)


########################
## Load Procedures
########################
# source("./R/estimation-scripts-semipar-pseudomse-v2.R")
source("./R/simulation-helpers.R")
source("./R/mediation-funs-postsel.R")

########################
## SuperLearners
########################

# create_gam = create.Learner("SL.gam",
#                             tune=list(
#                               deg.gam=2
#                             ))
create_mgcv <- create.Learner("SL.mgcv", tune=list(
  m=2,
  k=c(2,3,5),
  bs=c("ts")
))

cont_lib = c(
  "SL.glm",
  "SL.glmnet",
  # "SL.ridge",
  # "SL.xgboost",
  # "SL.bartMachine",
  # "SL.step.interaction",
  # "SL.ksvm",
  "SL.earth",
  "SL.randomForest",
  create_mgcv$names,
  "SL.mean"
)

bin_lib = c(
  "SL.glm",
  "SL.glmnet",
  # "SL.xgboost",
  # "SL.bartMachine",
  # "SL.step.interaction",
  # "SL.ksvm",
  "SL.earth",
  "SL.randomForest",
  create_mgcv$names,
  "SL.mean"
)

########################
## Simulations
########################

# Setup a lambda grid
lambda.len = 301
lambdas = n^(-0.75)*(2^seq(from=-2, to=10, length.out=lambda.len))

set.seed(841664)
folds <- caret::createFolds(y=1:n, k=V)
set.seed(841665)
split_folds <- caret::createFolds(y=1:n, k=splits)

cl <- parallel::makeCluster(24, "FORK", outfile=logfile)
doParallel::registerDoParallel(cl)
parallel::clusterSetRNGStream(cl, 2018)


tryCatch({
  simulations <- foreach(sim.idx=1:num_simulations) %dopar% {
    
    ########################
    ## Simulation mechanism
    ########################
    
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
    
    # mediators
    m <- (as.matrix(d) %*% alphas) + confounding_m + epsilon_m
    # m[,p] <- m[,-p] %>% rowMeans
    # m[,p] <- m[,1]^2 + m[,2]^2 + rowMeans(m[,4:(p-1)])
    # m <- m + epsilon_m
    
    # outcome confounding
    confounding_y <- apply(x, 1, function(row){
      psi_y(row)
    })
    
    set.seed(7243957 + sim.idx)
    # outcome equation
    y <- (NDE * d) + (m %*% betas) + (confounding_y) + rnorm(n, mean=0, sd=2)
    
    # for troubleshooting only
    true_em <- (as.matrix(d.prob) %*% alphas) + confounding_m
    true_ey <- as.numeric(NDE*d.prob + true_em %*% as.matrix(betas) + confounding_y)
    
    # if(sim.idx %% 10 == 0)
    #   print(paste("Completed simulation", sim.idx))
    
    # save the true centering variables for future use
    # em = replicate(p,confounding_m) + replicate(p,propensity)*alphas
    # ey = NDE*propensity + (em %*% betas) + confounding_y
    
    return(
      tibble(
        data.frame(
          x=x, d=d, m=m, y=y,
          propensity=d.prob
        ),
        confounding_m = true_em,
        confounding_y = true_ey
      )
    )
  }
  
  print("Created data.")
  parallel::clusterExport(cl, "simulations")
  
  # source("./R/simulation-helpers.R")
  source("./R/estimate-nuisance.R")
  
  
  print("Estimated mu-hats with SuperLearner.")
  parallel::clusterExport(cl, "mu.hats")
  
  inference_switch <- function(idx, ..., split.folds, boot.sl) {
    if(idx == 1){
      cv_sl_estimates_w_boot(...)
    } else if(idx==2) {
      cv_sl_estimates_w_mult(..., split.folds=split.folds, boot.sl=boot.sl)
    } else {
      stop("Error: Invalid inference procedure!")
    }
  }
  parallel::clusterExport(cl, "inference_switch")
  
  inference.results <- foreach(fun.idx=1:2) %do% {
    mix.results <- foreach(sim=simulations, mus=mu.hats) %dopar% {
      inference_switch(fun.idx, sim, mu.hats=mus,
                       weight.version="mixture", weight.gam=weight_gam,
                       lambdas=lambdas, folds=folds,
                       candidate.mediators=candidate.mediators,
                       x.cols=x.cols,
                       treat.col="d", outcome.col="y",
                       num_bootstraps=num_bootstraps,
                       boot.seed=2349871,
                       split.folds=split_folds,
                       boot.sl=FALSE)
    }
    print("Finished lasso simulations: mixture.")
    
    prd.results <- foreach(sim=simulations, mus=mu.hats) %dopar% {
      inference_switch(fun.idx, sim, mu.hats=mus,
                       weight.version="product", weight.gam=weight_gam,
                       lambdas=lambdas, folds=folds,
                       candidate.mediators=candidate.mediators,
                       x.cols=x.cols,
                       treat.col="d", outcome.col="y",
                       num_bootstraps=num_bootstraps,
                       boot.seed=2349871,
                       split.folds=split_folds,
                       boot.sl=FALSE)
    }
    print("Finished lasso simulations: product.")
    
    adp.results <- foreach(sim=simulations, mus=mu.hats) %dopar% {
      inference_switch(fun.idx, sim, mu.hats=mus,
                       weight.version="adaptive", weight.gam=weight_gam,
                       lambdas=lambdas, folds=folds,
                       candidate.mediators=candidate.mediators,
                       x.cols=x.cols,
                       treat.col="d", outcome.col="y",
                       num_bootstraps=num_bootstraps,
                       boot.seed=2349871,
                       split.folds=split_folds,
                       boot.sl=FALSE)
    }
    print("Finished lasso simulations: adaptive.")
    
    print(paste0("Finished inference technique: ", fun.idx))
    
    list(mix.results=mix.results, prd.results=prd.results, adp.results=adp.results)
  }
  
  oracle.results <- foreach(sim=simulations, mus=mu.hats) %dopar% {
    cv_sl_estimates_no_sel(sim, mu.hats=mus,
                           candidate.mediators=oracle.mediators,
                           x.cols=x.cols,
                           treat.col="d", outcome.col="y",
                           num_bootstraps=num_bootstraps,
                           boot.seed=2349871)
  }
  print("Finished oracle results.")
  
  full.results <- foreach(sim=simulations, mus=mu.hats) %dopar% {
    cv_sl_estimates_no_sel(sim, mu.hats=mus,
                           candidate.mediators=candidate.mediators,
                           x.cols=x.cols,
                           treat.col="d", outcome.col="y",
                           num_bootstraps=num_bootstraps,
                           boot.seed=2349871)
  }
  print("Finished full results.")
  
},finally={
  # stopCluster(cl)
  save.image(file=savefile)
})


print("Scenario Complete!")
stop = Sys.time()
print(stop - start)



#######################
## Print Results
#######################
prd.df <- result_dataframe(inference.results[[1]]$prd.results)
prd.df %>% summary
mean(prd.df$num_missed > 0)
prd.df <- result_dataframe(inference.results[[2]]$prd.results)
prd.df %>% summary
mean(prd.df$num_missed > 0)

mix.df <- result_dataframe(inference.results[[1]]$mix.results)
mix.df %>% summary
mean(mix.df$num_missed > 0)
mix.df <- result_dataframe(inference.results[[2]]$mix.results)
mix.df %>% summary
mean(mix.df$num_missed > 0)


adp.df <- result_dataframe(inference.results[[1]]$adp.results)
adp.df %>% summary
mean(adp.df$num_missed > 0)
adp.df <- result_dataframe(inference.results[[2]]$adp.results)
adp.df %>% summary
mean(adp.df$num_missed > 0)

# mix.df <- result_dataframe(mix.results)
# prd.df <- result_dataframe(prd.results)
# adp.df <- result_dataframe(adp.results)
