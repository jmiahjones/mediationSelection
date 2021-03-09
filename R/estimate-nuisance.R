# estimate-nuisance.R
library(nnls)
library(quadprog)
library(nloptr)

set.seed(841665)
split_folds <- caret::createFolds(y=1:n, k=splits)

# estimate_nuisance <- function(){
  
  
  ######################## SuperLearners ########################
  
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
    "SL.lm",
    # "SL.glmnet",
    # "SL.earth",
    # "SL.randomForest",
    # create_mgcv$names,
    "SL.mean"
    # "SL.ridge",
    # "SL.xgboost",
    # "SL.bartMachine",
    # "SL.step.interaction",
    # "SL.ksvm",
  )
  
  bin_lib = c(
    "SL.glm",
    # "SL.glmnet",
    # "SL.earth",
    # "SL.randomForest",
    # create_mgcv$names,
    "SL.mean"
    # "SL.xgboost",
    # "SL.bartMachine",
    # "SL.step.interaction",
    # "SL.ksvm",
  )
  
  
  
  print("Beginning SL...")
  set.seed(234987)
  
  # exports <- c(
  #   "m.cols",
  #   "x.cols",
  #   "bin_lib", 
  #   "cont_lib", 
  #   "split_folds"
  # )
  mu.hats <- foreach(sim=simulations, sim.idx=seq_along(simulations)
                     # .export=exports
  ) %dopar% {
    if(sim.idx %% 10 == 0) print(paste0("Estimating ", sim.idx, "..."))
    train_sl_mu(
      sim, m.cols=m.cols,
      x.cols=x.cols,
      treat.col="d", outcome.col="y", 
      bin_lib=bin_lib, cont_lib=cont_lib, 
      folds=split_folds,
      estD=T, 
      estM=T,
      estY=T
    )
  }
  
  # print("SL: MD")
  # set.seed(234987)
  
  # mu.hats.md <- foreach(sim=simulations) %dopar% {
  #   mu.hat.md <- train_sl_mu(
  #     sim, m.cols=m.cols,
  #     x.cols=x.cols,
  #     treat.col="d", outcome.col="y", 
  #     bin_lib=bin_lib, cont_lib=cont_lib, 
  #     folds=split_folds, 
  #     estD=T, 
  #     estM=T,
  #     estY=F
  #   )
  #   return(mu.hat.md)
  # }
  
  
  # Combine the MD and Y estimates in one
  # mu.hats <- foreach(mu.y=mu.hats.y, mu.md=mu.hats.md) %do% {c(mu.y, mu.md)}
  # 
  # rm(mu.hats.y, mu.hats.md)
  # return(mu.hats)
# }
