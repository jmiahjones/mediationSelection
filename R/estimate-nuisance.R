# estimate-nuisance.R
library(nnls)
library(quadprog)
library(nloptr)

splits = ifelse(n<2000, 4L, 2L) # number of sample splits for cross-fitting
set.seed(841665)
split_folds <- caret::createFolds(y=1:n, k=splits)

#### My SuperLearners ####
SL.mgcv <- function (Y, X, newX, family, obsWeights, 
                     bs="tp", m=3, k=6, select=FALSE,
                     cts.num = 4, 
                     ...) 
{
  .require("mgcv")
  stopifnot(packageVersion("mgcv") >= 1.8)
  if ("gam" %in% loadedNamespaces()) 
    warning("mgcv and gam packages are both in use. You might see an error because both packages use the same function names.")
  cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(
      paste("Y~", 
            paste(
              paste0("s(", 
                     colnames(X[, cts.x, drop = FALSE]), 
                     ", bs='", bs, "', m=", m, ", k=", k,
                     ")"
              ), collapse = "+"
            ), 
            "+", 
            paste(
              colnames(X[,!cts.x, drop = FALSE]), 
              collapse = "+"
            )
      )
    )
  }
  else {
    gam.model <- as.formula(
      paste(
        "Y~", 
        paste(
          paste0(
            "s(", colnames(X[, cts.x, drop = FALSE]), 
            ", bs='", bs, "', m=", m, ", k=", k,
            ")"
          ), 
          collapse = "+"
        )
      )
    )
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(
      paste("Y~", 
            paste(colnames(X), collapse = "+")
      )
    )
  }
  fit.gam <- mgcv::gam(gam.model, data = X, family = family, 
                       select=select,
                       control = mgcv::gam.control(maxit = 50), 
                       weights = obsWeights)
  pred <- mgcv::predict.gam(fit.gam, newdata = newX, type = "response")
  
  fit <- list(object = fit.gam)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.mgcv")
  return(out)
}

predict.SL.mgcv <- function (object, newdata, ...) 
{
  .require("mgcv")
  stopifnot(packageVersion("mgcv") >= 1.8)
  pred <- mgcv::predict.gam(object = object$object, newdata = newdata, 
                            type = "response")
  return(pred)
}

my.SL.randomForest <- function(Y, X, newX, family, 
                               mtry = ifelse(family$family == "gaussian", 
                                             max(floor(ncol(X)/3), 1), 
                                             floor(sqrt(ncol(X)))), 
                               ntree = 1000,
                               nodesize = ifelse(family$family == "gaussian", 5, 1), 
                               maxnodes = NULL, importance = FALSE, ...) 
{
  require("randomForest")
  if (family$family == "gaussian") {
    fit.rf <- randomForest::randomForest(y=Y, x = X, 
                                         ntree = ntree, xtest = newX, keep.forest = TRUE, 
                                         mtry = mtry, nodesize = nodesize, maxnodes = maxnodes, 
                                         importance = importance)
    pred <- fit.rf$test$predicted
    fit <- list(object = fit.rf)
  }
  if (family$family == "binomial") {
    fit.rf <- randomForest::randomForest(y = as.factor(Y), 
                                         x = X, ntree = ntree, xtest = newX, keep.forest = TRUE, 
                                         mtry = mtry, nodesize = nodesize, maxnodes = maxnodes, 
                                         importance = importance)
    pred <- fit.rf$test$votes[, 2]
    fit <- list(object = fit.rf)
  }
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.randomForest")
  return(out)
}

predict.my.SL.randomForest <- function(object, newdata, family, ...) 
{
  require("randomForest")
  if (family$family == "gaussian") {
    pred <- predict(object$object, newdata = newdata, type = "response")
  }
  if (family$family == "binomial") {
    pred <- predict(object$object, newdata = newdata, type = "vote")[, 
                                                                     2]
  }
  pred
}


# SL.glmnet
# function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, 
#           nlambda = 100, useMin = TRUE, loss = "deviance", ...) 
# {
#   .SL.require("glmnet")
#   if (!is.matrix(X)) {
#     X <- model.matrix(~-1 + polym(., degree=degree), X)
#     newX <- model.matrix(~-1 + ., newX)
#   }
#   fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights, 
#                              lambda = NULL, type.measure = loss, nfolds = nfolds, 
#                              family = family$family, alpha = alpha, nlambda = nlambda, 
#                              ...)
#   pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin, 
#                                                                     "lambda.min", "lambda.1se"))
#   fit <- list(object = fitCV, useMin = useMin)
#   class(fit) <- "SL.glmnet"
#   out <- list(pred = pred, fit = fit)
#   return(out)
# }



# estimate_nuisance <- function(){
  
  
  ######################## SuperLearner Library ########################
  
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
    "SL.glmnet",
    "SL.earth",
    "SL.randomForest",
    create_mgcv$names
    # "SL.mean"
    # "SL.ridge",
    # "SL.xgboost",
    # "SL.bartMachine",
    # "SL.step.interaction",
    # "SL.ksvm",
  )
  
  bin_lib = c(
    "SL.glm",
    "SL.glmnet",
    "SL.earth",
    "SL.randomForest",
    create_mgcv$names
    # "SL.mean"
    # "SL.xgboost",
    # "SL.bartMachine",
    # "SL.step.interaction",
    # "SL.ksvm",
  )
  
  
  ######################## SuperLearning ########################
  print("Beginning SL...")
  set.seed(234987)
  
  exports <- c(
    "SL.mgcv",
    bin_lib,
    cont_lib
  )
  mu.hats <- foreach(sim=simulations, sim.idx=seq_along(simulations),
                     .export=exports
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
