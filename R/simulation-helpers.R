# simulation-helpers.R
# Helper functions for simulations

true_sl_mu <- function(sel_df, candidate.mediators) {
  
  p = length(candidate.mediators)
  
  mu.dx = sel_df$propensity
  
  mu.mxis = lapply(1:p, function(i){
    mu.mxi <- sel_df$confounding_m[,i]
    return(mu.mxi)
  })
  
  mu.yx = sel_df$confounding_y
  
  return(
    list(mu.mxis=mu.mxis,
         mu.dx=mu.dx,
         mu.yx=mu.yx)
  )
  
}

my_sl_predict <- function(obj, Y=NULL){
  require(SuperLearner)
  obj_class <- class(obj)
  if(!obj_class %in% c("SuperLearner", "CV.SuperLearner")){
    #stop("Internal error: Only superlearners allowed for obj.")
    # treat these as already-computed predictions for each observation
    out <- obj
  }
  
  if(obj_class == "SuperLearner"){
    if(all(coef(obj) == 0))
    {
      if(is.null(Y)){
        stop("All metalearner coefficients are zero!")
      } else {
        out <- rep(mean(Y), length(Y))
      }
    } else {
      out <- predict(obj, onlySL=TRUE)$pred
    }
  } else if(obj_class == "CV.SuperLearner") {
    out <- obj$SL.predict
  }
  return(out)
}

result_dataframe <- function(results){
  true_M_str <- paste0("m.", true_M)
  foreach(this_res=results, .combine=rbind) %do% {
    data.frame(
      coverage_NDE = 1*((this_res$boot_ci[1,1] <= NDE) & (NDE <= this_res$boot_ci[2,1])),
      NDE_len = this_res$boot_ci[2,1] - this_res$boot_ci[1,1],
      coverage_NIE = 1*((this_res$boot_ci[1,2] <= NIE) & (NIE <= this_res$boot_ci[2,2])),
      NIE_len = this_res$boot_ci[2,2] - this_res$boot_ci[1,2],
      NIE_err = this_res$NIE_hat - NIE,
      NDE_err = this_res$NDE_hat - NDE,
      num_noise = length(setdiff(this_res$sel_M, true_M_str)),
      num_missed = length(setdiff(true_M_str, this_res$sel_M)),
      # which_correct = sapply(true_M, function(x) x %in% this_res$sel_M),
      picked_lam = this_res$lambda
    ) %>% 
      mutate(is_missed = 1*(num_missed > 0))
  }
}

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