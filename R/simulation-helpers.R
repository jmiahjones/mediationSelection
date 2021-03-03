# simulation-helpers.R
# Helper functions for simulations

##### SL Helper ####
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

##### SL-Training #####
train_sl_mu <- function(sel_df, candidate.mediators, 
                        x.cols, treat.col, outcome.col,
                        bin_lib, cont_lib, folds, cores=NULL, 
                        parallel_outfile=NULL,
                        save_pred_only=TRUE,
                        estD=T, 
                        estM=T,
                        estY=T,
                        verbose=F) {
  .require(foreach)
  .require(quadprog)
  .require(SuperLearner)
  # p = length(candidate.mediators)
  V <- length(folds)
  SL.CV.control <- list(V=V, validRows=folds, shuffle=FALSE)
  
  X <- sel_df %>% dplyr::select(all_of(x.cols))
  
  # tests for D
  if(estD){
    D <- clean_D(sel_df, treat.col)
  }
  
  if(estM){
    M <- dplyr::select(sel_df, all_of(candidate.mediators))
  }
  
  if(estY){
    Y <- sel_df %>% dplyr::pull(outcome.col)
  }
  
  
  out <-  .train_sl_mu(X, 
                       if(estD) D else NULL, 
                       if(estM) M else NULL,
                       if(estY) Y else NULL,
                       SL.CV.control,
                       bin_lib, cont_lib, 
                       save_pred_only,
                       cores,
                       parallel_outfile, 
                       estD, estM, estY,
                       verbose
  )
  return(out)
}

.train_sl_mu <- function(X, D, M, Y,
                         SL.CV.control,
                         bin_lib, cont_lib, 
                         save_pred_only,
                         cores=NULL, 
                         parallel_outfile=NULL,
                         estD=T, 
                         estM=T,
                         estY=T,
                         verbose){
  
  p <- ncol(M)
  save_coefs <- !save_pred_only
  
  out <- list()
  if(estM){
    print("Estimating Mu.mx")
    
    parallelB = !is.null(cores) & is.numeric(cores)
    if(parallelB){
      .require(doParallel)
      .require(foreach)
      cl = makeCluster(cores, "FORK", outfile=parallel_outfile)
      registerDoParallel(cl)
    }
    
    mu.mxis = foreach(i=1:p, 
                      .combine=ifelse(save_coefs, foreach:::defcombine,
                                      cbind)
    ) %dopar% {
      
      if(verbose)
        print(paste0(Sys.time(), ": Starting ", i, " of ", p))
      Mi <- dplyr::pull(M, i)
      
      mu.mxi <- try(
        CV.SuperLearner(Y=Mi, X=X, SL.library=cont_lib, family=gaussian(),
                        cvControl=SL.CV.control,
                        saveAll=F,
                        method="method.NNLS")
      )
      
      if(inherits(mu.mxi, "try-error")){
        mu.mxi <- CV.SuperLearner(Y=Mi, X=X, SL.library=cont_lib, family=gaussian(),
                                  cvControl=SL.CV.control,
                                  saveAll=F,
                                  method="method.CC_LS")
        
        if(inherits(mu.mxi, "try-error")){
          warning(paste0("SuperLearner for M|X Failed: ", mi))
          mu.mxi <- NULL
        }
      }
      
      mu.mxi.pred = my_sl_predict(mu.mxi)
      
      if(save_coefs){
        mu.mxi.coef = mu.mxi$coef
        outm <- list(mu.mxi=mu.mxi.pred,
                     mu.mxi.coef=mu.mxi.coef)
      } else {
        outm <- mu.mxi.pred
      }
      
      return(outm)
    }
    
    if(save_coefs){
      mu.mxi.preds = foreach(x=mu.mxis, .combine=cbind) %do% {x$mu.mxi}
      mu.mxi.coefs = lapply(mu.mxis, function(x) x$mu.mxi.coef)
      out <- c(
        out, 
        list(mu.mxis=mu.mxi.preds,
             mu.mxi.coefs=mu.mxi.coefs
        )
      )
    } else {
      out <- c(
        out,
        list(mu.mxis=mu.mxis)
      )
    }
    
    rm(mu.mxis)
    
    if(parallelB){
      stopCluster(cl)
    }
    
  }
  
  
  if(estD){
    print("Estimating Mu.dx")
    
    mu.dx = try(
      CV.SuperLearner(Y=D, X=X, 
                      SL.library=bin_lib, family=binomial(),
                      cvControl=SL.CV.control,
                      saveAll=F,
                      method="method.NNloglik")
    )
    
    if(inherits(mu.dx, "try-error")){
      mu.dx <- CV.SuperLearner(Y=D, X=X, 
                               SL.library=bin_lib, family=binomial(),
                               cvControl=SL.CV.control,
                               saveAll=F,
                               method="method.CC_nloglik")
    }
    if(inherits(mu.dx, "try-error")){
      stop("SuperLearner for D|X Failed!")  
    }
    
    mu.dx.pred = my_sl_predict(mu.dx)
    if(save_coefs){
      mu.dx.coef = mu.dx$coef
      out <- c(
        out,
        list(mu.dx=mu.dx.pred,
             mu.dx.coef=mu.dx.coef
        )
      )
    } else {
      out <- c(
        out,
        list(mu.dx=mu.dx.pred)
      )
    }
    
    
    rm(mu.dx)
  }
  
  if(estY){
    
    print("Estimating Mu.yx")
    
    mu.yx <- try( 
      CV.SuperLearner(Y=Y, X=X,
                      SL.library=cont_lib, family=gaussian(),
                      cvControl=SL.CV.control,
                      saveAll=F,
                      method="method.NNLS")
    )
    
    if(inherits(mu.yx, "try-error")){
      mu.yx = CV.SuperLearner(Y=Y, X=X,
                              SL.library=cont_lib, family=gaussian(),
                              cvControl=SL.CV.control,
                              saveAll=F,
                              method="method.CC_LS")
    }
    mu.yx.pred = my_sl_predict(mu.yx)
    if(save_coefs){
      mu.yx.coef = mu.yx$coef
      out <- c(
        out,
        list(
          mu.yx=mu.yx.pred,
          mu.yx.coef=mu.yx.coef
        )
      )
    } else {
      out <- c(
        out,
        list(mu.yx=mu.yx.pred)
      )
    }
    
    rm(mu.yx)
    
  }
  
  return(out)
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


# DEPRECATED #####

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