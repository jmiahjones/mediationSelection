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
train_sl_mu <- function(df, m.cols, 
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
  
  V <- length(folds)
  SL.CV.control <- list(V=V, validRows=folds, shuffle=FALSE)
  
  X <- df %>% dplyr::pull(x.cols) %>% as.data.frame
  
  # tests for D
  if(estD){
    D <- clean_D(df, treat.col)
  }
  
  if(estM){
    M <- df %>% dplyr::pull(m.cols)
  }
  
  if(estY){
    Y <- df %>% dplyr::pull(outcome.col)
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
    } else {
      registerDoSEQ()
    }
    
    mu.mxis = foreach(i=1:p, 
                      .combine=ifelse(save_coefs, foreach:::defcombine,
                                      cbind)
    ) %dopar% {
      
      if(verbose)
        print(paste0(Sys.time(), ": Starting ", i, " of ", p))
      Mi <- M[,i]
      
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

