# estimate-nuisance.R
use_backfitting <- F
if(use_sl == "T"){
  
  print("SL: Y")
  if(file.exists(slfile_y)){
    load(slfile_y)
  } else {
    set.seed(234987)
    mu.hats.y <- foreach(sim=simulations) %dopar% {
      
      train_sl_mu(
        sim, candidate.mediators=candidate.mediators,
        x.cols=x.cols,
        treat.col="d", outcome.col="y", 
        bin_lib=bin_lib, cont_lib=cont_lib, 
        folds=split_folds,
        estD=F, 
        estM=F,
        estY=T
      )
    }
    save(mu.hats.y, file=slfile_y)
  }
  
  print("SL: MD")
  if(file.exists(slfile_md)){
    load(slfile_md)
  } else {
    set.seed(234987)
    
    mu.hats.md <- foreach(sim=simulations) %dopar% {
      mu.hat.md <- train_sl_mu(
        sim, candidate.mediators=candidate.mediators,
        x.cols=x.cols,
        treat.col="d", outcome.col="y", 
        bin_lib=bin_lib, cont_lib=cont_lib, 
        folds=split_folds, 
        estD=T, 
        estM=T,
        estY=F
      )
      if(!use_backfitting){
        return(mu.hat.md)
      } else {
        
        set.seed(234987)
        
        # backfit to get m residuals
        backfit_res <- backfit_psi_m(sim, mu.hats=mu.hat.md,
                                     candidate.mediators=candidate.mediators,
                                     treat.col="d")
        m_1 <- backfit_res$m_1
        alpha_tildes <- backfit_res$alpha_tildes
        sim[,candidate.mediators] <- m_1
        
        # retrain for psi
        psi.m.hat <- train_sl_mu(
          sim, candidate.mediators=candidate.mediators,
          x.cols=x.cols,
          treat.col="d", outcome.col="y", 
          bin_lib=bin_lib, cont_lib=cont_lib, 
          folds=split_folds, 
          estD=F, 
          estM=T,
          estY=F
        )
        
        # use propensity to get mu.m
        mu.mxis <- foreach(j=1:p) %do% {
          res <- psi.m.hat$mu.mxis[[j]] + alpha_tildes[j]*mu.hat.md$mu.dx 
          stopifnot(is.numeric(res))
          return(res)
        }
        psi.m.hat$mu.mxis <- mu.mxis
        # psi.m.hat now contains the conditional mean
        
        mu.hat.md.prime <- c(psi.m.hat,
                             list(
                               mu.dx=mu.hat.md$mu.dx,
                               mu.dx.coef=mu.hat.md$mu.dx.coef
                             ))
        return(mu.hat.md.prime)
      }
    }
    
    save(mu.hats.md, file=slfile_md)
  }
  
  # Combine the MD and Y estimates in one
  mu.hats <- foreach(mu.y=mu.hats.y, mu.md=mu.hats.md) %do% c(mu.y, mu.md)
  
} else {
  mu.hats <- foreach(sim=simulations) %dopar% {
    true_sl_mu(sim, candidate.mediators=candidate.mediators)
  }
  
}