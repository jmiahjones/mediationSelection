# estimate-nuisance.R
estimate_nuisance <- function(){
  
if(use_sl == T){
  
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
    # save(mu.hats.y, file=slfile_y)
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
      return(mu.hat.md)
    }
    
    # save(mu.hats.md, file=slfile_md)
  }
  
  # Combine the MD and Y estimates in one
  mu.hats <- foreach(mu.y=mu.hats.y, mu.md=mu.hats.md) %do% {c(mu.y, mu.md)}
  
} else {
  mu.hats <- foreach(sim=simulations) %dopar% {
    true_sl_mu(sim, candidate.mediators=candidate.mediators)
  }
  
}

return(mu.hats)
}