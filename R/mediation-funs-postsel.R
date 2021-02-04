# mediation-funs-alasso-boot.R
# Functions that perform the estimation of NDE and NIE using semiparametric confounding control.

clean_D <- function(sel_df, treat.col){
  D <- dplyr::pull(sel_df, treat.col)
  if(is.factor(D)){
    if(length(levels(D)) != 2){
      stop("Only accepting binary treatments!")
    }
    message("Assuming exposure is the second level.")
    active_lvl = levels(D)[2]
    D = 1*(D == active_lvl)
    
  } else {
    lvls = sort(unique(D))
    if(length(lvls) != 2) {
      stop("Only accepting binary treatments! Treatment does not have 2 unique values.")
    }
    if(!all.equal(lvls, c(0,1))){
      message("Assuming exposure is the largest value.")
      active_lvl = lvls[2]
      D = 1*(D == active_lvl)
    }
  }
  return(D)
}

robinsonize_system <- function(D, M, Y, mu.hats){
  .require(tibble)
  
  p <- ncol(M)
  mu.mxis <- mu.hats$mu.mxis
  mu.dx <- mu.hats$mu.dx
  
  d.prob = as.numeric(my_sl_predict(mu.dx))
  dc = D - d.prob
  
  m_0 = sapply(1:p, function(i){
    # estimate the counfounding effect on mediatior i
    
    Mi <- dplyr::pull(M, i)
    mu.mxi = mu.mxis[[i]]
    
    col <- Mi - my_sl_predict(mu.mxi)
    return(col)
  })
  stopifnot(ncol(m_0)==p)
  
  rob_tbl <- tibble::tibble(
    dc=dc,
    m_0=m_0
  )
  
  if(!is.null(Y)){
    mu.yx <- mu.hats$mu.yx
    y_0 = as.numeric(Y - my_sl_predict(mu.yx))
    rob_tbl <- rob_tbl %>% tibble::add_column(y_0=y_0)
  }
  
  return(
    rob_tbl
  )
}

.caretFoldToGlmnetFold <- function(folds, n){
  V <- length(folds)
  foldid <- rep(0, n)
  for(fold.idx in 1:V){
    holdout.idx <- folds[[fold.idx]]
    foldid[holdout.idx] <- fold.idx
  }
  if(any(foldid == 0))
    stop("Error in folds!")
  return(foldid)
}

.require <- function(...)
  require(..., quietly=TRUE)

#################################################################
### Methods for computing cross-validation
#################################################################

##### Cross-validated SuperLearner #####
# TODO finish re-factoring to make everything modular

# shorter version
# TODO add in Das bootstrap
cv_sl_estimates_w_boot <- function(sel_df, mu.hats,
                                   weight.version, weight.gam=1,
                                   lambdas, folds,
                                   candidate.mediators, x.cols, 
                                   treat.col, outcome.col,
                                   do.boot=TRUE,
                                   num_bootstraps=1000,
                                   boot.seed=NA) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  .require(dplyr)
  .require(SuperLearner)
  .require(foreach)
  
  n <- nrow(sel_df)
  
  # selection data
  X <- sel_df %>% dplyr::select(all_of(x.cols))
  D <- clean_D(sel_df, treat.col)
  M <- sel_df %>% dplyr::select(all_of(candidate.mediators))
  Y <- sel_df %>% dplyr::pull(outcome.col)
  
  p <- ncol(M)
  
  rob_system <- robinsonize_system(D, M, Y, mu.hats)
  dc <- rob_system$dc
  m_0 <- rob_system$m_0
  y_0 <- rob_system$y_0
  
  mfit <- lm(m_0 ~ dc)
  alpha_tildes = coef(mfit)[2, ]
  summ_mfit <- summary(mfit)
  if(do.boot){
    mhat <- predict(mfit)
    etahat <- mfit$residuals
    Sigma.m.breve.n <- vcov(mfit)[2*(1:p),2*(1:p)]
  }
  rm(mfit)
  
  # message(class(y_0))
  yfit = lm(y_0 ~ dc + m_0)
  beta_tildes = coef(yfit)[-(1:2)]
  sigma.check.sq.n <- mean(yfit$residuals^2)
  
  
  design_mat <- cbind(dc, m_0)
  foldid <- .caretFoldToGlmnetFold(folds, n)
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights.0 <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights.0 <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights.0 <- 1/(abs(beta_tildes))
  )
  
  weight.gam.len <- length(weight.gam)
  best.adapt.fit <- NULL
  best.kap <- NULL
  min.risk <- Inf
  for(weight.gam.idx in 1:weight.gam.len){
    kap <- weight.gam[weight.gam.idx]
    weights <- weights.0^(kap)
    adapt.fit <- glmnet::cv.glmnet(x=design_mat, y=y_0, foldid=foldid, 
                                   intercept=FALSE, 
                                   family = "gaussian", 
                                   penalty.factor = c(0, weights), 
                                   lambda = lambdas
    )
    this.risk <- min(adapt.fit$cvm)
    if(this.risk < min.risk){
      min.risk <- this.risk
      best.adapt.fit <- adapt.fit
      best.kap <- kap
    }
  }
  
  adapt.fit <- best.adapt.fit
  min.loss.lambda <- adapt.fit$lambda.min
  
  adapt.fit.coef <- coef(adapt.fit, s="lambda.min")
  is_selected <- adapt.fit.coef[-(1:2),] != 0
  selected_idx <- which(is_selected)
  alpha_hats = alpha_tildes[selected_idx]
  gamma_hat <- adapt.fit.coef[2,]
  beta_hats <- adapt.fit.coef[-(1:2),][selected_idx]
  NIE_hat <- sum(alpha_hats * beta_hats)
  
  return_list <- list(
    alpha_hats = alpha_hats, beta_hats = beta_hats,
    alpha_tildes=alpha_tildes, beta_tildes=beta_tildes,
    NIE_hat = NIE_hat, NDE_hat = gamma_hat,
    ATE_hat = NIE_hat + gamma_hat,
    sel_M = candidate.mediators[is_selected],
    lambda=min.loss.lambda,
    kap=best.kap
    )
  
  if(do.boot){
    # TODO Do Das Bootstrap
    yhat <- predict(adapt.fit, newx=design_mat, s="lambda.min")
    epshat <- y_0 - yhat
    n <- length(y_0)
    sigma.hat <- sqrt(mean(epshat^2))
    # Cn11.inv <- n*MASS::ginv(crossprod(design_mat[, c(T,is_selected)]))
    # Cn.inv <- n*MASS::ginv(crossprod(design_mat))
    # 
    # # TODO allow general gamma
    # eta_mat <- min.loss.lambda/(2*n)*
    #   # Cn11.inv %*%
    #   sweep(design_mat[, c(T,is_selected)],
    #         2,
    #         sign(beta_hats)/(abs(beta_hats)^2),
    #         "*")

    # Sigma.hat.n <- 1/n*crossprod(
    #   design_mat[, c(T,is_selected)] %*% Cn11.inv
    # )
    # Sigma.hat.n.half <- chol(Sigma.hat.n)
    # Sigma.hat.n.half.inv <- MASS::ginv(Sigma.hat.n.half)
    
    # Sigma.check.n <- 1/n*crossprod(
    #   design_mat[, c(T,is_selected)] %*% Cn11.inv
    # )
    # Sigma.check.n.half <- chol(Sigma.check.n)
    # Sigma.check.n.half.inv <- MASS::ginv(Sigma.check.n.half)
    # 
    # Sigma.tilde.n <- 1/n*crossprod(
    #   epshat*(design_mat[, c(T,is_selected)]) %*% Cn11.inv
    # )
    # Sigma.tilde.n.half <- chol(Sigma.tilde.n)
    # Sigma.tilde.n.half.inv <- MASS::ginv(Sigma.tilde.n.half)
    # 
    # # TODO: allow other gamma
    # b.breve.n <- Cn11.inv %*% ((min.loss.lambda/sqrt(n))*
    #                              sign(beta_hats)/abs(beta_tildes)) 
    # 
    # R.n <- (1/sqrt(sigma.check.sq.n))*Sigma.check.n.half.inv %*% (beta_hats + b.breve.n)
    
    # Sigma.breve.n <- 1/n*crossprod(design_mat[, c(T,is_selected)] + eta_mat)
    # Sigma.breve.n.half <- chol(Sigma.breve.n)
    # Sigma.breve.n.half.inv <- MASS::ginv(Sigma.breve.n.half)
    # 
    # # TODO: Figure out what kind of studentization to do
    # sigma.d.hat <- sqrt(mean(dc))
    # Sigma.m.breve.n.half <- chol(Sigma.m.breve.n)
    # Sigma.m.breve.n.half.inv <- MASS:ginv(Sigma.m.breve.n.half)
    
    
    B=num_bootstraps
    boot_results <- foreach(boot.idx=1:B, .combine=rbind) %do% {
      if(is.numeric(boot.seed))
        set.seed(boot.seed + boot.idx)
      G <- rbeta(n, shape1=1/2, shape2=3/2)
      muG <- .25 #.5/(.5 + 1.5)
      Gtil <- (G - muG)/muG
      
      H <- -1 + 2*rbinom(n, size=1, prob=0.5) # Rademacher rv
      etastar <- sweep(etahat, 1, H, "*")
      epsstar <- epshat*Gtil
      
      m_0_star <- mhat + etastar
      y_0_star <- yhat + epsstar
      
      mfit_star <- lm(m_0_star ~ dc)
      alpha_tildes_star = coef(mfit_star)[2, ]
      
      yfit = lm(y_0_star ~ design_mat - 1)
      beta_tildes_star = coef(yfit)[-1]
      sigma.check.sq.n.star <- mean((yfit$residuals*Gtil)^2)
      
      switch (weight.version,
              "mixture" = weights <- 1/(abs(beta_tildes_star) * (1 + abs(alpha_tildes_star))),
              "product" = weights <- 1/(abs(beta_tildes_star) * abs(alpha_tildes_star)),
              "adaptive"= weights <- 1/(abs(beta_tildes_star))
      )
      weights <- weights^(best.kap)
      
      adapt.fit.star <- glmnet::glmnet(x=design_mat, y=y_0_star,
                                       intercept=FALSE, 
                                       family = "gaussian", 
                                       penalty.factor = c(0, weights), 
                                       lambda = lambdas
      )
      # min.loss.lambda.star <- adapt.fit.star$lambda.min
      adapt.fit.star.coef <- coef(adapt.fit.star, s=min.loss.lambda)
      is_selected_star <- adapt.fit.star.coef[-(1:2),] != 0
      selected_idx_star <- which(is_selected_star)
      alpha_hats_star = alpha_tildes_star[selected_idx_star]
      gamma_hat_star <- adapt.fit.star.coef[2,]
      beta_hats_star <- adapt.fit.star.coef[-(1:2),][selected_idx_star]
      NIE_hat_star <- sum(alpha_hats_star * beta_hats_star)
      
      # TODO: Enable studentized pivots
      # y_hat_star <- predict(adapt.fit.star, newx=design_mat, s="lambda.min")
      # eps_hat_star <- y_0_star - y_hat_star
      # sigma_hat_star <- sqrt(mean((eps_hat_star*Gtil)^2))
      #
      # sigma_m_hat_star <- sqrt(mean((mfit_star$residuals*Gtil)^2))
      # 
      # Sigma.m.breve.star <- vcov(mfit_star)[2*(1:p),2*(1:p)]
      # 
      # Tstar_alpha <- alpha_hats_star - alpha_hats
      # Tstar_beta <- beta_hats_star - beta_hats
      # Tstar_gamma <- gamma_hat_star - gamma_hat
      # 
      # pivot_outcome_model <- sigma.hat/sigma_hat_star*Sigma.breve.n.half.inv %*%
      #   c(Tstar_gamma, Tstar_beta)
      # 
      # pivot_mediators <- sigma.m.hat/sigma_m_hat_star*Sigma.m.breve.n.half.inv %*%
      #   Tstar_alpha
      
      # boot_return <- list(
      #   alpha_star = alpha_hats_star,
      #   beta_star = beta_hats_star,
      #   gamma_star = gamma_hat_star,
      #   NIE_star = NIE_hat_star,
      #   is_selected_star = is_selected_star
      # )
      boot_return <- c(
        gamma_star = gamma_hat_star,
        NIE_star = NIE_hat_star
      )
      return(boot_return)
    }
    
    boot_ci <- apply(boot_results, 2, quantile, probs=c(0.025, 0.975))
    
    return_list <- c(
      return_list, 
      list(boot_ci=boot_ci)
    )
  }
  
  return(
    return_list
  )
  
}


#####################################################
# Minnier Multiplier Bootstrap
#####################################################
cv_sl_estimates_w_mult <- function(sel_df, mu.hats,
                                   weight.version, weight.gam=1,
                                   lambdas, folds,
                                   candidate.mediators, x.cols, 
                                   treat.col, outcome.col,
                                   num_bootstraps=1000L,
                                   do.boot=TRUE,
                                   boot.seed=NA,
                                   boot.sl=TRUE,
                                   split.folds=NULL) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  .require(dplyr)
  .require(SuperLearner)
  .require(foreach)
  
  if(boot.sl & is.null(split.folds))
    stop("Need split.folds for boostrapping SL!")
  
  n <- nrow(sel_df)
  
  # selection data
  X <- sel_df %>% dplyr::select(all_of(x.cols))
  D <- clean_D(sel_df, treat.col)
  M <- sel_df %>% dplyr::select(all_of(candidate.mediators))
  Y <- sel_df %>% dplyr::pull(outcome.col)
  
  p <- ncol(M)
  
  rob_system <- robinsonize_system(D, M, Y, mu.hats)
  dc <- rob_system$dc
  m_0 <- rob_system$m_0
  y_0 <- rob_system$y_0
  
  mfit <- lm(m_0 ~ dc)
  alpha_tildes = coef(mfit)[2, ]
  summ_mfit <- summary(mfit)
  rm(mfit)
  
  yfit = lm(y_0 ~ dc + m_0)
  beta_tildes = coef(yfit)[-(1:2)]
  
  
  design_mat <- cbind(dc, m_0)
  foldid <- .caretFoldToGlmnetFold(folds, n)
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights.0 <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights.0 <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights.0 <- 1/(abs(beta_tildes))
  )
  
  weight.gam.len <- length(weight.gam)
  best.adapt.fit <- NULL
  best.kap <- NULL
  min.risk <- Inf
  for(weight.gam.idx in 1:weight.gam.len){
    kap <- weight.gam[weight.gam.idx]
    weights <- weights.0^(kap)
    adapt.fit <- glmnet::cv.glmnet(x=design_mat, y=y_0, foldid=foldid, 
                                   intercept=FALSE, 
                                   family = "gaussian", 
                                   penalty.factor = c(0, weights), 
                                   lambda = lambdas
    )
    this.risk <- min(adapt.fit$cvm)
    if(this.risk < min.risk){
      min.risk <- this.risk
      best.adapt.fit <- adapt.fit
      best.kap <- kap
    }
  }
  
  adapt.fit <- best.adapt.fit
  min.loss.lambda <- adapt.fit$lambda.min
  
  adapt.fit.coef <- coef(adapt.fit, s="lambda.min")
  is_selected <- adapt.fit.coef[-(1:2),] != 0
  selected_idx <- which(is_selected)
  alpha_hats = alpha_tildes[selected_idx]
  gamma_hat <- adapt.fit.coef[2,]
  beta_hats <- adapt.fit.coef[-(1:2),][selected_idx]
  NIE_hat <- sum(alpha_hats * beta_hats)
  
  return_list <- list(
    alpha_hats = alpha_hats, beta_hats = beta_hats,
    alpha_tildes=alpha_tildes, beta_tildes=beta_tildes,
    NIE_hat = NIE_hat, NDE_hat = gamma_hat,
    ATE_hat = NIE_hat + gamma_hat,
    sel_M = candidate.mediators[is_selected],
    lambda=min.loss.lambda,
    kap=best.kap
  )
  
  if(do.boot){
    B=num_bootstraps
    boot_results <- foreach(boot_idx=1:B, .combine=rbind) %do% {
      if(is.numeric(boot.seed))
        set.seed(boot.seed + boot_idx)
      # G <- rbeta(n, shape1=1/2, shape2=3/2)
      # muG <- .25 #.5/(.5 + 1.5)
      # Gtil <- (G - muG)/muG
      # G <- Gtil
      
      if(boot.sl){
        H <- rexp(n=n, rate=1) # multiplier for SL recombination
        boot_sys <- multiplier_bootstrap_recombine_SL(
          split.folds, mu.hats,
          Y=Y, M=M, D=D, obsWeights=H
        )
        y_0_star <- boot_sys$y_0
        m_0_star <- boot_sys$m_0
        dc_star <- boot_sys$dc
        design_mat_star <- cbind(dc_star, m_0_star)
      } else {
        y_0_star <- y_0
        m_0_star <- m_0
        dc_star <- dc
        design_mat_star <- design_mat
      }
      G <- rgamma(n=n, shape=1, rate=1)
      
      mfit <- lm(m_0_star ~ dc_star, weights = G)
      alpha_tildes_star = coef(mfit)[2, ]
      
      yfit = lm(y_0_star ~ design_mat_star - 1, weights = G)
      beta_tildes_star = coef(yfit)[-1]
      
      switch (weight.version,
              "mixture" = weights <- 1/(abs(beta_tildes_star) * (1 + abs(alpha_tildes_star))),
              "product" = weights <- 1/(abs(beta_tildes_star) * abs(alpha_tildes_star)),
              "adaptive"= weights <- 1/(abs(beta_tildes_star))
      )
      weights <- weights^(best.kap)
      
      adapt.fit.star <- glmnet::glmnet(x=design_mat_star, y=y_0_star,
                                       weights = G,
                                       intercept=FALSE, 
                                       family = "gaussian", 
                                       penalty.factor = c(0, weights), 
                                       lambda = lambdas
      )
      # min.loss.lambda.star <- adapt.fit.star$lambda.min
      # adapt.fit.star.coef <- coef(adapt.fit.star, s="lambda.min")
      adapt.fit.star.coef <- coef(adapt.fit.star, s=min.loss.lambda)
      is_selected_star <- adapt.fit.star.coef[-(1:2),] != 0
      selected_idx_star <- which(is_selected_star)
      alpha_hats_star = alpha_tildes_star[selected_idx_star]
      gamma_hat_star <- adapt.fit.star.coef[2,]
      beta_hats_star <- adapt.fit.star.coef[-(1:2),][selected_idx_star]
      NIE_hat_star <- sum(alpha_hats_star * beta_hats_star)
      
      boot_return <- c(
        gamma_star = gamma_hat_star,
        NIE_star = NIE_hat_star
      )
      return(boot_return)
    }
    
    boot_ci <- apply(boot_results, 2, quantile, probs=c(0.025, 0.975))
    
    return_list <- c(
      return_list, 
      list(boot_ci=boot_ci)
    )
  }
  
  return(
    return_list
  )
  
}

#####################################################
# Bootstrap w/out Selection
#####################################################
cv_sl_estimates_no_sel <- function(sel_df, mu.hats,
                                   # weight.version, weight.gam=1,
                                   # lambdas, folds,
                                   candidate.mediators, x.cols, 
                                   treat.col, outcome.col,
                                   num_bootstraps=1000L,
                                   do.boot=TRUE,
                                   boot.seed=NA) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  .require(dplyr)
  .require(SuperLearner)
  .require(foreach)
  
  n <- nrow(sel_df)
  
  # selection data
  X <- sel_df %>% dplyr::select(all_of(x.cols))
  D <- clean_D(sel_df, treat.col)
  M <- sel_df %>% dplyr::select(all_of(candidate.mediators))
  Y <- sel_df %>% dplyr::pull(outcome.col)
  
  p <- ncol(M)
  
  rob_system <- robinsonize_system(D, M, Y, mu.hats)
  dc <- rob_system$dc
  m_0 <- rob_system$m_0
  y_0 <- rob_system$y_0
  
  mfit <- lm(m_0 ~ dc - 1)
  if(p>1){
    alpha_hats = coef(mfit)[1, ]
  } else {
    alpha_hats = coef(mfit)[1]
  }
  
  yfit = lm(y_0 ~ dc + m_0 - 1)
  gamma_hat <- coef(yfit)[1]
  beta_hats = coef(yfit)[-1]
  NIE_hat <- sum(alpha_hats * beta_hats)
  
  
  
  return_list <- list(
    alpha_hats = alpha_hats, beta_hats = beta_hats,
    NIE_hat = NIE_hat, NDE_hat = gamma_hat,
    ATE_hat = NIE_hat + gamma_hat,
    sel_M = candidate.mediators,
    lambda=0,
    kap=NULL
  )
  
  if(do.boot){
    B=num_bootstraps
    boot_results <- foreach(boot_idx=1:B, .combine=rbind) %do% {
      if(is.numeric(boot.seed))
        set.seed(boot.seed + boot_idx)
      
      G <- rexp(n=n, rate=1)
      
      mfit_star <- lm(m_0 ~ dc - 1, weights=G)
      if(p>1){
        alpha_hats_star = coef(mfit_star)[1, ]
      } else {
        alpha_hats_star = coef(mfit_star)[1]
      }
      
      yfit_star = lm(y_0 ~ dc + m_0 - 1, weights=G)
      gamma_hat_star <- coef(yfit_star)[1]
      beta_hats_star = coef(yfit_star)[-1]
      NIE_hat_star <- sum(alpha_hats * beta_hats)
      
      boot_return <- c(
        NDE_star = gamma_hat_star,
        NIE_star = NIE_hat_star
      )
      return(boot_return)
    }
    
    boot_ci <- apply(boot_results, 2, quantile, probs=c(0.025, 0.975))
    
    return_list <- c(
      return_list, 
      list(boot_ci=boot_ci)
    )
  }
  
  return(
    return_list
  )
  
}



#####################################################
# UPoSI Simultaneous Coverage
#####################################################
cv_sl_estimates_w_uposi <- function(sel_df, mu.hats,
                                   weight.version, lambdas, folds,
                                   candidate.mediators, x.cols, 
                                   treat.col, outcome.col,
                                   do.uposi=TRUE) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  .require(dplyr)
  .require(SuperLearner)
  
  
  n <- nrow(sel_df)
  
  # selection data
  X <- sel_df %>% dplyr::select(all_of(x.cols))
  D <- clean_D(sel_df, treat.col)
  M <- sel_df %>% dplyr::select(all_of(candidate.mediators))
  Y <- sel_df %>% dplyr::pull(outcome.col)
  
  p <- ncol(M)
  
  rob_system <- robinsonize_system(D, M, Y, mu.hats)
  dc <- rob_system$dc
  m_0 <- rob_system$m_0
  y_0 <- rob_system$y_0
  
  mfit <- lm(m_0 ~ dc)
  alpha_tildes = coef(mfit)[2, ]
  summ_mfit <- summary(mfit)
  rm(mfit)
  
  yfit = lm(y_0 ~ dc + m_0)
  beta_tildes = coef(yfit)[-(1:2)]
  
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights <- 1/(abs(beta_tildes))
  )
  
  
  design_mat <- cbind(dc, m_0)
  foldid <- .caretFoldToGlmnetFold(folds, n)
  
  adapt.fit <- glmnet::cv.glmnet(x=design_mat, y=y_0, foldid=foldid, 
                                 intercept=FALSE, 
                                 family = "gaussian", 
                                 penalty.factor = c(0, weights), 
                                 lambda = lambdas
  )
  min.loss.lambda <- adapt.fit$lambda.min
  
  adapt.fit.coef <- coef(adapt.fit, s="lambda.min")
  is_selected <- adapt.fit.coef[-(1:2),] != 0
  selected_idx <- which(is_selected)
  alpha_hats = alpha_tildes[selected_idx]
  
  yfit <- lm(y_0 ~ dc + m_0[,selected_idx])
  yfit.coef <- coef(yfit)
  gamma_hat <- yfit.coef[2,]
  beta_hats <- yfit.coef[-(1:2),]
  NIE_hat <- sum(alpha_hats * beta_hats)
  
  return_list <- list(
    alpha_hats = alpha_hats, beta_hats = beta_hats,
    alpha_tildes=alpha_tildes, beta_tildes=beta_tildes,
    NIE_hat = NIE_hat, NDE_hat = gamma_hat,
    ATE_hat = NIE_hat + gamma_hat,
    sel_M = candidate.mediators[is_selected],
    lambda=min.loss.lambda,
    summ_Y = summary(yfit),
    summ_allM = summ_mfit
  )
  
  if(do.uposi){
    # TODO Do bootstrap
    yhat <- predict(adapt.fit, newx=design_mat, s="lambda.min")
    epshat <- y_0 - yhat
    n <- length(y_0)
    sigma.hat <- sqrt(mean(epshat^2))
    Cn11.inv <- n*MASS::ginv(crossprod(design_mat[, c(T,is_selected)]))
    Cn.inv <- n*MASS::ginv(crossprod(design_mat))
    
    # TODO allow general gamma
    eta_mat <- min.loss.lambda/(2*n)*
      # Cn11.inv %*% 
      sweep(design_mat[, c(T,is_selected)], 
            2, 
            sign(beta_hats)/(abs(beta_hats)^2),
            "*")
    
    Sigma.hat.n <- 1/n*crossprod(
      epshat*(
        design_mat[, c(T,is_selected)] + eta_mat
      ) %*% Cn11.inv
    )
    Sigma.hat.n.half <- chol(Sigma.hat.n)
    Sigma.hat.n.half.inv <- MASS::ginv(Sigma.hat.n.half)
    
    Sigma.breve.n <- 1/n*crossprod(design_mat[, c(T,is_selected)] + eta_mat)
    Sigma.breve.n.half <- chol(Sigma.breve.n)
    Sigma.breve.n.half.inv <- MASS::ginv(Sigma.breve.n.half)
    
    # TODO: Figure out what kind of studentization to do
    sigma.d.hat <- sqrt(mean(dc))
    Sigma.m.breve.n.half <- chol(Sigma.m.breve.n)
    Sigma.m.breve.n.half.inv <- MASS:ginv(Sigma.m.breve.n.half)
    
    
    B=1000
    boot_results <- foreach(boot.idx=1:B, .combine=rbind) %do% {
      if(exists("boot.seed"))
        set.seed(boot.seed + boot.idx)
      G <- rbeta(n, shape1=1/2, shape2=3/2)
      muG <- .25 #.5/(.5 + 1.5)
      Gtil <- (G - muG)/muG
      
      etastar <- sweep(etahat, 1, Gtil, "*")
      epsstar <- epshat*Gtil
      
      m_0_star <- mhat + etastar
      y_0_star <- yhat + epsstar
      
      mfit_star <- lm(m_0_star ~ dc)
      alpha_tildes_star = coef(mfit_star)[2, ]
      
      yfit = lm(y_0_star ~ design_mat - 1)
      beta_tildes_star = coef(yfit)[-1]
      
      switch (weight.version,
              "mixture" = weights <- 1/(abs(beta_tildes_star) * (1 + abs(alpha_tildes_star))),
              "product" = weights <- 1/(abs(beta_tildes_star) * abs(alpha_tildes_star)),
              "adaptive"= weights <- 1/(abs(beta_tildes_star))
      )
      
      adapt.fit.star <- glmnet::cv.glmnet(x=design_mat, y=y_0_star, foldid=foldid, 
                                          intercept=FALSE, 
                                          family = "gaussian", 
                                          penalty.factor = c(0, weights), 
                                          lambda = lambdas
      )
      min.loss.lambda.star <- adapt.fit.star$lambda.min
      adapt.fit.star.coef <- coef(adapt.fit.star, s="lambda.min")
      is_selected_star <- adapt.fit.star.coef[-(1:2),] != 0
      selected_idx_star <- which(is_selected_star)
      alpha_hats_star = alpha_tildes_star[selected_idx_star]
      gamma_hat_star <- adapt.fit.star.coef[2,]
      beta_hats_star <- adapt.fit.star.coef[-(1:2),][selected_idx_star]
      NIE_hat_star <- sum(alpha_hats_star * beta_hats_star)
      
      # TODO: Enable studentized pivots
      # y_hat_star <- predict(adapt.fit.star, newx=design_mat, s="lambda.min")
      # eps_hat_star <- y_0_star - y_hat_star
      # sigma_hat_star <- sqrt(mean((eps_hat_star*Gtil)^2))
      # 
      # sigma_m_hat_star <- sqrt(mean((mfit_star$residuals*Gtil)^2))
      # 
      # Sigma.m.breve.star <- vcov(mfit_star)[2*(1:p),2*(1:p)]
      # 
      # Tstar_alpha <- alpha_hats_star - alpha_hats
      # Tstar_beta <- beta_hats_star - beta_hats
      # Tstar_gamma <- gamma_hat_star - gamma_hat
      # 
      # pivot_outcome_model <- sigma.hat/sigma_hat_star*Sigma.breve.n.half.inv %*% 
      #   c(Tstar_gamma, Tstar_beta)
      # 
      # pivot_mediators <- sigma.m.hat/sigma_m_hat_star*Sigma.m.breve.n.half.inv %*%
      #   Tstar_alpha
      
      # boot_return <- list(
      #   alpha_star = alpha_hats_star,
      #   beta_star = beta_hats_star,
      #   gamma_star = gamma_hat_star,
      #   NIE_star = NIE_hat_star,
      #   is_selected_star = is_selected_star
      # )
      boot_return <- c(
        gamma_star = gamma_hat_star,
        NIE_star = NIE_hat_star
      )
      return(boot_return)
    }
    
    boot_ci <- apply(boot_results, 2, quantile, probs=c(0.025, 0.975))
    
    return_list <- c(
      return_list, 
      list(boot_ci=boot_ci)
    )
  }
  
  return(
    return_list
  )
  
}

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
    
    mu.mxis = foreach(i=1:p) %dopar% {
      
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
      mu.mxi.coef = mu.mxi$coef
      outm <- list(mu.mxi=mu.mxi.pred,
                   mu.mxi.coef=mu.mxi.coef)
      if(!save_pred_only){
        mu.mxi.lib.pred = mu.mxi$library.predict
        mu.mxi.method = mu.mxi$method
        outm <- c(outm,
                  list(
                    mu.mxi.lib.pred=mu.mxi.lib.pred,
                    mu.mxi.method=mu.mxi.method
                  ))
      }
      
      return(outm)
    }
    mu.mxi.preds = lapply(mu.mxis, function(x) x$mu.mxi)
    mu.mxi.coefs = lapply(mu.mxis, function(x) x$mu.mxi.coef)
    out <- c(
      out, 
      list(mu.mxis=mu.mxi.preds,
           mu.mxi.coefs=mu.mxi.coefs
      )
    )
    if(!save_pred_only){
      mu.mxi.lib.preds = lapply(mu.mxis, function(x) x$mu.mxi.lib.pred)
      mu.mxi.methods = lapply(mu.mxis, function(x) x$mu.mxi.method)
      
      out <- c(
        out, 
        list(
          mu.mxi.lib.preds=mu.mxi.lib.preds,
          mu.mxi.methods=mu.mxi.methods
        )
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
    mu.dx.coef = mu.dx$coef
    out <- c(
      out,
      list(mu.dx=mu.dx.pred,
           mu.dx.coef=mu.dx.coef
      )
    )
    
    if(!save_pred_only){
      mu.dx.lib.pred = mu.dx$library.predict
      mu.dx.method = mu.dx$method
      out <- c(
        out,
        list(
          mu.dx.lib.pred=mu.dx.lib.pred,
          mu.dx.method=mu.dx.method
        )
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
    mu.yx.coef = mu.yx$coef
    out <- c(
      out,
      list(
        mu.yx=mu.yx.pred,
        mu.yx.coef=mu.yx.coef
      )
    )
    
    if(!save_pred_only){
      mu.yx.lib.pred = mu.yx$library.predict
      mu.yx.method = mu.yx$method
      
      out <- c(
        out,
        list(
          mu.yx.lib.pred=mu.yx.lib.pred,
          mu.yx.method=mu.yx.method
        )
      )
    }
    
    rm(mu.yx)
    
  }
  
  return(out)
}

my_sl_method_compute <- function(SL.method, library.predict, outcome, obsWeights,
                                 new.library.predict=library.predict){
  coef <- SL.method$computeCoef(library.predict, outcome, 
                                libraryNames=colnames(library.predict),
                                verbose=F, obsWeights=obsWeights,
                                control=list(saveFitLibrary = TRUE, trimLogit = 0.001)
  )$coef
  
  preds <- SL.method$computePred(new.library.predict, coef,
                                 control=list(saveFitLibrary = TRUE, trimLogit = 0.001)
  )
}

my_sl_split_compute <- function(split.folds, SL.method, library.predict, outcome, obsWeights){
  .require(foreach)
  stopifnot(is.list(split.folds))
  preds <- foreach(fold=split.folds, .combine=c) %do% {
    my_sl_method_compute(SL.method, library.predict[fold,], outcome[fold], obsWeights[fold],
                         library.predict[-fold,])
  }
  
  idxs <- do.call(c, split.folds)
  sorted.preds <- rep(0, length(split_folds))
  sorted.preds[idxs] <- preds
  return(sorted.preds)
}

multiplier_bootstrap_recombine_SL <- function(split.folds, mu.hats,
                                              Y, M, D, obsWeights){
  .require(foreach)
  .require(SuperLearner)
  .require(dplyr)
  .require(tibble)
  
  mu.yx.method <- mu.hats$mu.yx.method
  mu.yx.lib.pred <- mu.hats$mu.yx.lib.pred
  mu.dx.method <- mu.hats$mu.dx.method
  mu.dx.lib.pred <- mu.hats$mu.dx.lib.pred
  mu.mxi.methods <- mu.hats$mu.mxi.methods
  mu.mxi.lib.preds <- mu.hats$mu.mxi.lib.preds
  
  # Y
  Y.pred <- my_sl_split_compute(split.folds, 
                                SL.method=mu.yx.method, 
                                library.predict=mu.yx.lib.pred, 
                                outcome=Y,
                                obsWeights=obsWeights)
  new_y_0 <- Y - Y.pred
  
  # D
  D.pred <- my_sl_split_compute(split.folds, 
                                SL.method=mu.dx.method, 
                                library.predict=mu.dx.lib.pred, 
                                outcome=D,
                                obsWeights=obsWeights)
  new_dc <- D - D.pred
  
  # M
  new_m_0 <- foreach(j=1:ncol(M), .combine=cbind) %do% {
    Mj.pred <- my_sl_split_compute(split.folds, 
                                  SL.method=mu.mxi.methods[[j]],
                                  library.predict=mu.mxi.lib.preds[[j]], 
                                  outcome=dplyr::pull(M, j),
                                  obsWeights=obsWeights)
    dplyr::pull(M, j) - Mj.pred
  }
  
  tibble::tibble(y_0=new_y_0, dc=new_dc, m_0=new_m_0)
}

backfit_psi_m <- function(sel_df, mu.hats, candidate.mediators,
                          treat.col){
  
  p = length(candidate.mediators)
  
  mu.mxis <- mu.hats$mu.mxis
  mu.dx <- mu.hats$mu.dx
  mu.yx <- mu.hats$mu.yx
  
  # tests
  D <- clean_D(sel_df, treat.col)
  M <- sel_df %>% dplyr::select(all_of(candidate.mediators))
  stopifnot(p == ncol(M))
  
  rob_system <- robinsonize_system(D, M, Y=NULL, mu.hats)
  dc <- rob_system$dc
  m_0 <- rob_system$m_0
  
  mfit <- lm(m_0 ~ dc)
  alpha_tildes = coef(mfit)[2, ]
  
  m_1 = sapply(1:p, function(i){ dplyr::pull(M, i) - D*alpha_tildes[i] })
  if(!is.data.frame(m_1)) {
    m_1 <- data.frame(m_1)
    colnames(m_1) <- candidate.mediators
  }
  return(list(m_1=m_1,
              alpha_tildes=alpha_tildes))
}