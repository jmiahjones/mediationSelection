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
  
  # m_0 = sapply(1:p, function(i){
  m_0 <- foreach(i=1:p, .combine=cbind) %do% {
    # estimate the counfounding effect on mediatior i
    
    Mi <- dplyr::pull(M, i)
    mu.mxi = mu.mxis[[i]]
    
    col <- Mi - mu.mxi
    return(col)
  # })
  }
  stopifnot(ncol(m_0)==p)
  
  rob_tbl <- tibble::tibble(
    dc=dc,
    m_0=m_0
  )
  
  if(!is.null(Y)){
    mu.yx <- mu.hats$mu.yx
    y_0 = as.numeric(Y - mu.yx)
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

######## With Selection #####################
cv_sl_estimates_w_sel <- function(y_0, m_0, dc,
                                  weight.version, weight.gam=1,
                                  lambdas, folds,
                                  opt=NULL,
                                  do.boot=TRUE,
                                  num_bootstraps=1000,
                                  boot.seed=NA) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  .require(dplyr)
  .require(tibble)
  .require(foreach)
  
  assertthat::assert_that(
    is.matrix(m_0)
  )
  n <- nrow(m_0)
  p <- ncol(m_0)
  
  assertthat::assert_that(
    n == length(dc), 
    n == length(y_0),
    p>1,
    is.numeric(dc) || (is.matrix(dc) && ncol(dc) == 1),
    is.numeric(y_0) || (is.matrix(y_0) && ncol(y_0) == 1)
  )
  
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
  
  # TODO: Remove hard-coded true m set
  num_missed <- length(setdiff(1:3, selected_idx))
  num_noise <- length(setdiff(selected_idx, 1:3))
  
  sel_size <- length(selected_idx)
  return_tbl <- tibble(
    n=n,
    p=p,
    opt,
    model_version = weight.version,
    lambda=min.loss.lambda,
    NIE_hat = NIE_hat, 
    NDE_hat = gamma_hat,
    ATE_hat = NIE_hat + gamma_hat,
    sel_size = sel_size,
    kap=best.kap,
    sel_info=list(tibble(
      alpha_hats=alpha_hats, beta_hats=beta_hats,
      sel_M = selected_idx
    )),
    num_missed = num_missed,
    num_noise = num_noise
  )
  
  if(do.boot){
    
    yhat <- predict(adapt.fit, newx=design_mat, s="lambda.min")
    epshat <- y_0 - yhat
    
    das_boot_ci <- das_boot_fun(
      design_mat,
      y_0, m_0, dc,
      yhat, mhat,
      etahat,
      epshat,
      weight.version,
      min.loss.lambda,
      best.kap,
      lambdas,
      B=num_bootstraps,
      boot.seed=boot.seed
    )
    
    minnier_boot_ci <- minnier_boot_fun(
      design_mat,
      y_0, m_0, dc,
      weight.version,
      min.loss.lambda,
      best.kap,
      lambdas,
      B=num_bootstraps,
      boot.seed=boot.seed
    )
    
    naive_boot <- naive_boot_fun(y_0, m_0, dc,
                                 selected_idx,
                                 B=num_bootstraps,
                                 boot.seed=boot.seed
    )
    
    naive_delta <- naive_delta_inf(y_0, m_0, dc, selected_idx)
    
    boot_tbl <- rbind(
      # selection inf methods
      tibble(
        target=c("NDE","NIE"),
        inf_method="minnier",
        lower=minnier_boot_ci[,1],
        upper=minnier_boot_ci[,2],
      ),
      tibble(
        target=c("NDE","NIE"),
        inf_method="das",
        lower=das_boot_ci[,1],
        upper=das_boot_ci[,2]
      ),
      # naive methods
      tibble(
        target=c("NDE","NIE"),
        inf_method="naive_boot",
        lower=naive_boot[,1],
        upper=naive_boot[,2],
      ),
      tibble(
        target=c("NDE","NIE"),
        inf_method="naive_delta",
        lower=naive_delta[,1],
        upper=naive_delta[,2]
      )
    ) %>% 
      tidyr::pivot_wider(
        names_from="target",
        values_from=c("lower", "upper")
      )
    
    return_tbl <- tibble(return_tbl, boot_tbl)
  }
  
  return(
    return_tbl
  )
  
}


########## No Selection ####################
cv_sl_estimates_no_sel <- function(y_0, m_0, dc,
                                   model_name,
                                   opt=NULL,
                                   num_bootstraps=1000L,
                                   do.boot=TRUE,
                                   boot.seed=NA) {
  
  .require(dplyr)
  .require(tibble)
  .require(foreach)
  
  assertthat::assert_that(
    is.matrix(m_0)
  )
  n <- nrow(m_0)
  p <- ncol(m_0)
  
  assertthat::assert_that(
    n == length(dc), 
    n == length(y_0),
    p>1,
    is.numeric(dc) || (is.matrix(dc) && ncol(dc) == 1),
    is.numeric(y_0) || (is.matrix(y_0) && ncol(y_0) == 1)
  )
  
  
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
  
  # TODO: Remove hard-code
  num_missed <- 0
  num_noise <- p-3
  
  return_tbl <- tibble(
    n=n,
    p=p,
    opt,
    model_version = model_name,
    lambda = NA,
    NIE_hat = NIE_hat, 
    NDE_hat = gamma_hat,
    ATE_hat = NIE_hat + gamma_hat,
    sel_size = NA,
    kap=NA,
    sel_info=list(tibble(
      alpha_hats=alpha_hats, beta_hats=beta_hats,
      sel_M = NA
    )),
    num_missed = num_missed,
    num_noise = num_noise
  )
  
  if(do.boot){
    
    naive_boot <- naive_boot_fun(y_0, m_0, dc,
                                 1:p,
                                 B=num_bootstraps,
                                 boot.seed=boot.seed
    )
    
    naive_delta <- naive_delta_inf(y_0, m_0, dc, 1:p)
    
    boot_tbl <- rbind(
      # naive methods
      tibble(
        target=c("NDE","NIE"),
        inf_method="naive_boot",
        lower=naive_boot[,1],
        upper=naive_boot[,2],
      ),
      tibble(
        target=c("NDE","NIE"),
        inf_method="naive_delta",
        lower=naive_delta[,1],
        upper=naive_delta[,2]
      )
    ) %>% 
      tidyr::pivot_wider(
        names_from="target",
        values_from=c("lower", "upper")
      )
    
    return_tbl <- tibble(return_tbl, boot_tbl)
  }
  
  return(
    return_tbl
  )
  
}


####### Inference Helpers ##################################

# Lasso-Boostrap Helpers ###################################
das_boot_fun <- function(design_mat,
                         y_0, m_0, dc,
                         yhat, mhat,
                         etahat,
                         epshat,
                         weight.version,
                         min.loss.lambda,
                         best.kap,
                         lambdas,
                         B,
                         boot.seed=NULL){
  
  # TODO Do Studentization from Das
  n <- length(y_0)
  
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
    
    boot_return <- c(
      gamma_star = gamma_hat_star,
      NIE_star = NIE_hat_star
    )
    return(boot_return)
  }
  
  boot_ci <- foreach(i=1:ncol(boot_results), .combine=rbind) %do%{
    quantile(boot_results[,i], probs=c(0.025, 0.975))
  }
  
  return(boot_ci)
  
}


minnier_boot_fun <- function(design_mat,
                             y_0, m_0, dc,
                             weight.version,
                             min.loss.lambda,
                             best.kap,
                             lambdas,
                             B,
                             boot.seed=NULL){
  
  
  boot_results <- foreach(boot_idx=1:B, .combine=rbind) %do% {
    if(is.numeric(boot.seed))
      set.seed(boot.seed + boot_idx)
    
    G <- rgamma(n=n, shape=1, rate=1)
    
    mfit <- lm(m_0 ~ dc, weights = G)
    alpha_tildes_star = coef(mfit)[2, ]
    
    yfit = lm(y_0 ~ design_mat - 1, weights = G)
    beta_tildes_star = coef(yfit)[-1]
    
    switch (weight.version,
            "mixture" = weights <- 1/(abs(beta_tildes_star) * (1 + abs(alpha_tildes_star))),
            "product" = weights <- 1/(abs(beta_tildes_star) * abs(alpha_tildes_star)),
            "adaptive"= weights <- 1/(abs(beta_tildes_star))
    )
    weights <- weights^(best.kap)
    
    adapt.fit.star <- glmnet::glmnet(x=design_mat, y=y_0,
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
    
    boot_result <- c(
      gamma_star = gamma_hat_star,
      NIE_star = NIE_hat_star
    )
    return(boot_result)
  }
  
  boot_ci <- foreach(i=1:ncol(boot_results), .combine=rbind) %do%{
    quantile(boot_results[,i], probs=c(0.025, 0.975))
  }
  
  return(boot_ci)
  
}

# Naive Helpers #########################################

naive_boot_fun <- function(y_0, m_0, dc, 
                           sel_M_idxs,
                           B,
                           boot.seed=NULL){
  
  boot_results <- foreach(boot_idx=1:B, .combine=rbind) %do% {
    if(is.numeric(boot.seed))
      set.seed(boot.seed + boot_idx)
    
    p_sel <- length(sel_M_idxs)
    m_sel <- as.matrix(m_0[,sel_M_idxs])
    
    G <- rmultinom(1, n, prob=rep(1/n, n))
    
    mfit <- lm(m_sel ~ dc, weights = G)
    
    if(p_sel > 1){
      alpha_hats_star <- coef(mfit)[2, ]
    } else {
      alpha_hats_star <- coef(mfit)[2]
    }
    
    yfit <- lm(y_0 ~ dc + m_sel - 1, weights = G)
    beta_hats_star <- coef(yfit)[-1]
    gamma_hat_star <- coef(yfit)[1]
    
    NIE_hat_star <- sum(alpha_hats_star * beta_hats_star)
    
    boot_result <- c(
      gamma_star = gamma_hat_star,
      NIE_star = NIE_hat_star
    )
    return(boot_result)
  }
  
  boot_ci <- foreach(i=1:ncol(boot_results), .combine=rbind) %do%{
    quantile(boot_results[,i], probs=c(0.025, 0.975))
  }
  
  return(boot_ci)
  
}


naive_delta_inf <- function(y_0, m_0, dc, sel_M_idxs, ret_est=F){
  p_sel <- length(sel_M_idxs)
  m_sel <- as.matrix(m_0[,sel_M_idxs])
  even_idx_sel <- 2*(1:p_sel)
  sel_naive_fit_y <- lm(y_0 ~ dc + m_sel)
  sel_naive_fit_m <- lm(m_sel ~ dc)
  nde_se <- sqrt(vcov(sel_naive_fit_y)[1,1])
  beta_hat <- coef(sel_naive_fit_y)[-(1:2)]
  if(p_sel > 1){
    alpha_hat <- coef(sel_naive_fit_m)[2,]
  } else {
    alpha_hat <- coef(sel_naive_fit_m)[2]
  }
  
  beta_vcov <- vcov(sel_naive_fit_y)[-(1:2),-(1:2)]
  alpha_vcov <- vcov(sel_naive_fit_m)[even_idx_sel,even_idx_sel]
  # browser()
  nie_var <- crossprod(alpha_hat, beta_vcov %*% alpha_hat) +
    crossprod(beta_hat, alpha_vcov %*% beta_hat)
  nie_se <- sqrt(nie_var)
  nde_hat <- coef(sel_naive_fit_y)[2] %>% unname
  nde_ci <- nde_hat + c(-1,1)*1.96*nde_se
  
  nie_hat <- sum(alpha_hat * beta_hat)
  nie_ci <- nie_hat + c(-1,1)*c(1.96*nie_se)
  
  ci_mat <- rbind(nde_ci, nie_ci)
  colnames(ci_mat) <- c("NDE", "NIE")
  if(ret_est){
    return(
      list(est=c(NDE=nde_hat, NIE=nie_hat),
           ci=ci_mat)
    )
  } else {
    return(ci_mat)
  }
}

