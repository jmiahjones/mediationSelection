# pseudo-mse-with-semi-scripts.R
# Functions that perform the estimation of NDE and NIE using semiparametric confounding control.
# Includes both the double-robust linear technique, as well as the svm method.

#################################################################
### Methods for computing cross-validation
#################################################################



## DR Method
dr_lin_estimates <- function(sel_data, est_data, weight.version, lambdas) { # double robust approach
  
  X <- dplyr::select(sel_data, c("x.1", "x.2", "x.3", "x.4", "x.5")) %>% as.matrix()
  XD <- dplyr::select(sel_data, c("x.1", "x.2", "x.3", "x.4", "x.5", "d")) %>% as.matrix()
  M <- dplyr::select(sel_data, candidate.mediators) %>% as.matrix()
  
  # for the estimation set, we don't train any predictors in lm, so keep these as data frames
  X.est <- dplyr::select(est_data, c("x.1", "x.2", "x.3", "x.4", "x.5")) %>% as.matrix()
  XD.est <- dplyr::select(est_data, c("x.1", "x.2", "x.3", "x.4", "x.5", "d")) %>% as.matrix()
  M.est <- dplyr::select(est_data, candidate.mediators) %>% as.matrix()
  
  
  # remove confounding effect from mediators
  mu.mx = lm(M~X)
  m_0 = M - predict(mu.mx)
  
  # remove confounding effect from treatment
  mu.dx = glm(d~., as.data.frame(XD), family = "binomial")
  dc = sel_data$d - predict(mu.dx, type="response")
  
  # first-stage: alphas
  alpha_tildes = coef(lm(m_0 ~ dc + 0))
  
  # first-stage: betas
  # remove treatment effect from outcome and mediators
  mu.yxd = lm(sel_data$y~XD)
  y_1 = sel_data$y - predict(mu.yxd)
  
  mu.mxd = lm(M~XD)
  m_1 = M - predict(mu.mxd)
  
  beta_tildes = coef(lm(y_1 ~ m_1 + 0))
  
  # perform variable selection
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights <- 1/(abs(beta_tildes))
  )
  
  # regress mediators w/ d,x removed on y w/ d,x removed
  weight.lasso.fit <- glmnet::glmnet(m_1, y_1, intercept=FALSE,
                                     family = "gaussian", 
                                     penalty.factor = weights,
                                     lambda = lambdas)
  
  
  # now do the centering on the estimation data
  # for the alphas
  m_0.est = M.est - mypredict.lm(mu.mx, X.est)
  dc.est = est_data$d - predict(mu.dx, as.data.frame(X.est), type="response")
  
  alpha_tildes.est = coef(lm(m_0.est ~ dc.est + 0))
  
  # for the betas
  y_1.est = est_data$y - mypredict.lm(mu.yxd, XD.est)
  
  result <- lapply(lambdas, function(lam)
  {
    weight.lasso.coef <- coef(weight.lasso.fit, lam)
    is_selected <- weight.lasso.coef[-1,] != 0
    
    if(!any(is_selected)){
      alpha_hats = NA
      beta_hats = NA
      NIE_hat = NA
      theta_hat = NA
    } else {
      selected_idx <- which(is_selected)
      
      ## first stage after variable selection
      ## must use the estimation data set
      
      # center the new treatments
      alpha_hats = alpha_tildes.est[selected_idx]
      
      m_1.est = M.est - mypredict.lm(mu.mxd, XD.est)
      m_1.est = m_1.est[, selected_idx]
      
      beta_hats = coef(lm(y_1.est ~ m_1.est + 0))
      
      ## second stage, get theta
      y.prime = est_data$y - m_0.est[, selected_idx] %*% matrix(beta_hats, ncol=1)
      
      
      ### train mu.yprime.x on the larger selection data set
      y.prime.sel = sel_data$y - m_0[, selected_idx] %*% matrix(beta_hats, ncol=1)
      mu.yprime.x = lm(y.prime.sel ~ X)
      ### end of training
      
      y_2 = y.prime - mypredict.lm(mu.yprime.x, X.est)
      
      theta_hat <- coef(lm(y_2 ~ dc.est + 0))
      NIE_hat <- sum(alpha_hats * beta_hats)
    }
    
    
    return(
      list(alpha_hats = alpha_hats, beta_hats = beta_hats,
           NIE_hat = NIE_hat, NDE_hat = theta_hat,
           is_selected = is_selected)
    )
  })
  
  theta_hats = purrr::map_dbl(result, function(lambda) lambda$NDE_hat)
  NIE_hats = purrr::map_dbl(result, function(lambda) lambda$NIE_hat)
  return(
    list(full_result = result,
         theta_hats = theta_hats,
         NIE_hats = NIE_hats,
         is_selected = result$is_selected)
  )
}

#################################################################
### Methods for computing unbiased, noisy estimates
#################################################################
# function using linear models
est_full_drlin_nde <- function(data, candidate.mediators){
  
  X <- dplyr::select(data, c("x.1", "x.2", "x.3", "x.4", "x.5")) %>% as.matrix()
  XD <- dplyr::select(data, c("x.1", "x.2", "x.3", "x.4", "x.5", "d")) %>% as.matrix()
  M <- dplyr::select(data, candidate.mediators) %>% as.matrix()
  
  # remove confounding effect from mediators
  mu.mx = lm(M~X)
  m_0 = M - predict(mu.mx)
  
  # remove confounding effect from treatment
  mu.dx = glm(d~., as.data.frame(XD), family = "binomial")
  dc = data$d - predict(mu.dx, type="response")
  
  # first-stage: alphas
  alpha_hats = coef(lm(m_0 ~ dc + 0))
  
  # first-stage: betas
  # remove treatment effect from outcome and mediators
  mu.yxd = lm(data$y~XD)
  y_1 = data$y - predict(mu.yxd)
  
  # m_1 = m_0 - matrix(data$d, ncol=1) %*% matrix(alpha_hats, nrow=1)
  mu.mxd = lm(M~XD)
  m_1 = M - predict(mu.mxd)
  
  beta_hats = coef(lm(y_1 ~ m_1 + 0))
  
  
  y.prime = data$y - m_0 %*% matrix(beta_hats, ncol=1)
  mu.yprime.x = lm(y.prime ~ X)
  y_2 = y.prime - predict(mu.yprime.x)
  
  theta_hat <- coef(lm(y_2 ~ dc + 0))
  
  # return(list(alpha_hats, beta_hats, theta_hat))
  return(theta_hat)
  
  # NIE_hat = sum(alpha_hats * beta_hats)
  # return(NIE_hat)
}


##### Cross-validated SuperLearner #####

# shorter version
cv_sl_estimates <- function(sel_df, mu.hats,
                            weight.version, lambdas, folds,
                            candidate.mediators, x.cols, treat.col, outcome.col) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  require(dplyr)
  require(SuperLearner)
  
  
  p = length(candidate.mediators)
  V <- length(folds)
  
  mu.mxis <- mu.hats$mu.mxis
  mu.dx <- mu.hats$mu.dx
  mu.yx <- mu.hats$mu.yx
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  m_0 = sapply(1:p, function(i){
    # estimate the counfounding effect on mediatior i
    mi = candidate.mediators[i]
    Mi <- dplyr::pull(sel_df, mi)
    mu.mxi = mu.mxis[[i]]
    
    col <- Mi - my_sl_predict(mu.mxi)
    return(col)
  })
  stopifnot(ncol(m_0)==p)
  
  d.prob = as.numeric(my_sl_predict(mu.dx))
  dc = sel_df[,treat.col] - d.prob
  
  alpha_tildes = coef(lm(m_0 ~ dc))[2, ]
  
  y_0 = sel_df[,outcome.col] - my_sl_predict(mu.yx)
  yfit = lm(y_0 ~ dc + m_0)
  beta_tildes = coef(yfit)[-(1:2)]
  
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights <- 1/(abs(beta_tildes))
          # "both"    = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes)))
  )
  
  print("Running Adaptive Lasso")
  
  relaxo.loss <- function(x, y, intercept=FALSE, family, penalty.factor, lambda,
                          holdout.x, holdout.y){
    lasso.fit <- glmnet::glmnet(x, y, intercept=intercept,
                                family = family, 
                                penalty.factor = penalty.factor,
                                lambda = lambda)
    
    beta_0 <- rep(0, ncol(x)+1)
    lapply(lambda, function(lam){
      lasso.coef <- coef(lasso.fit, s=lam)
      is_selected <- lasso.coef[-1,1] != 0
      beta <- beta_0
      beta[c(T, is_selected)] <- coef(lm(y~x[,is_selected]))
      
      error <- holdout.y - beta[1] - as.numeric(holdout.x %*% matrix(beta[-1], ncol=1))
      loss <- error^2
      return(loss)
    }) %>% do.call(cbind, .)
  }
  
  design_mat <- cbind(dc, m_0)
  losses <- vector("list", V)
  for(fold.idx in 1:V){
    holdout.idx <- folds[[fold.idx]]
    ## perform variable selection
    # regress mediators w/ d,x removed on y w/ d,x removed
    holdout.mat <- design_mat[holdout.idx,]
    holdout.y <- y_0[holdout.idx]
    train.mat <- design_mat[-holdout.idx,]
    train.y <- y_0[-holdout.idx]
    losses[[fold.idx]] <- relaxo.loss(x=train.mat, y=train.y, intercept=FALSE, 
                                      family = "gaussian", 
                                      penalty.factor = c(0, weights), 
                                      lambda = lambdas, holdout.x=holdout.mat,
                                      holdout.y = holdout.y)
  }
  loss <- do.call(rbind, losses) %>% colMeans()
  
  min.loss.lambda <- lambdas[which.min(loss)]
  
  adapt.fit <- glmnet::glmnet(x=design_mat, y=y_0, intercept=F, 
                              family="gaussian", penalty.factor = c(0, weights),
                              lambda=min.loss.lambda*(2^(-1:3)))
  
  adapt.fit.coef <- coef(adapt.fit, s=min.loss.lambda)
  is_selected <- adapt.fit.coef[-(1:2),] != 0
  
  total_effect = coef(lm(y_0 ~ dc))[-1]
  
  # print("Getting Estimated Models")
  
  if(!any(is_selected)){
    alpha_hats = NA
    beta_hats = NA
    NIE_hat = NA
    gamma_hat = total_effect
  } else {
    selected_idx <- which(is_selected)
    
    ## first stage after variable selection
    ## must use the estimation data set
    
    # center the new treatments
    m_0_selected = m_0[, selected_idx]
    # alpha_hats = coef(lm(m_0_selected ~ dc))[2, ]
    alpha_hats = alpha_tildes[selected_idx]
    
    yfit = lm(y_0 ~ dc + m_0_selected)
    gamma_hat <- coef(yfit)[2]
    beta_hats <- coef(yfit)[-(1:2)]
    
    NIE_hat <- sum(alpha_hats * beta_hats)
  }
  
  
  return(
    list(alpha_hats = alpha_hats, beta_hats = beta_hats,
         NIE_hat = NIE_hat, NDE_hat = gamma_hat,
         sel_M = candidate.mediators[is_selected],
         losses=loss,
         lambda=min.loss.lambda)
  )
  
}

#####################################################
# Cross-validation: Compute estimates at each lambda
#####################################################
# shorter version
cv_sl_estimates_tuning <- function(sel_df, mu.hats,
                                   weight.version, lambdas, folds,
                                   candidate.mediators, x.cols, treat.col, outcome.col) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  require(dplyr)
  require(SuperLearner)
  
  
  p = length(candidate.mediators)
  V <- length(folds)
  
  mu.mxis <- mu.hats$mu.mxis
  mu.dx <- mu.hats$mu.dx
  mu.yx <- mu.hats$mu.yx
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  m_0 = sapply(1:p, function(i){
    # estimate the counfounding effect on mediatior i
    mi = candidate.mediators[i]
    Mi <- dplyr::pull(sel_df, mi)
    mu.mxi = mu.mxis[[i]]
    
    col <- Mi - my_sl_predict(mu.mxi)
    return(col)
  })
  stopifnot(ncol(m_0)==p)
  
  d.prob = as.numeric(my_sl_predict(mu.dx))
  dc = sel_df[,treat.col] - d.prob
  
  alpha_tildes = coef(lm(m_0 ~ dc))[2, ]
  
  y_0 = sel_df[,outcome.col] - my_sl_predict(mu.yx)
  yfit = lm(y_0 ~ dc + m_0)
  beta_tildes = coef(yfit)[-(1:2)]
  
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights <- 1/(abs(beta_tildes))
          # "both"    = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes)))
  )
  
  print("Running Adaptive Lasso")
  
  relaxo.loss <- function(x, y, intercept=FALSE, family, penalty.factor, lambda,
                          holdout.x, holdout.y){
    lasso.fit <- glmnet::glmnet(x, y, intercept=intercept,
                                family = family, 
                                penalty.factor = penalty.factor,
                                lambda = lambda)
    
    beta_0 <- rep(0, ncol(x)+1)
    lapply(lambda, function(lam){
      lasso.coef <- coef(lasso.fit, s=lam)
      is_selected <- lasso.coef[-1,1] != 0
      beta <- beta_0
      relax.fit <- lm(y~x[,is_selected])
      epsilons <- residuals(relax.fit)
      beta[c(T, is_selected)] <- coef(relax.fit)
      error <- holdout.y - beta[1] - as.numeric(holdout.x %*% matrix(beta[-1], ncol=1))
      # loss <- error^2
      
      ## ASSUME: D is in the first column in X, so X[,-1] is the centered M
      # M <- X[,-1]
      loss <- sapply(1:p, function(j){
        abs(alpha_tildes[j] * sum(error*holdout.x[,j+1]))
      }) %>% sum
      return(loss)
    }) %>% do.call(cbind, .)
  }
  
  design_mat <- cbind(dc, m_0)
  losses <- vector("list", V)
  for(fold.idx in 1:V){
    holdout.idx <- folds[[fold.idx]]
    ## perform variable selection
    # regress mediators w/ d,x removed on y w/ d,x removed
    holdout.mat <- design_mat[holdout.idx,]
    holdout.y <- y_0[holdout.idx]
    train.mat <- design_mat[-holdout.idx,]
    train.y <- y_0[-holdout.idx]
    losses[[fold.idx]] <- relaxo.loss(x=train.mat, y=train.y, intercept=FALSE, 
                                      family = "gaussian", 
                                      penalty.factor = c(0, weights), 
                                      lambda = lambdas, holdout.x=holdout.mat,
                                      holdout.y = holdout.y)
  }
  loss <- do.call(rbind, losses) %>% colSums()
  
  min.loss.lambda <- lambdas[which.min(loss)]
  total_effect = coef(lm(y_0 ~ dc))[-1]
  
  results_per_lambda <- lapply(lambdas, function(lam){
    adapt.fit <- glmnet::glmnet(x=design_mat, y=y_0, intercept=F, 
                              family="gaussian", penalty.factor = c(0, weights),
                              lambda=lam*(2^(-1:2)))
  
  	adapt.fit.coef <- coef(adapt.fit, s=lam)
    is_selected <- adapt.fit.coef[-(1:2),] != 0
    
    if(!any(is_selected)){
      alpha_hats = NA
      beta_hats = NA
      NIE_hat = NA
      gamma_hat = total_effect
    } else {
      selected_idx <- which(is_selected)
      
      ## first stage after variable selection
      ## must use the estimation data set
      
      # center the new treatments
      m_0_selected = m_0[, selected_idx]
      # alpha_hats = coef(lm(m_0_selected ~ dc))[2, ]
      alpha_hats = alpha_tildes[selected_idx]
      
      yfit = lm(y_0 ~ dc + m_0_selected)
      gamma_hat <- coef(yfit)[2]
      beta_hats <- coef(yfit)[-(1:2)]
      
      NIE_hat <- sum(alpha_hats * beta_hats)
    }
    
    
    return(
      list(alpha_hats = alpha_hats, beta_hats = beta_hats,
           NIE_hat = NIE_hat, NDE_hat = gamma_hat,
           sel_M = candidate.mediators[is_selected],
           lambda=min.loss.lambda)
    )
  })
  
  return(
    list(min.loss.lambda=min.loss.lambda,
         results=results_per_lambda,
         loss=loss
        )
  )
}


linear_estimates <- function(data, candidate.mediators, x.cols, treat.col, outcome.col) 
{ # standard linear approach
  
  x_terms = paste(x.cols, collapse=" + ")
  m_terms = paste(candidate.mediators, collapse=" + ")
  
  y_formula <- paste0(outcome.col, " ~ ", treat.col, " + ", m_terms, " + ", x_terms)
  m_forms <- paste(candidate.mediators, "~", treat.col, " + ", x_terms)
  
  fit <- lm(y_formula, data=data)
  NDE_hat <- coef(fit)[treat.col]
  names(NDE_hat) <- NULL
  
  alpha_hats <- sapply(m_forms, function(form) coef(lm(form, data=data))["d"])
  beta_hats <- coef(fit)[candidate.mediators]
  
  NIE_hat <- sum(alpha_hats * beta_hats)
  return(list(NDE_hat=NDE_hat, NIE_hat=NIE_hat, beta_hat=beta_hats, alpha_hat=alpha_hats))
}

adaptive_estimates <- function(data, candidate.mediators, x.cols, treat.col, 
                               outcome.col) 
{ # standard linear approach
  require(glmnet)
  x_terms = paste(x.cols, collapse=" + ")
  m_terms = paste(candidate.mediators, collapse=" + ")
  
  y_formula <- paste0(outcome.col, " ~ ", treat.col, " + ", m_terms, " + ", x_terms)
  
  fit <- lm(y_formula, data=data)
  beta_tildes <- coef(fit)[candidate.mediators]
  pen <- c(0, 1/abs(beta_tildes), rep(0, length(x.cols)))
  # browser()
  y <- data[,outcome.col] %>% as.numeric
  x <- model.matrix(as.formula(paste(y_formula, "+ 0")), data)
  adap.fit <- cv.glmnet(x,y, family="gaussian", penalty.factor=pen)
  
  coefs <- coef(adap.fit, "lambda.min") %>% as.numeric
  coefs <- coefs[2:(length(candidate.mediators)+1)]
  sel.idx <- which(coefs!=0)
  
  y_sel_formula <- paste0(outcome.col, " ~ ", treat.col, " + ", 
                          paste(candidate.mediators[sel.idx], collapse=" + "), 
                          " + ", x_terms)
  
  sel.fit <- lm(y_sel_formula, data=data)
  NDE_hat <- coef(sel.fit)[treat.col]
  names(NDE_hat) <- NULL
  
  m_forms <- paste(candidate.mediators[sel.idx], "~", treat.col, " + ", x_terms)
  
  alpha_hats <- sapply(m_forms, function(form) coef(lm(form, data=data))[treat.col])
  beta_hats <- coef(sel.fit)[candidate.mediators[sel.idx]]
  
  NIE_hat <- sum(alpha_hats * beta_hats)
  return(list(NDE_hat=NDE_hat, NIE_hat=NIE_hat, beta_hat=beta_hats, alpha_hat=alpha_hats,
              sel_M=candidate.mediators[sel.idx]))
}

lasso_estimates <- function(data, candidate.mediators, x.cols, treat.col, 
                            outcome.col) 
{ # standard linear approach
  require(glmnet)
  x_terms = paste(x.cols, collapse=" + ")
  m_terms = paste(candidate.mediators, collapse=" + ")
  
  y_formula <- paste0(outcome.col, " ~ ", treat.col, " + ", m_terms, " + ", x_terms)
  
  # browser()
  y <- data[,outcome.col] %>% as.numeric
  x <- model.matrix(as.formula(paste(y_formula, "+ 0")), data)
  adap.fit <- cv.glmnet(x,y, family="gaussian")
  
  coefs <- coef(adap.fit, "lambda.min") %>% as.numeric
  coefs <- coefs[2:(length(candidate.mediators)+1)]
  sel.idx <- which(coefs!=0)
  
  y_sel_formula <- paste0(outcome.col, " ~ ", treat.col, " + ", 
                          paste(candidate.mediators[sel.idx], collapse=" + "), 
                          " + ", x_terms)
  
  sel.fit <- lm(y_sel_formula, data=data)
  NDE_hat <- coef(sel.fit)[treat.col]
  names(NDE_hat) <- NULL
  
  m_forms <- paste(candidate.mediators[sel.idx], "~", treat.col, " + ", x_terms)
  
  alpha_hats <- sapply(m_forms, function(form) coef(lm(form, data=data))[treat.col])
  beta_hats <- coef(sel.fit)[candidate.mediators[sel.idx]]
  
  NIE_hat <- sum(alpha_hats * beta_hats)
  return(list(NDE_hat=NDE_hat, NIE_hat=NIE_hat, beta_hat=beta_hats, alpha_hat=alpha_hats,
              sel_M=candidate.mediators[sel.idx]))
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


# shorter version
sl_estimates <- function(sel_df, mu.hats,
                         weight.version=NA, lambdas=NA,
                         candidate.mediators, x.cols, treat.col, outcome.col,
                         folds) {
  
  require(dplyr)
  require(SuperLearner)
  
  p = length(candidate.mediators)
  V <- length(folds)
  
  mu.mxis <- mu.hats$mu.mxis
  mu.dx <- mu.hats$mu.dx
  mu.yx <- mu.hats$mu.yx
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  m_0 = sapply(1:p, function(i){
    # estimate the counfounding effect on mediatior i
    mi = candidate.mediators[i]
    Mi <- dplyr::pull(sel_df, mi)
    mu.mxi = mu.mxis[[i]]
    
    col <- Mi - my_sl_predict(mu.mxi)
    return(col)
  })
  stopifnot(ncol(m_0)==p)
  
  d.prob = as.numeric(my_sl_predict(mu.dx))
  
  dc = sel_df[,treat.col] - d.prob
  
  alpha_hats = coef(lm(m_0 ~ dc))[2, ]
  
  y_0 = sel_df[,outcome.col] - my_sl_predict(mu.yx)
  
  yfit = lm(y_0 ~ dc + m_0)
  gamma_hat <- coef(yfit)[2]
  beta_hats <- coef(yfit)[-(1:2)]
  
  NIE_hat <- sum(alpha_hats * beta_hats)
  
  
  return(
    list(alpha_hats = alpha_hats, beta_hats = beta_hats,
         NIE_hat = NIE_hat, NDE_hat = gamma_hat,
         sel_M = candidate.mediators,
         fits=list(y=coef(mu.yx), d=coef(mu.dx), m=lapply(mu.mxis, coef))
    )
  )
  
  
}

train_sl_mu <- function(sel_df, candidate.mediators, 
                        x.cols, treat.col, outcome.col,
                        bin_lib, cont_lib, folds) {
  
  p = length(candidate.mediators)
  V <- length(folds)
  SL.CV.control <- list(V=V, validRows=folds, shuffle=FALSE)
  saveAllFlag <- FALSE
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  
  print("Estimating Mu.mx")
  
  mu.mxis = lapply(1:p, function(i){
    mi = candidate.mediators[i]
    Mi <- dplyr::pull(sel_df, mi)
    
    mu.mxi <- CV.SuperLearner(Y=Mi, X=X, SL.library=cont_lib, family=gaussian(),
                              cvControl=SL.CV.control,
                              saveAll=saveAllFlag,
                              method="method.CC_LS")
    
    return(
      mu.mxi
    )
  })
  
  
  print("Estimating Mu.dx")
  
  mu.dx = try(
    CV.SuperLearner(Y=dplyr::pull(sel_df, treat.col), X=X, 
                    SL.library=bin_lib, family=binomial(),
                    cvControl=SL.CV.control,
                    saveAll=saveAllFlag,
                    method="method.NNloglik")
  )
  
  if(inherits(mu.dx, "try-error")){
    mu.dx <- CV.SuperLearner(Y=dplyr::pull(sel_df, treat.col), X=X, 
                             SL.library=bin_lib, family=binomial(),
                             cvControl=SL.CV.control,
                             saveAll=saveAllFlag,
                             method="method.CC_LS")
  }
  
  
  print("Estimating Mu.yx")
  
  mu.yx = CV.SuperLearner(Y=sel_df[, outcome.col], X=X,
                          SL.library=cont_lib, family=gaussian(),
                          cvControl=SL.CV.control,
                          saveAll=saveAllFlag,
                          method="method.CC_LS")
  
  return(
    list(mu.mxis=mu.mxis,
         mu.dx=mu.dx,
         mu.yx=mu.yx)
  )
  
}

true_sl_mu <- function(sel_df, candidate.mediators, 
                       x.cols, treat.col, outcome.col,
                       psi_d, psi_m, psi_y,
                       true_alpha, true_beta, true_gamma) {
  
  p = length(candidate.mediators)
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  mu.dx = apply(X, 1, psi_d)
  psi_m_xs = apply(X, 1, psi_m)
  psi_y_xs = apply(X, 1, psi_y)
  stopifnot(is.numeric(mu.dx),
            is.numeric(psi_m_xs),
            is.numeric(psi_y_xs))
  
  mu.mxis = lapply(1:p, function(i){
    mi = candidate.mediators[i]
    Mi <- dplyr::pull(sel_df, mi)
    
    mu.mxi <- true_alpha[i] * mu.dx + psi_m_xs
    
    return(mu.mxi)
  })

  mu.yx = (true_gamma + sum(true_alpha*true_beta))*mu.dx +
          sum(true_beta)*psi_m_xs +
          psi_y_xs
  
  return(
    list(mu.mxis=mu.mxis,
         mu.dx=mu.dx,
         mu.yx=mu.yx)
  )
  
}
