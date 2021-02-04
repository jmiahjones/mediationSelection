
muyd.err <- foreach(sim.idx=1:num_simulations, .combine=rbind) %do% {
  sim <- simulations[[sim.idx]]
  
  muy.err <- mu.hats[[sim.idx]]$mu.yx - sim$confounding_y
  mud.err <- mu.hats[[sim.idx]]$mu.dx - sim$propensity
  # mum.err <- mu.hats[[sim.idx]]$mu.mxis[[1]] - sim$confounding_m[,1]
  
  muy.norm <- sqrt(mean(muy.err^2))
  mud.norm <- sqrt(mean(mud.err^2))
  
  c(muy.norm, mud.norm)
}


all.mum.err <- foreach(sim.idx=1:num_simulations, .combine=rbind) %:%
  foreach(j=1:p, .combine=cbind) %do%
  {
    sim <- simulations[[sim.idx]]
    
    mum.err <- mu.hats[[sim.idx]]$mu.mxis[[j]] - sim$confounding_m[,j]
    (mum.err^2) %>% mean %>% sqrt
  }

all.err <- cbind(muyd.err, all.mum.err)
dim(all.err)

emp.err <- colMeans(all.err) %>% matrix(ncol=1)
rownames(emp.err) <- c("Y", "D", paste0("M", 1:p))

n^(.25)*emp.err
emp.err[1]

foreach(sim.idx=1:num_simulations, .combine=rbind) %do%{
  mu.hats[[sim.idx]]$mu.mxi.coefs[[1]] %>% colMeans
} %>% summary


cl <- makeCluster(5, "FORK")
doParallel::registerDoParallel(cl)
perf <- foreach(sim.idx=1:num_simulations, .combine=rbind) %dopar%{
  sim <- simulations[[sim.idx]]
  rob_sys <- robinsonize_system(sim$d, 
                                select(sim, candidate.mediators), 
                                sim$y,
                                mu.hats[[sim.idx]])
  dc <- rob_sys$dc
  m_0 <- rob_sys$m_0
  y_0 <- rob_sys$y_0
  
  fit <- lm(y_0 ~ dc + m_0)
  se <- sandwich::vcovHC(fit)[2,2] %>% sqrt
  nde_hat <- coef(fit)[2] %>% unname
  fitted <- fitted(fit)
  resid <- resid(fit)
  ci_boot <- sapply(1:1000, function(b){
    g <- 2*rbinom(n, size=1, prob=.5) - 1
    resid_star <- g*resid
    y_0_star <- fitted + resid_star
    qr.coef(fit$qr, y_0_star)[2]
  }) %>% quantile(probs=c(.025, .975))
  cover_boot <- 1*(ci_boot[1] <= NDE & ci_boot[2] >= NDE) %>% unname
  ci <- nde_hat + c(-1,1)*1.96*se
  cover <- 1*(ci[1] <= NDE & ci[2] >= NDE) %>% unname
  c(nde_hat, se, cover_boot, cover, ci_boot, ci)
}
stopCluster(cl)
stopImplicitCluster()

foreach(res=inference.results[[1]]$prd.results,.combine=mean) %do%{
  sum(res$beta_tildes[1:3]*res$alpha_tildes[1:3])
}

lm.mu.err <- foreach(sim.idx=1:num_simulations, .combine=c) %do%
{
  sim <- simulations[[sim.idx]]
  
  fit1 <- lm(m.1 ~ I(x.1^2) + x.2,
              data=sim,
              subset = split_folds[[1]])

  fit2 <- lm(m.1 ~ I(x.1^2) + x.2,
              data=sim,
              subset = split_folds[[2]])
  
  pred <- rep(0, n)
  pred[split_folds[[1]]] <- predict(fit2, sim[split_folds[[1]],])
  pred[split_folds[[2]]] <- predict(fit1, sim[split_folds[[2]],])
  
  mum.err <- pred - sim$confounding_m[,1]
  (mum.err^2) %>% mean %>% sqrt
}
mean(lm.mu.err)*n^(.25)


# gam_mu <- foreach(sim=simulations[1:100]) %dopar%{
#   train_sl_mu(
#     sim, candidate.mediators=candidate.mediators[1],
#     x.cols=x.cols,
#     treat.col="d", outcome.col="y", 
#     bin_lib=bin_lib, cont_lib=create_gam$names, 
#     folds=split_folds,
#     estD=F, 
#     estM=T,
#     estY=F
#   )
# }


library(dplyr)
library(doParallel)
library(foreach)
bs <- "ts"
m <- 2
k <- 2
library(mgcv)
# library(randomForest)
cl <- makeCluster(10, "FORK")
doParallel::registerDoParallel(cl)
ns <- c(500, 1000, 1500, 2000)
split_list <- c(2,4,8)
splits <- 2
n=2000
psi_m <- function(x) { 
  x1 <- x[1]
  x2 <- x[2]
  1*x1^2 + x2
}

gam.mu.err <- foreach(sim.idx=1:100, .combine=rbind) %:%
  # foreach(n=ns, .combine=c) %dopar%
  foreach(splits=split_list, .combine=c) %dopar%
{
  sim <- replicate(5, rnorm(n, sd=.5))
  colnames(sim) <- paste0("x.", 1:5)
  mu <- sim %>% apply(1, psi_m)
  sim <- as.data.frame(sim)
  sim$m.new <- mu + rnorm(n)
  
  split_ids <- cut(1:n, splits) %>% as.integer
  pred <- rep(0, n)
  for(val.idx in 1:splits){
    val.rows <- which(split_ids == val.idx)
    fit <- gam(m.new ~ s(x.1, bs=bs, m=m, k=k) + 
                 s(x.2, bs=bs, m=m, k=k) +
                 s(x.3, bs=bs, m=m, k=k) + 
                 s(x.4, bs=bs, m=m, k=k) +
                 s(x.5, bs=bs, m=m, k=k),
               data=sim[-val.rows,],
               family=gaussian())
    pred[val.rows] <- predict(fit, sim[val.rows,])
  }
  
  mum.err <- pred - mu
  (mum.err^2) %>% mean %>% sqrt
}
stopCluster(cl)
plot(split_list, colMeans(gam.mu.err))
# plot(ns, colMeans(gam.mu.err)*n^(.25))


fiterr <- lm(I(colMeans(gam.mu.err)) ~ I(ns^(-1/4)))
lines(ns, fiterr$fitted.values)
# mean(gam.mu.err)*n^(.25)
# coef(fiterr)




library(doParallel)
library(mgcv)
library(randomForest)
cl <- makeCluster(20, "FORK")
doParallel::registerDoParallel(cl)

rf.mu.err <- foreach(sim.idx=1:num_simulations, .combine=c) %dopar%
{
  sim <- simulations[[sim.idx]]
  
  # fit1 <- lm(m.1 ~ I(x.1^2) + x.2,
  #             data=sim,
  #             subset = split_folds[[1]])
  # 
  # fit2 <- lm(m.1 ~ I(x.1^2) + x.2,
  #             data=sim,
  #             subset = split_folds[[2]])
  
  fit <- randomForest(m.1 ~ x.1 + x.2 + x.3 + x.4 + x.5,
                      data=sim,
                      mtry=1,
                      nodesize=5L,
                      ntree=1000L)
  
  # use OOB predictions
  pred <- fit$predicted
  
  # pred <- rep(0, n)
  # pred[split_folds[[1]]] <- predict(fit2, sim[split_folds[[1]],])
  # pred[split_folds[[2]]] <- predict(fit1, sim[split_folds[[2]],])
  
  mum.err <- pred - sim$confounding_m[,1]
  (mum.err^2) %>% mean %>% sqrt
}
stopCluster(cl)
mean(rf.mu.err)*n^(.25)



library(doParallel)
library(kernlab)
cl <- makeCluster(20, "FORK")
doParallel::registerDoParallel(cl)
ksvm_mu <- foreach(sim=simulations[1:100], .combine=c) %dopar%{
  pred <- train_sl_mu(
    sim, candidate.mediators=candidate.mediators[1],
    x.cols=x.cols,
    treat.col="d", outcome.col=candidate.mediators[1],
    bin_lib=bin_lib, cont_lib="SL.ksvm",
    folds=split_folds,
    estD=F,
    estM=F,
    estY=T
  )$mu.yx
  mum.err <- pred - sim$confounding_m[,1]
  (mum.err^2) %>% mean %>% sqrt
}
stopCluster(cl)
mean(ksvm_mu)*n^(.25)




# check custom mgcv function
library(doParallel)
library(mgcv)
create_mgcv <- create.Learner("SL.mgcv", tune=list(
  m=2,
  k=c(1,2,4,6),
  bs=c("ts")
))
cl <- makeCluster(20, "FORK")
doParallel::registerDoParallel(cl)
mgcv_mu <- foreach(sim=simulations) %dopar%{
  sl <- train_sl_mu(
    sim, candidate.mediators=candidate.mediators[1],
    x.cols=x.cols,
    treat.col="d", outcome.col=candidate.mediators[1],
    bin_lib=bin_lib, cont_lib=create_mgcv$names,
    folds=split_folds,
    estD=F,
    estM=F,
    estY=T
  )
  pred <- sl$mu.yx
  mum.err <- pred - sim$confounding_m[,1]
  l2 <- (mum.err^2) %>% mean %>% sqrt
  # discrete_sl <- do.call(c, sl$whichDiscreteSL)
  # c(l2, as.numeric(substr(discrete_sl, 9,10)))
  list(err=l2, coef=sl$mu.yx.coef)
}
stopCluster(cl)
mean(sapply(mgcv_mu, function(x) x$err))*n^(.25)

coefs <- foreach(m=mgcv_mu, .combine=rbind) %do% {m$coef}
summary(coefs)


# check bias
prd.bias <- foreach(sim.idx=1:num_simulations, .combine=rbind) %do% {
  res <- inference.results[[1]]$prd.results[[sim.idx]]
  
  c(
    NIE=res$NIE_hat - NIE,
    NDE=res$NDE_hat - NDE,
    alpha=res$alpha_tildes - alphas,
    beta=res$beta_tildes - betas
    # mean(res$beta_hats - betas[as.numeric(substr(res$sel_M, 3,4))])
  )
}
prd.bias %>%
  # abs %>%
  colMeans



mu.hats[[1]]$mu.mxi.coefs[[1]]
mu.hats[[1]]$mu.yx.coef
mu.hats[[1]]$mu.dx.coef

mu.err[,2] %>% mean

mu.true <- foreach(sim.idx=1:num_simulations, .combine=rbind) %do% {
  sim <- simulations[[sim.idx]]
  
  muy <- sim$confounding_y
  mud <- sim$propensity
  mum <- sim$confounding_m[,1]
  cbind(muy, mud, mum)
}


sim.idx=466
sim <- simulations[[sim.idx]]
muy.err <- mu.hats[[sim.idx]]$mu.yx - sim$confounding_y
mud.err <- mu.hats[[sim.idx]]$mu.dx - sim$propensity
mum.err <- mu.hats[[sim.idx]]$mu.mxis[[1]] - sim$confounding_m[,1]
mu.err <- cbind(muy.err, mud.err, mum.err)
sqrt(colMeans(mu.err^2))
sqrt(colMeans(mu.err^2))*n^(.25)
mc <- foreach(j=1:p, .combine=cbind) %do% {
  mu.mj <- mu.hats[[sim.idx]]$mu.mxis[[j]]
  mc <- pull(sim, candidate.mediators[j]) - mu.mj
}

yc <- sim$y - mu.hats[[sim.idx]]$mu.yx
dc <- sim$d - mu.hats[[sim.idx]]$mu.dx
lm(yc ~ dc + mc)
coef(lm(yc ~ dc + mc))
prd.results[[sim.idx]]$alpha_tildes
prd.results[[sim.idx]]$alpha_tildes - alphas
alphas
prd.results[[sim.idx]]$beta_tildes - betas
max(abs(prd.results[[sim.idx]]$beta_tildes - betas))
max(abs(prd.results[[sim.idx]]$alpha_tildes - alphas))
mu.hats[[sim.idx]]$mu.mxi.coefs[[1]]
sim <- simulations[[sim.idx]]

# check mu errors
# compare NDEs
# nde_compare <- foreach(sim.idx=1:num_simulations, .combine=rbind) %do% {
merr <- foreach(sim.idx=1:num_simulations, .combine=rbind) %do% {
  sim <- simulations[[sim.idx]]
  # true_mu <- true_sl_mu(sim, candidate.mediators=candidate.mediators)
  mctrue <- foreach(j=1:p, .combine=cbind) %do% {
    mu.mj <- sim$confounding_m[,j]
    dplyr::pull(sim, candidate.mediators[j]) - mu.mj
  }
  yctrue <- sim$y - sim$confounding_y
  dctrue <- sim$d - sim$propensity
  
  
  mc <- foreach(j=1:p, .combine=cbind) %do% {
    mu.mj <- mu.hats[[sim.idx]]$mu.mxis[[j]]
    mc <- pull(sim, candidate.mediators[j]) - mu.mj
  }
  yc <- sim$y - mu.hats[[sim.idx]]$mu.yx
  dc <- sim$d - mu.hats[[sim.idx]]$mu.dx
  alpha_hats <- coef(lm(mc ~ dc - 1))[2,]
  alpha_hats2 <- mu.hats[[sim.idx]]$mu.dx
  
  # c(
  #   coef(lm(yc ~ dc + mc + mu.hats[[sim.idx]]$mu.dx))
  #   coef(lm(yc ~ dc + mctrue))[2],
  #   coef(lm(yc ~ dctrue + mc))[2],
  #   coef(lm(yc ~ dctrue + mctrue))[2],
  #   coef(lm(yctrue ~ dc + mc))[2],
  #   coef(lm(yctrue ~ dctrue + mctrue))[2]
  # )
}
summary(nde_compare - NDE)


m.err <- foreach(i=1:num_simulations, .combine=rbind) %:% foreach(j=1:p, .combine=c) %do%{
  sqrt(mean((mu.hats[[i]]$mu.mxis[[j]] - simulations[[i]]$confounding_m[,j])^2))
}
n^(1/4)*colMeans(m.err)






ols_nie <- foreach(i=1:num_simulations, .combine=c) %do%{
  alpha_til <- inference.results[[1]]$prd.results[[i]]$alpha_tildes
  beta_til <- inference.results[[1]]$prd.results[[i]]$beta_tildes
  sum(alpha_til*beta_til)
}
mean(ols_nie - NIE)

