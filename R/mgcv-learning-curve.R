num_cores <- 10

library(dplyr)
library(doParallel)
library(foreach)
bs <- "ts"
m <- 2
k <- 2
library(mgcv)

cl <- makeCluster(num_cores, "FORK")
doParallel::registerDoParallel(cl)
ns <- c(500, 1000, 1500, 2000, 3000, 5000)
split_list <- c(2,4,8)
splits <- 2
n=2000
psi_m <- function(x) { 
  x1 <- x[1]
  x2 <- x[2]
  1*x1^2 + x2
}


gam.mu.err <- foreach(sim.idx=1:100, .combine=rbind) %:%
  foreach(n=ns, .combine=c) %dopar%
  # foreach(splits=split_list, .combine=c) %dopar%
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
# plot(split_list, colMeans(gam.mu.err), xlab="Number of Splits", ylab="L2 Mu Error")
plot(ns/splits, colMeans(gam.mu.err)*n^(.25), xlab="Training Sample Size", 
     ylab="n^(1/4)*L2 Mu Error",
     main="MGCV Learning Curve", sub=paste0("Splits: ", splits))