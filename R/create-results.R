args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  stop("No commandline arguments found!")
}

stopifnot(length(args) == 3)
proc_num <- as.numeric(args[1])
total_num <- as.numeric(args[2])
no_cores <- as.numeric(args[3])


library(dplyr)
library(foreach)
library(doParallel)

tmp.dir <- "./cache/temp/"
if(!dir.exists(tmp.dir))
  dir.create(tmp.dir)

temp_save <- paste0(tmp.dir, "cache-", proc_num, ".RData")

weight_names <- c("prd", "mix", "adp")
method_names <- c("das", "minnier")

filenames <- dir("./cache/", pattern = "*selection.RData")
filenames <- paste0("./cache/", filenames)

num_to_do <- length(filenames)
file_idxs <- which(
  as.integer(
    cut(seq.int(num_to_do), breaks=total_num)
  ) == proc_num)


naive_delta_inf <- function(dc, m_0, y_0, sel_M_idxs, ret_est=F){
  m_sel <- as.matrix(m_0[,sel_M_idxs])
  p_sel <- ncol(m_sel)
  even_idx_sel <- 2*(1:p_sel)
  sel_naive_fit_y <- lm(y_0 ~ dc + m_sel)
  sel_naive_fit_m <- lm(m_sel ~ dc)
  nde_se <- sqrt(vcov(sel_naive_fit_y)[1,1])
  beta_hat <- coef(sel_naive_fit_y)[-(1:2)]
  alpha_hat <- coef(sel_naive_fit_m)[2,]
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
  
  ci_mat <- cbind(nde_ci, nie_ci)
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


cl <- makeCluster(no_cores, "FORK")
registerDoParallel(cl)
# registerDoSEQ()


temp_res <- foreach(file.idx=file_idxs, .combine=rbind) %do% {
  print(paste0("Starting ", file.idx, " of ", length(filenames), "..."))
  file <- filenames[file.idx]
  load(file)
  source("./R/simulation-helpers.R")
  
  # method.idx: method_names
  res_df <- foreach(method.idx=2, .combine=rbind) %:%
    foreach(weight=weight_names, .combine=rbind) %dopar%
    {
      boot_name=method_names[method.idx]
      method_results <- inference.results[[method.idx]]
      these_results <- eval(parse(text=paste0("method_results$", weight, ".results")))
      res_df <- result_dataframe(these_results)
      
      summ_df <- res_df %>%
        summarize(across(everything(),
                         list(min=min, median=median, mean=mean, max=max),
                         .names="{.fn}_{.col}")
        ) %>%
        tidyr::gather(key="stat", value="val", everything())
      
      stat_col <- stringr::str_split(summ_df$stat, "\\_", n=2, simplify=T)
      summ_df <- summ_df %>% mutate(stat=stat_col[,1], col=stat_col[,2]) %>%
        mutate(
          n=n,
          nsim=num_simulations,
          scenario=abbv_scn,
          coefs=coef_setting,
          weight_name=weight,
          weight_gamma=weight_gam_str, # change to str for multiple gam
          boot_method=boot_name,
          SL=use_sl,
          suffix=suffix_arg
        )
      return(summ_df)
    }
  
  naive_df <- foreach(sim.idx=1L:num_simulations, .combine=rbind) %dopar% {
    # print(paste0("Beginning weight: ", weight, "."))
    if(sim.idx %% 100 == 1)
      print(paste0("Beginning ", sim.idx, " of ", num_simulations, "."))

    # selection data
    sel_df <- simulations[[sim.idx]]
    X <- sel_df %>% dplyr::select(all_of(x.cols))
    D <- clean_D(sel_df, "d")
    M <- sel_df %>% dplyr::select(all_of(candidate.mediators))
    Y <- sel_df %>% dplyr::pull(y)

    p <- ncol(M)

    rob_system <- robinsonize_system(D, M, Y, mu.hats[[sim.idx]])
    dc <- rob_system$dc
    m_0 <- rob_system$m_0
    y_0 <- rob_system$y_0

    coverages <- foreach(weight=weight_names, .combine=c) %do%
    {
      # print(sim.idx)
      method_results <- inference.results[[1]]
      these_results <- eval(parse(text=paste0("method_results$", weight, ".results")))
      sel_M <- these_results[[sim.idx]]$sel_M

      M_idxs <- which(candidate.mediators %in% sel_M)
      naive_cis <- naive_delta_inf(dc, m_0, y_0, M_idxs)

      NDE_coverage <- naive_cis[1,1] <= NDE & NDE <= naive_cis[2,1]
      NIE_coverage <- naive_cis[1,2] <= NIE & NIE <= naive_cis[2,2]

      c(NDE_coverage, NIE_coverage)
    }

    data.frame(
      stat="mean",
      val=coverages,
      # col=c("coverage_NDE", "coverage_NIE"),
      col=rep(c("coverage_NDE", "coverage_NIE"), times=3L),
      n=n,
      nsim=num_simulations,
      scenario=abbv_scn,
      coefs=coef_setting,
      weight_name=rep(weight_names, each=2L), # weight_name=weight,
      weight_gamma=weight_gam_str, # change to str for multiple gam
      boot_method="naive",
      SL=use_sl,
      suffix=suffix_arg
    )
  }

  naive_df %>% group_by_at(vars(!val)) %>% 
    summarize(val=mean(val)) %>% ungroup %>%
    rbind(res_df) %>% 
    return
  # %>% filter(col=="coverage_NDE", boot_method=="naive") %>% head
  # return(res_df)
}

stopImplicitCluster()

save(temp_res, file=temp_save)
