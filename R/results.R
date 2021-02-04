
library(dplyr)
library(foreach)
library(doParallel)

plot_dir <- "./plots/"
all_res_save <- "./cache/all_res.RData"
weight_names <- c("prd", "mix", "adp")
method_names <- c("das", "minnier")

# n <- 2000
filenames <- dir("./cache/", pattern = "*selection.RData")
filenames <- paste0("./cache/", filenames) #[grep(paste0(n, "*"), filenames)]
#TODO Remove
filenames <- filenames[1:24]

naive_delta_inf <- function(dc, m_0, y_0, sel_M_idxs, ret_est=F){
  m_sel <- m_0[,sel_M_idxs]
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
  

# cl_num <- min(length(weight_names)*length(method_names), 10)
cl_num <- min(length(weight_names), 10)
cl <- makeCluster(cl_num, "FORK")
registerDoParallel(cl)
# registerDoSEQ()

all_res <- foreach(file=filenames, file.idx=seq_along(filenames), .combine=rbind) %do% {
  print(paste0("Starting ", file.idx, " of ", length(filenames), "..."))
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
  
  # naive_df <- foreach(sim.idx=1L:num_simulations, .combine=rbind) %do% {
  #   # print(paste0("Beginning weight: ", weight, "."))
  #   if(sim.idx %% 100 == 1)
  #     print(paste0("Beginning ", sim.idx, " of ", num_simulations, "."))
  #   
  #   # selection data
  #   sel_df <- simulations[[sim.idx]]
  #   X <- sel_df %>% dplyr::select(all_of(x.cols))
  #   D <- clean_D(sel_df, "d")
  #   M <- sel_df %>% dplyr::select(all_of(candidate.mediators))
  #   Y <- sel_df %>% dplyr::pull(y)
  #   
  #   p <- ncol(M)
  #   
  #   rob_system <- robinsonize_system(D, M, Y, mu.hats[[sim.idx]])
  #   dc <- rob_system$dc
  #   m_0 <- rob_system$m_0
  #   y_0 <- rob_system$y_0
  #   
  #   coverages <- foreach(weight=weight_names, .combine=rbind) %dopar%
  #   {
  #     # print(sim.idx)
  #     method_results <- inference.results[[1]]
  #     these_results <- eval(parse(text=paste0("method_results$", weight, ".results")))
  #     sel_M <- these_results[[sim.idx]]$sel_M
  #     
  #     # TODO: put in naive bootstrap
  #     # naive_cis <- foreach(b=1:1000L, .combine=rbind) %do%{
  #     #   set.seed(45701+b)
  #     #   boot.idx <- sample.int(n, replace=T)
  #     #   yb <- y_0[boot.idx]
  #     #   mb <- m_0[boot.idx,]
  #     #   db <- dc[boot.idx]
  #     #
  #     #   yfit_b <- lm(yb ~ db + mb)
  #     #   beta_b <- coef(yfit_b)[-(1:2)]
  #     #   NDE_b <- coef(yfit_b)[2]
  #     #   mfit_b <- lm(mb ~ db)
  #     #   if(p>1){
  #     #     alpha_b <- coef(mfit_b)[2,]
  #     #   } else {
  #     #     alpha_b <- coef(mfit_b)[2]
  #     #   }
  #     #   NIE_b <- sum(alpha_b * beta_b)
  #     #   c(NDE_b, NIE_b)
  #     # } %>% apply(2, quantile, probs=c(0.025, 0.975))
  # 
  #     M_idxs <- which(candidate.mediators %in% sel_M)
  #     naive_cis <- naive_delta_inf(dc, m_0, y_0, M_idxs)
  # 
  #     NDE_coverage <- naive_cis[1,1] <= NDE & NDE <= naive_cis[2,1]
  #     NIE_coverage <- naive_cis[1,2] <= NIE & NIE <= naive_cis[2,2]
  # 
  #     c(NDE_coverage, NIE_coverage)
  #   }
  # 
  #   data.frame(
  #     stat="mean",
  #     val=coverages,
  #     # col=c("coverage_NDE", "coverage_NIE"),
  #     col=rep(c("coverage_NDE", "coverage_NIE"), times=3L),
  #     n=n,
  #     nsim=num_simulations,
  #     scenario=abbv_scn,
  #     coefs=coef_setting,
  #     weight_name=rep(weight_names, each=2L), # weight_name=weight,
  #     weight_gamma=weight_gam_str, # change to str for multiple gam
  #     boot_method="naive",
  #     SL=use_sl,
  #     suffix=suffix_arg
  #   )
  # }
  # 
  # naive_df %>% group_by(-val) %>% summarize(val=mean(val)) %>% 
  #   ungroup %>% rbind(res_df) 
  #%>% filter(col=="coverage_NDE", boot_method=="naive") %>% head
  return(res_df)
}

stopImplicitCluster()

save.image(file=all_res_save)

# all_res %>% filter(stat=="mean", col=="coverage_NDE", weight_name=="prd", boot_method=="minnier") %>% View
# all_res %>% filter(stat=="mean", col=="NDE_err", weight_name=="prd") %>% View
# all_res %>% filter(stat=="mean", col=="coverage_NIE", weight_name=="prd", boot_method=="minnier") %>% View
# all_res %>% filter(stat=="mean", col=="NIE_err", weight_name=="prd") %>% View
# res_df %>% View
all_res %>% filter(stat=="mean", col=="is_missed", weight_name=="adp",
                   boot_method=="minnier") %>% View

# # Root-n bias
# library(ggplot2)
# plot_list <- vector("list", 6L)
# plot_list[[1]] <- all_res %>% 
#   filter(stat=="mean", col %in% c("NDE_err", "NIE_err"), 
#          weight_name=="prd",
#          boot_method=="das") %>% 
#   mutate(val=sqrt(n)*abs(val), col=substr(col, 1,3),
#          weight_gamma=factor(weight_gamma)
#   ) %>% 
#   ggplot(aes(x=scenario, y=val, color=coefs, shape=SL)) + 
#   geom_point(size=4, alpha=0.7, position=position_jitter(width=.25, height=0)) +
#   labs(title="Root-n Absolute Bias due to SuperLearner", 
#        subtitle="Product Weights",
#        x="Psi Function Scenario", y="Root-n Absolute Bias", 
#        color="Coef Size", shape="SL Used T/F") +
#   facet_grid(cols=vars(col))
# 
# plot_list[[2]] <- all_res %>% 
#   filter(stat=="mean", col %in% c("NDE_err", "NIE_err"), 
#          weight_name=="mix",
#          boot_method=="das") %>% 
#   mutate(val=sqrt(n)*abs(val), col=substr(col, 1,3),
#          weight_gamma=factor(weight_gamma)
#   ) %>% 
#   ggplot(aes(x=scenario, y=val, color=coefs, shape=SL)) + 
#   geom_point(size=4, alpha=0.7, position=position_jitter(width=.25, height=0)) +
#   labs(title="Root-n Absolute Bias due to SuperLearner", 
#        subtitle="Mix Weights",
#        x="Psi Function Scenario", y="Root-n Absolute Bias", 
#        color="Coef Size", shape="SL Used T/F") +
#   facet_grid(cols=vars(col))
# 
# 
# plot_list[[3]] <- all_res %>% 
#   filter(stat=="mean", col %in% c("NDE_err", "NIE_err"), 
#          weight_name=="adp",
#          boot_method=="das") %>% 
#   mutate(val=sqrt(n)*abs(val), col=substr(col, 1,3),
#          weight_gamma=factor(weight_gamma)
#   ) %>% 
#   ggplot(aes(x=scenario, y=val, color=coefs, shape=SL)) + 
#   geom_point(size=4, alpha=0.5, position=position_jitter(width=.25, height=0)) +
#   labs(title="Root-n Absolute Bias due to SuperLearner", 
#        subtitle="Adaptive Lasso Weights",
#        x="Psi Function Scenario", y="Root-n Absolute Bias", 
#        color="Coef Size", shape="SL Used T/F") +
#   facet_grid(cols=vars(col))
# 
# 
# # Coverage Rates:
# plot_list[[4]] <- all_res %>% 
#   filter(stat=="mean", col %in% c("coverage_NDE", "coverage_NIE"), 
#          weight_name=="prd"
#   ) %>%
#   mutate(col=substr(col, 10,13)
#          # weight_gamma=factor(weight_gamma)
#   ) %>% 
#   ggplot(aes(x=scenario, y=val, color=boot_method, shape=SL)) + 
#   geom_hline(yintercept=0.95) +
#   geom_point(size=3, alpha=0.7, position=position_jitter(width=.25, height=0)) +
#   labs(title="Coverage Rates for Das, Minnier, and Naive Bootstraps", 
#        subtitle="Product Weights",
#        x="Psi Function Scenario", y="Coverage Rate", 
#        color="Boot Method", shape="SL Used T/F") +
#   facet_grid(cols=vars(col), rows=vars(coefs))
# 
# plot_list[[5]] <- all_res %>% 
#   filter(stat=="mean", col %in% c("coverage_NDE", "coverage_NIE"), 
#          weight_name=="mix"
#   ) %>%
#   mutate(col=substr(col, 10,13)
#          # weight_gamma=factor(weight_gamma)
#   ) %>% 
#   ggplot(aes(x=scenario, y=val, color=boot_method, shape=SL)) + 
#   geom_hline(yintercept=0.95) +
#   geom_point(size=3, alpha=0.7, position=position_jitter(width=.25, height=0)) +
#   labs(title="Coverage Rates for Das, Minnier, and Naive Bootstraps", 
#        subtitle="Mixture Weights",
#        x="Psi Function Scenario", y="Coverage Rate", 
#        color="Boot Method", shape="SL Used T/F") +
#   facet_grid(cols=vars(col), rows=vars(coefs))
# 
# plot_list[[6]] <- all_res %>% 
#   filter(stat=="mean", col %in% c("coverage_NDE", "coverage_NIE"), 
#          weight_name=="adp"
#   ) %>%
#   mutate(col=substr(col, 10,13)
#          # weight_gamma=factor(weight_gamma)
#   ) %>% 
#   ggplot(aes(x=scenario, y=val, color=boot_method, shape=SL)) + 
#   geom_hline(yintercept=0.95) +
#   geom_point(size=3, alpha=0.7, position=position_jitter(width=.25, height=0)) +
#   labs(title="Coverage Rates for Das, Minnier, and Naive Bootstraps", 
#        subtitle="Adaptive Lasso Weights",
#        x="Psi Function Scenario", y="Coverage Rate", 
#        color="Boot Method", shape="SL Used T/F") +
#   facet_grid(cols=vars(col), rows=vars(coefs))
# 
# 
# pdf(paste0(plot_dir, "bias.pdf"))
# for(i in 1:3){
#   print(plot_list[[i]])
# }
# dev.off()
# 
# pdf(paste0(plot_dir, "coverage.pdf"))
# for(i in 4:6){
#   print(plot_list[[i]])
# }
# dev.off()

all_res %>% 
  filter(stat=="median", boot_method=="minnier", col=="num_missed") %>% 
  mutate(n=factor(n), coefs=factor(coefs)) %>% 
  ggplot(aes(x=n, y=val, color=weight_name, shape=SL)) +
  geom_point(position = position_dodge(width=.2), size=4, alpha=0.7) +
  facet_grid(rows=vars(coefs), cols=vars(scenario)) +
  labs(title="Median Number of Mediators Missed")

all_res %>% 
  filter(stat=="mean", boot_method=="minnier", col=="is_missed") %>% 
  mutate(n=factor(n), coefs=factor(coefs)) %>% 
  ggplot(aes(x=n, y=val, color=weight_name, shape=SL)) +
  geom_point(position = position_dodge(width=.4), size=4, alpha=0.7) +
  facet_grid(rows=vars(coefs), cols=vars(scenario)) +
  labs(title="Percent of Scenarios where any Mediator was Excluded")

#### Create Variable Selection Performance Table ####
all_res %>%
  filter((stat=="mean" & boot_method=="minnier" & col=="is_missed") |
         (stat=="mean" & boot_method=="minnier" & col=="num_noise")
         ) %>%
  # filter(n==2000, coefs=="large", SL=="T") %>%
  mutate(scenario=toupper(scenario),
         val=if_else(col == "is_missed", 1-val, val),
         col=if_else(col == "is_missed", "PC", "MN")) %>% 
  tidyr::pivot_wider(id_cols=c(coefs, n, weight_name, SL), 
                     names_from=c(scenario, col),
                     names_glue="{scenario}_{col}",
                     values_from=val) %>% 
  arrange(coefs, n, desc(weight_name), SL) %>% 
  knitr::kable(format="latex", digits=2L)
  

library(ggplot2)
small_plot_list <- vector("list", 4)

small_plot_list[[1]] <- all_res %>% 
  filter(stat=="mean", col %in% c("NDE_err", "NIE_err"), 
         boot_method=="minnier",
         coefs=="large") %>% 
  mutate(val=sqrt(n)*abs(val), col=substr(col, 1,3)) %>% 
  mutate(n=factor(n)) %>% 
  ggplot(aes(x=n, y=val, color=weight_name, shape=SL)) + 
  geom_point(size=4, alpha=0.7, position=position_jitter(width=.05, height=0)) +
  labs(title="Root-n Absolute Bias due to SuperLearner", 
       subtitle="Large Setting",
       x="n", y="Root-n Absolute Bias", 
       color="Weight", shape="SL Used T/F") +
  facet_grid(rows=vars(col), cols=vars(scenario))

small_plot_list[[2]] <- all_res %>% 
  filter(stat=="mean", col %in% c("NDE_err", "NIE_err"), 
         boot_method=="minnier",
         coefs=="small") %>% 
  mutate(val=sqrt(n)*abs(val), col=substr(col, 1,3)
  ) %>% 
  mutate(n=factor(n)) %>% 
  ggplot(aes(x=n, y=val, color=weight_name, shape=SL)) + 
  geom_point(size=4, alpha=0.7, position=position_jitter(width=.05, height=0)) +
  labs(title="Root-n Absolute Bias due to SuperLearner", 
       subtitle="Small Setting",
       x="n", y="Root-n Absolute Bias", 
       color="Weight", shape="SL Used T/F") +
  facet_grid(rows=vars(col), cols=vars(scenario))


small_plot_list[[3]] <- all_res %>% 
  filter(stat=="mean", 
         col %in% c("coverage_NDE", "coverage_NIE"), 
         boot_method=="minnier",
         coefs=="large") %>% 
  mutate(col=substr(col, 10,13)) %>%
  mutate(n=factor(n)) %>% 
  ggplot(aes(x=n, y=val, color=weight_name, shape=SL)) + 
  geom_point(size=4, alpha=0.7, position=position_jitter(width=.05, height=0)) +
  geom_hline(yintercept=.95) +
  labs(title="Coverage Rates for Minnier Bootstrap", 
       subtitle="Large Setting",
       x="n", y="Coverage", 
       color="Weight", shape="SL Used T/F") +
  facet_grid(rows=vars(col), cols=vars(scenario))

small_plot_list[[4]] <- all_res %>% 
  filter(stat=="mean", 
         col %in% c("coverage_NDE", "coverage_NIE"), 
         boot_method=="minnier",
         coefs=="small") %>% 
  mutate(col=substr(col, 10,13)) %>%
  mutate(n=factor(n)) %>% 
  ggplot(aes(x=n, y=val, color=weight_name, shape=SL)) + 
  geom_point(size=4, alpha=0.7, position=position_jitter(width=.05, height=0)) +
  geom_hline(yintercept=.95) +
  labs(title="Coverage Rates for Minnier Bootstrap", 
       subtitle="Small Setting",
       x="n", y="Coverage", 
       color="Weight", shape="SL Used T/F") +
  facet_grid(rows=vars(col), cols=vars(scenario))


pdf(paste0(plot_dir, "bias.pdf"))
for(i in 1:2){
  print(small_plot_list[[i]])
}
dev.off()

pdf(paste0(plot_dir, "coverage.pdf"))
for(i in 3:4){
  print(small_plot_list[[i]])
}
dev.off()