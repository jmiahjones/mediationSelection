library(ggplot2)
library(qs)
library(dplyr)

suffix <- "ml0-9-1.qs"
files <- dir("./results/", "*.qs", full.names = T)
files <- files[grep(suffix, files)]

all_results <- lapply(files, qread, nthreads=2) %>% bind_rows
all_results <- all_results %>% mutate(
  NIE = if_else(coef_setting == "small", 3*16*n^(-1/4), 2.4),
  NDE = 2,
  # coverage_NDE = (lower_NDE <= NDE) & (NDE <= upper_NDE),
  # coverage_NIE = (lower_NIE <= NIE) & (NIE <= upper_NIE),
  scenario = toupper(scenario)
  # inf_method = if_else(inf_method == "minner", "minnier", inf_method) # fix typo
)

# all_results <- qread("./results/all-results-ml0-9-0.qs", nthreads=2)

# check_df <- all_results %>% 
#   filter(model_version=="full") %>% 
#   group_by(n, p, num_simulations, 
#            scenario, coef_setting, use_sl, suffix, 
#            cores, model_version) %>% 
#   summarize(bias_NDE = abs(mean(err_NDE)), 
#             bias_NIE = abs(mean(err_NIE))) %>% 
#   mutate_at(vars(starts_with("bias_")), ~sqrt(n)*.)
# 
# check_df %>% filter(scenario == "lll") %>% 
#   filter(coef_setting == "large") %>% 
#   tidyr::pivot_longer(starts_with("bias_"),
#                                  names_prefix = "bias_",
#                                  names_to="target", values_to="bias") %>% 
#   ggplot(aes(x=n, y=bias, color=use_sl)) +
#   geom_point(size=2.5, alpha=0.7) +
#   scale_x_log10() +
#   facet_grid(rows=vars(target))
# 
# 
# 
# check_df <- all_results %>% 
#   filter(scenario=="lll", model_version=="product", use_sl) %>% 
#   dplyr::group_by(
#     n, p, num_simulations,
#     scenario, coef_setting, use_sl, suffix,
#     cores, model_version, inf_method) %>%
#   dplyr::summarise(coverage_NDE = mean(coverage_NDE),
#                    coverage_NIE = mean(coverage_NIE),
#                    .groups = "keep")
# 
# check_df %>%
#   tidyr::pivot_longer(starts_with("coverage_"),
#                       names_prefix = "coverage_",
#                       names_to="target", values_to="coverage") %>%
#   ggplot(aes(x=n, y=coverage, color=inf_method)) +
#   geom_hline(yintercept=0.95) +
#   geom_point(size=2.5, alpha=0.7) +
#   scale_x_log10() +
#   facet_grid(rows=vars(target))



bias_df <- all_results %>% filter(inf_method=="naive_boot") %>%
  group_by(n, p, num_simulations, 
           scenario, coef_setting, use_sl, suffix, 
           cores, model_version) %>% 
  summarize(bias_NDE = abs(mean(err_NDE)), 
            bias_NIE = abs(mean(err_NIE))) %>% 
  tidyr::pivot_longer(starts_with("bias_"),
                      names_prefix = "bias_",
                      names_to="target", values_to="bias")


plot_bias <- function(bias_df, coefsize, geom, normalize=T){
  
  assertthat::assert_that(coefsize %in% c("Small", "Large"))
  assertthat::assert_that(geom %in% c("point", "line"))
  assertthat::assert_that(is.logical(normalize))
  lower_coefsize <- tolower(coefsize)
  
  plot_df <- bias_df %>% 
    filter(coef_setting == lower_coefsize)
  
  if(normalize){
    plot_df <- plot_df %>% mutate_at(vars(starts_with("bias_")), ~sqrt(n)*.)
  }
  
  the_geom <- if(geom == "point"){
    geom_point(alpha=0.7, size=2)
  } else {
    geom_line()
  }
  
  the_scale <- if(geom == "point"){
    scale_shape_discrete(labels=c("True","Estimated"), 
                         name=paste0("\u03bc Method"))
  } else {
    scale_linetype_discrete(labels=c("True","Estimated"), 
                            name=paste0("\u03bc Method"))
  }
  
  plot_df %>% 
    mutate(xaxis=n) %>%
    ggplot(aes(x=xaxis, y=bias, color=model_version, 
               shape=use_sl, linetype=use_sl)) +
    the_geom +
    scale_x_log10(breaks=500*2^(0:3)) +
    facet_grid(rows=vars(target), cols=vars(scenario)) +
    theme_bw() +
    the_scale +
    labs(title=paste0(if_else(normalize, "\u221an ", ""),
                      "Absolute Bias vs Sample Size by Estimation Procedure"),
         subtitle=paste0(coefsize, " Setting"),
         x="n", y=paste0(if_else(normalize, "\u221an ", ""),
                         "Absolute Bias"), 
         color="Estimator")
}

plot_bias(bias_df, "Large", "line")
plot_bias(bias_df, "Small", "line")

plot_coverage <- function(cover_df, coefsize, geom) {
  
  assertthat::assert_that(coefsize %in% c("Small", "Large"))
  assertthat::assert_that(geom %in% c("point", "line"))
  lower_coefsize <- tolower(coefsize)
  
  plot_df <- cover_df %>% 
    filter(coef_setting == lower_coefsize)
  
  the_geom <- if(geom == "point"){
    geom_point(alpha=0.7, size=2)
  } else {
    geom_line()
  }
  
  the_scale <- if(geom == "point"){
    scale_shape_discrete(labels=c("True","Estimated"), 
                         name=paste0("\u03bc Method"))
  } else {
    scale_linetype_discrete(labels=c("True","Estimated"), 
                            name=paste0("\u03bc Method"))
  }
  
  plot_df %>% 
    ggplot(aes(x=n, y=coverage, color=model_version,
               shape=use_sl, linetype=use_sl)) +
    geom_hline(yintercept=0.95) +
    the_geom +
    scale_x_log10(breaks=500*2^(0:3)) +
    facet_grid(rows=vars(target), cols=vars(scenario)) +
    theme_bw() +
    the_scale +
    labs(title="Coverage Rates for Minnier Bootstrap", 
         subtitle=paste0(coefsize, " Setting"),
         x="n", y="Coverage",
         color="Estimator")
}


cover_df <- all_results %>% 
  dplyr::group_by(
    n, p, num_simulations,
    scenario, coef_setting, use_sl, suffix,
    cores, model_version, inf_method
  ) %>%
  dplyr::summarise(coverage_NDE = mean(coverage_NDE),
                   coverage_NIE = mean(coverage_NIE),
                   .groups = "keep") %>% 
  tidyr::pivot_longer(starts_with("coverage_"),
                      names_prefix = "coverage_",
                      names_to="target", values_to="coverage")

cover_df$inf_method %>% unique
cover_df %>% filter(inf_method=="minnier") %>% plot_coverage("Large", "line")
cover_df %>% filter(inf_method=="minnier") %>% plot_coverage("Small", "line")



#### Create Variable Selection Performance Table ####
# all_results %>% head %>% tidyr::unnest(sel_info) %>% View
# 
# library(foreach)
# sel_perf <- foreach(i=seq_along(all_results$sel_info), .combine=rbind) %do% {
#   num_missed <- length(
#     setdiff(
#       1:3, all_results$sel_info[[i]]$sel_M
#     )
#   )
#   
#   num_noise <- length(
#     setdiff(
#       all_results$sel_info[[i]]$sel_M, 1:3
#     )
#   )
#   c(num_missed, num_noise)
# }
# colnames(sel_perf) <- c("num_missed", "num_noise")
# all_results <- cbind(all_results, sel_perf)


all_results %>% filter(inf_method=="naive_boot") %>% 
  filter(model_version %in% c("product", "mixture", "adaptive")) %>% 
  group_by(coef_setting, n, model_version, use_sl) %>% 
  summarize(PC = mean(num_missed == 0), 
            MN = median(num_noise)) %>% 
  # mutate(scenario=toupper(scenario),
  #        val=if_else(col == "is_missed", 1-val, val),
  #        col=if_else(col == "is_missed", "PC", "MN")) %>% 
  # tidyr::pivot_wider(id_cols=c(coefs, n, weight_name, SL), 
  #                    names_from=c(scenario, col),
  #                    names_glue="{scenario}_{col}",
  #                    values_from=val) %>% 
  arrange(coef_setting, n, desc(model_version), use_sl) %>% 
  knitr::kable(format="latex", digits=2L)


# # save the result in qs
# library(qs)
# qs::qsave(all_results, "./results/all-results-ml0-9-0.qs")

