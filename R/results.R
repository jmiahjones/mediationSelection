library(ggplot2)
library(qs)
library(dplyr)


suffix <- "ml0-9-2.qs"
files <- dir("./results/", "*.qs", full.names = T)
files <- files[grep(paste0("result-\\d+-.+", suffix), files, perl=T)]

all_results <- lapply(files, qread, nthreads=2) %>% bind_rows
all_results <- all_results %>% mutate(
  # NIE = if_else(coef_setting == "small", 48*n^(-1/2), 2.4),
  # NDE = 2,
  # coverage_NDE = (lower_NDE <= NDE) & (NDE <= upper_NDE),
  # coverage_NIE = (lower_NIE <= NIE) & (NIE <= upper_NIE),
  scenario = toupper(scenario),
  inf_method = if_else(inf_method == "minner", "minnier", inf_method) # fix typo
)
all_results <- all_results %>% 
  mutate(model_version = factor(model_version, levels=c(
    "product", "mixture", "adaptive", "full", "oracle"
  )))

model_colors <- RColorBrewer::brewer.pal(nlevels(all_results$model_version), 
                                         "Set1")
names(model_colors) <- levels(all_results$model_version)

  

bias_df <- all_results %>% filter(inf_method=="naive_boot") %>%
  group_by(n, p, num_simulations, 
           scenario, coef_setting, use_sl, suffix, 
           cores, model_version) %>% 
  summarize(bias_NDE = abs(mean(err_NDE)), 
            bias_NIE = abs(mean(err_NIE))) %>% 
  tidyr::pivot_longer(starts_with("bias_"),
                      names_prefix = "bias_",
                      names_to="target", values_to="bias")


plot_bias <- function(bias_df, coefsize, geom, model_colors,
                      normalize=T){
  
  assertthat::assert_that(coefsize %in% c("Small", "Large"))
  assertthat::assert_that(geom %in% c("point", "line"))
  assertthat::assert_that(is.logical(normalize))
  lower_coefsize <- tolower(coefsize)
  
  plot_df <- bias_df %>% 
    filter(coef_setting == lower_coefsize)
  
  if(normalize){
    plot_df <- plot_df %>% mutate(bias=sqrt(n)*bias)
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
    scale_color_manual(values=model_colors) +
    the_scale +
    labs(title=paste0(if_else(normalize, "\u221an ", ""),
                      "Absolute Bias vs Sample Size by Estimation Procedure"),
         subtitle=paste0(coefsize, " Setting"),
         x="n", y=paste0(if_else(normalize, "\u221an ", ""),
                         "Absolute Bias"), 
         color="Estimator")
}

png("plots/bias_large.png", width=840, height=550)
plot_bias(bias_df, "Large", "line", model_colors)
dev.off()
png("plots/bias_small.png", width=840, height=550)
plot_bias(bias_df, "Small", "line", model_colors)
dev.off()

plot_coverage <- function(cover_df, coefsize, geom, model_colors) {
  
  assertthat::assert_that(coefsize %in% c("Small", "Large"))
  assertthat::assert_that(geom %in% c("point", "line", "linerange"))
  lower_coefsize <- tolower(coefsize)
  
  plot_df <- cover_df %>% 
    filter(coef_setting == lower_coefsize)
  
  the_geom <- if(geom == "point"){
    geom_point(alpha=0.7, size=2)
  } else if(geom == "linerange"){
    geom_line()
  } else {
    geom_line()
  }
  
  bars <- if(geom == "linerange") { 
    geom_errorbar(position=position_dodge(width=0.1),
                  width=0.1) 
  } else {NULL}
  
  the_scale <- if(geom == "point"){
    scale_shape_discrete(labels=c("True","Estimated"), 
                         name=paste0("\u03bc Method"))
  } else {
    scale_linetype_discrete(labels=c("True","Estimated"), 
                            name=paste0("\u03bc Method"))
  }
  
  plot_df %>% 
    mutate(
      ns = floor(coverage * num_simulations),
      nf = num_simulations - ns,
      wilson.se = (1.96/(num_simulations+1.96^2))*
        sqrt(ns*nf/num_simulations + (1.96^2)/4),
      wilson.center = (ns + 0.5*1.96^2)/(num_simulations + 1.96^2),
      ci.coverage.lower = wilson.center - wilson.se,
      ci.coverage.upper = wilson.center + wilson.se
    ) %>% 
    ggplot(aes(x=n, y=coverage, color=model_version,
               shape=use_sl, linetype=use_sl,
               ymin=ci.coverage.lower, ymax=ci.coverage.upper
               )) +
    geom_hline(yintercept=0.95) +
    the_geom +
    bars +
    scale_x_log10(breaks=500*2^(0:3)) +
    facet_grid(rows=vars(target), cols=vars(scenario)) +
    theme_bw() +
    scale_color_manual(values=model_colors) +
    the_scale +
    labs(title="Coverage Rates for Proposed Bootstrap", 
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

png("plots/coverage_large.png", width=840, height=550)
cover_df %>% filter(inf_method=="minnier" 
                    # | (inf_method=="naive_boot" & model_version=="oracle")
                    ) %>% 
  plot_coverage("Large", "linerange", model_colors)
dev.off()
png("plots/coverage_small.png", width=840, height=550)
cover_df %>% filter(inf_method=="minnier") %>% 
  plot_coverage("Small", "linerange", model_colors)
dev.off()



##############################################################


# biasvar_df <- all_results %>% filter(inf_method=="naive_boot") %>%
#   group_by(n, p, num_simulations, 
#            scenario, coef_setting, use_sl, suffix, 
#            cores, model_version) %>% 
#   summarize(bias_NDE = abs(mean(err_NDE)), 
#             bias_NIE = abs(mean(err_NIE)),
#             sd_NDE = (sd(err_NDE)), 
#             sd_NIE = (sd(err_NIE))
#             ) %>% 
#   tidyr::pivot_longer(starts_with("bias_"),
#                       names_prefix = "bias_",
#                       names_to="target", values_to="bias") %>% 
#   tidyr::pivot_longer(starts_with("sd_"),
#                       names_prefix = "sd_",
#                       names_to="target2", values_to="sd") %>% 
#   filter(target == target2) %>% dplyr::select(!one_of("target2"))
# 
# 
# plot_biasvar <- function(bias_df, coefsize, geom, model_colors,
#                          normalize=T){
#   
#   assertthat::assert_that(coefsize %in% c("Small", "Large"))
#   assertthat::assert_that(geom %in% c("point", "line"))
#   assertthat::assert_that(is.logical(normalize))
#   lower_coefsize <- tolower(coefsize)
#   
#   plot_df <- bias_df %>% 
#     filter(coef_setting == lower_coefsize)
#   
#   if(normalize){
#     plot_df <- plot_df %>% mutate_at(vars(starts_with("bias_"),
#                                           starts_with("sd_")), ~sqrt(n)*.)
#   }
#   
#   the_geom <- if(geom == "point"){
#     geom_point(alpha=0.7, size=2)
#   } else {
#     geom_line()
#   }
#   
#   the_scale <- if(geom == "point"){
#     scale_shape_discrete(labels=c("True","Estimated"), 
#                          name=paste0("\u03bc Method"))
#   } else {
#     scale_linetype_discrete(labels=c("True","Estimated"), 
#                             name=paste0("\u03bc Method"))
#   }
#   
#   plot_df %>% 
#     mutate(xaxis=n) %>%
#     ggplot(aes(x=sd, y=bias, color=model_version, 
#                shape=use_sl, linetype=use_sl)) +
#     the_geom +
#     # scale_x_log10(breaks=500*2^(0:3)) +
#     facet_grid(rows=vars(target), cols=vars(scenario)) +
#     theme_bw() +
#     scale_color_manual(values=model_colors) +
#     the_scale +
#     labs(title=paste0(if_else(normalize, "\u221an ", ""),
#                       "Absolute Bias vs Sample Size by Estimation Procedure"),
#          subtitle=paste0(coefsize, " Setting"),
#          x="SD", y=paste0(if_else(normalize, "\u221an ", ""),
#                          "Absolute Bias"), 
#          color="Estimator")
# }
# 
# plot_biasvar(biasvar_df, "Large", "point", model_colors, normalize = F)
