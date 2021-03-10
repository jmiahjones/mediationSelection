library(qs)
library(foreach)

################# Call Main ########################
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  warning("No commandline arguments found. Debugging...")
  args <- c(
    100, # n
    20, # num_simulations
    "lll", # abbv_scn
    "large", # coef_setting
    "cv", # weight_gam
    "T", # use_sl
    "debug", # suffix_arg
    2 # cores
  )
}

DEFAULT_CORES <- 20

stopifnot(length(args) %in% (7:8))
n <- as.numeric(args[1])
num_simulations <- as.numeric(args[2])
abbv_scn <- args[3]
coef_setting <- args[4]
weight_gam <- as.numeric(args[5])
use_sl <- args[6]
suffix_arg <- args[7]
cores <- as.numeric(args[8])

assertthat::assert_that(use_sl %in% c("F","T"))
use_sl_str <- use_sl
use_sl <- use_sl == "T"

if(is.na(weight_gam)){
  # weight_gam was not numeric, so assume we want to cv it
  weight_gam_str <- gsub(" ", "", args[5])
  weight_gam <- c(1,2,3)
} else {
  weight_gam_str <- weight_gam
}

if(is.na(cores)){
  cores <- DEFAULT_CORES
}

savefile_suffix <- paste(n, num_simulations, abbv_scn, coef_setting,
                         weight_gam_str, use_sl_str, suffix_arg, sep="-")

save_file <- paste0("./results/result-",
                    savefile_suffix,
                    ".qs"
                    )

source("./R/var-selection-sim-cv.R")
print(paste0("Cores: ", cores))

sl_loop <- if(use_sl) c(T, F) else F
main_result <- foreach(slB=sl_loop, .combine=rbind) %do%{
  main(
    n,
    num_simulations,
    abbv_scn,
    coef_setting,
    weight_gam,
    slB,
    suffix_arg,
    cores
  )
}

# save the result in RDS
qs::qsave(main_result, save_file)

