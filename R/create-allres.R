library(foreach)

tmp.dir <- "./cache/temp/"
filenames <- dir(tmp.dir, pattern = "*.RData")
filenames <- paste0(tmp.dir, filenames)
all_res_save <- "./cache/all_res.RData"


all_res <- foreach(file=filenames, .combine=rbind, .errorhandling="stop") %do% { load(file); return(temp_res) }

save(all_res, file=all_res_save)
print("Complete! Verify no errors before erasing cache...")
