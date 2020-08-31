# min.trees = 1000#setting this in Rmd file
# max.trees = 99000
if (testing == 1){
  eta = 0.001
  max_depth = 1
  nruns = 2
  nrounds = 3000
  } else {
  eta = c(0.0001, 0.001, 0.01)
  max_depth = c(1, 2, 3, 4)
  nruns = nruns
  nrounds = 100000
}
TRAIT <- read.csv("input/mammal_imputation_data.csv", header = T)

sample_size = dim(TRAIT)[1]#if there's a small sample size then set n.minobsinnode to exclude 10
if (sample_size>1000){
  n.minobsinnode = c(2,5,10)
} else {
  n.minobsinnode = c(2,5)
}
#Would recommend setting the working directory to a new folder for imputation so that all of the saved objects here do not mix with other files
packages <- c("tidyverse", "magrittr", "Hmisc", "gghighlight", "caret", "gbm", "Matrix", "pdp", "caTools", "ROCR", "foreach", "dismo", "doSNOW", "parallel", "patchwork", "reshape2", "gtable")
sapply(packages, library, character.only = T)

# detectCores()
# cl <- makeCluster(6, "SOCK")
# registerDoSNOW(cl)
setwd(wd_set)#set this outside this file

col_names_inds = which(names(TRAIT) %in% lab_list)
#only included variables that may be possible to impute, so I got rid of infant mortality rate, etc. with over 100 NAs.
balt <- foreach(LAB = lab_list, .combine = "rbind", .errorhandling = "pass") %do% { #kept an errorhandling of pass to deal with anything that probably will come up
  alt <- TRAIT[which(!is.na(TRAIT[, LAB])), ] #subset to get rid of NAs that will gum up modeling
  # VAR <- colnames(alt)[c(8:40, 43:68)] #all variables considered
  VAR <- colnames(alt)[col_names_inds] #all variables considered
  VAR <- VAR[-which(VAR == LAB)] #remove label
  GRID <- gridSearch(alt, label = LAB, eta = eta, max_depth = max_depth, n.minobsinnode = n.minobsinnode, vars = VAR, k_split = 0.8, distribution = "gaussian", nrounds = nrounds)
  PRMS <- GRID[[1]] %>% filter(n.trees > min.trees & n.trees < max.trees) #filter parameters considered
  PRMS <- PRMS[which.max(PRMS$eval_test), 1:3] #choose according to maximized test eval
  save(GRID, file = paste(LAB, "_grid_results.Rdata", sep = ""))
  BOOT <- bootstrapGBM(DF = alt, label = LAB, vars = VAR, k_split = 0.8, distribution = "gaussian", eta = PRMS[1], max_depth = PRMS[2], nrounds = nrounds, nruns = nruns, bootstrap = "observed", method = "cv", cv.folds = 5, n.minobsinnode = PRMS[3], file_label=LAB, id_field = id_field)
  save(BOOT, file = paste(LAB, "_bootstrap_results.Rdata", sep = ""))
  NULB <- bootstrapGBM(DF = alt, label = LAB, vars = VAR, k_split = 0.8, distribution = "gaussian", eta = PRMS[1], max_depth = PRMS[2], nrounds = nrounds, nruns = nruns, bootstrap = "null", method = "cv", cv.folds = 5, n.minobsinnode = PRMS[3], file_label=paste(LAB, "null", sep = "_"), id_field = id_field)
  save(NULB, file = paste(LAB, "_nullbootstrap_results.Rdata", sep = ""))
  NAMES <- names(which(apply(BOOT[[1]][, -1:-6], 2, mean) < 1)) #which values are under a 1 importance
  VAR <- VAR[!VAR %in% NAMES] #remove values from consideration to rerun bootstrap
  BOOTR <- bootstrapGBM(DF = alt, label = LAB, vars = VAR, k_split = 0.8, distribution = "gaussian", eta = PRMS[1], max_depth = PRMS[2], nrounds = nrounds, nruns = nruns, bootstrap = "observed", method = "cv", cv.folds = 5, n.minobsinnode = PRMS[3] , file_label= paste(LAB, "min", sep = "_"), id_field =id_field)
  save(BOOTR, file = paste(LAB, "_bootstrap_results_min1.Rdata", sep = ""))
  #write.csv(BOOTR[[2]], "femmd_model_bootstrapmin_results2.csv", row.names = F)
  NULBR <- bootstrapGBM(DF = alt, label = LAB, vars = VAR, k_split = 0.8, distribution = "gaussian", eta = PRMS[1], max_depth = PRMS[2], nrounds = nrounds, nruns = nruns, bootstrap = "null", method = "cv", cv.folds = 5, n.minobsinnode = PRMS[3], file_label=paste(LAB, "null", "min", sep = "_"), id_field = id_field)
  save(NULBR, file = paste(LAB, "_nullbootstrap_results_min1.Rdata", sep = ""))
  cbind.data.frame(var = LAB, nsamp = nrow(alt), eta = PRMS[1], max_depth = PRMS[2], n.minobsinnode = PRMS[3], bootEval = mean(BOOT[[1]]$eval_test), nullBootEval = mean(NULB[[1]]$eval_test), min1BootEval = mean(BOOTR[[1]]$eval_test), min1NullBootEval = mean(NULBR[[1]]$eval_test))
}