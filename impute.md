impute
================
Han lab
8/31/2020

### install and load packages needed for study

``` r
print(Sys.time())
neededPackages <- c( "gbm", "caret", "Matrix", "pdp", "caTools", "ROCR", "dplyr", "foreach", "dismo", "doSNOW", "parallel", 
                     "tidyverse", "magrittr", "Hmisc", "gghighlight", "dismo",  "patchwork", "reshape2", "gtable")

                  
# pkgTest is a helper function to load packages and install packages only when they are not installed yet.

pkgTest <- function(x)
{
  if(x %in% rownames(installed.packages()) == FALSE) {
    install.packages(x, dependencies= TRUE, repos = "http://cran.us.r-project.org")    
  }
  library(x, character.only = TRUE)
}

for(package in neededPackages){pkgTest(package)}
```

``` r
testing = 1#if testing == 1 then only do a small subset of the variables
lab_list = c("female_maturity_d", "male_maturity_d", "weaning_d", "development_d", 
                        "log_litterclutch_size_n", "litters_or_clutches_per_y", "log_inter_litterbirth_interval_y",
                        "log_birthhatching_weight_g", "log_weaning_weight_g", "metabolicRate_W", "temperature_K", 
                        "longevity_y", "log_male_body_mass_g", "adult_svl_cm", "mass_specific_production", "log_range_size", "tnc_ecoregion_breadth")
max.trees = 100000
if (testing == 1){
  lab_list = c("development_d","male_maturity_d" , "longevity_y" )#pick this one for which imputation has high accuracy
  min.trees = 0
} else {
  min.trees = 1000
}

output_csv = "output_csv/"

folder = "output/"
folder_graphs = "graphs/"
#remove graph files
files_remove = list.files(folder_graphs)
file.remove(paste0(folder_graphs, files_remove))
```

    ## [1] TRUE TRUE TRUE

``` r
accuracy_threshold = 0.8#corrected test eval must be above this for us to include the imputed values where missing
NAMES = lab_list
```

\#\#cores

\#source
gridSearch.R

``` r
source("https://raw.githubusercontent.com/fischhoff/fishbase2/master/gridSearch.R")
```

\#source
bootstrapGBM.R

``` r
source("https://raw.githubusercontent.com/fischhoff/fishbase2/master/bootstrapGBM.R")
```

\#\#try out bootstrapGBM using Pima dataset
\#\#<https://www.kaggle.com/kumargh/pimaindiansdiabetescsv?select=pima-indians-diabetes.csv>
\#\#label = “X1”

``` r
source("code/pima.R")
```

    ## 'data.frame':    1534 obs. of  5 variables:
    ##  $ predictions   : num  0.0515 0.6617 0.0195 0.7175 0.1331 ...
    ##  $ bootstrap_run : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ id_field_vals : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ label         : chr  "X1" "X1" "X1" "X1" ...
    ##  $ original_value: int  0 1 0 1 0 1 0 1 1 0 ...

``` r
temp = OUT[[3]]
names(temp)
```

    ## [1] "predictions"    "bootstrap_run"  "id_field_vals"  "label"         
    ## [5] "original_value"

\#source imputation code

``` r
# wd_set = "\\09ies\Han\COVID19\imputation\impute\output"
id_field = "Species"
wd_set = "output"
files_remove = list.files(paste0(wd_set, "/"))
file.remove(paste0(wd_set, "/", files_remove))
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
source("code/MammalImputationLoop.R")
# source("MammalImputationLoop_v2.R")
print(Sys.time())
```

    ## [1] "2020-08-31 12:23:35 EDT"

\#\#compare predicted to actual

``` r
PREDI <- foreach(x = NAMES) %do% {
load(paste0(folder, x, "_bootstrap_results.Rdata"))
load(paste0(folder, x, "_bootstrap_results_min1.Rdata"))
boot <- BOOT[[3]] %>% group_by(id_field_vals) %>% summarise(mean_pred = mean(predictions),
                                                            sd_pred = sd(predictions),
                                                            orig_vals = original_value[1]) %>%
  mutate(type = "all_vars")

BOOTR[[3]] %>% group_by(id_field_vals) %>% summarise(mean_pred = mean(predictions),
                                                            sd_pred = sd(predictions),
                                                            orig_vals = original_value[1]) %>%
  mutate(type = "min1_vars") %>% rbind(boot) %>% ggplot() +
  
   geom_crossbar(aes(y = mean_pred, ymin = mean_pred - sd_pred, ymax = mean_pred + sd_pred, x = id_field_vals, color = type, group = type), position = position_dodge(width = 1)) +
  scale_color_manual(values = c("darkcyan", "mediumpurple")) +
  geom_point(aes(x = id_field_vals, y = orig_vals), shape = 21, fill = "darkorange", color = "darkorange3", size = 2.5, alpha = 0.1) +
  labs(x = "Record", y = x) +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = "transparent", size = 1.2), panel.grid.major = element_line(color = "grey80"))+
  theme(axis.text.x = element_text(angle = 90))
}
PREDI
```

    ## list()

\#\#look at performance

``` r
a  = 1
out_perf = NULL
for (a in 1:length(lab_list)){
  obs = paste0(folder, lab_list[a], "_bootstrap_results.Rdata")
  load(obs)
  P <- BOOT[[1]]
  source("code/R_get_performance.R")
  tmp_sum$predictor_variables_included = "all"

  ##now get null data
  obs = paste0(folder, lab_list[a], "_nullbootstrap_results.Rdata")
  load(obs)  
  P <- NULB[[1]]#get performance for null data 
  eval_test_null = mean(P$eval_test)
  tmp_sum$eval_test_null = eval_test_null
  #add data
  out_perf = rbind(out_perf, tmp_sum)
  
  #now do same for model w/ only variables w/ importance over 1 percent
  obs = paste0(folder, lab_list[a], "_bootstrap_results_min1.Rdata")
  load(obs)
  P <- BOOTR[[1]]
  source("code/R_get_performance.R")
  tmp_sum$predictor_variables_included = "vars w import > 1"
  
    ##now get null data -- min1
  # obs = paste0(folder, lab_list[a], "_nullbootstrap_results_min1.Rdata")
  # load(obs)
  # P <- NULBR[[1]]#get performance for null data -- min1
  # eval_test_null = mean(P$eval_test)
  tmp_sum$eval_test_null = NA#eval_test_null
  out_perf = rbind(out_perf, tmp_sum)
  
}

out_perf$eval_test_corrected = out_perf$eval_test - abs(out_perf$eval_test_null)
write.csv(out_perf, file= paste0(output_csv,"performance_imputation.csv"), row.names = FALSE)

perf_high = subset(out_perf, eval_test_corrected > accuracy_threshold)
unique(perf_high$variable)
```

    ## character(0)

\#\#look at only deviance curve for parameters selected based on grid
search – add to the plot the name of the predictor variable that’s being
imputed.

``` r
NAMES <- lab_list
####Deviance Curves
GRIDS <- foreach(x = NAMES) %do% {
load(paste0(folder, x, "_grid_results.Rdata"))
# PRMS <- GRID[[1]] %>% filter(n.trees > min.trees & n.trees < max.trees) #filter parameters considered
PRMS <- subset(GRID[[1]], n.trees > min.trees & n.trees < max.trees) #filter parameters considered
# PRMS <- PRMS[which.max(PRMS$eval_test), 1:3] #choose according to maximized test eval  
PRMS <- PRMS[which.max(PRMS$eval_test), ] #choose according to maximized test eval  
GRID[[2]] = subset(GRID[[2]], group == PRMS$group)#get just this part of the deviance curve data
patchwork::wrap_plots(lapply(1:length(unique(GRID[[2]]$group)), function(i) GRID[[2]] %>% filter(group == unique(GRID[[2]]$group)[i]) %>% ggplot() +
  geom_line(aes(x = index, y = train), color = "black", size = 1) +
  geom_line(aes(x = index, y = valid), color = "green", size = 1) +
    geom_vline(xintercept = GRID[[2]] %>% filter(group == unique(GRID[[2]]$group)[i]) %>% dplyr::select(best.iter) %>% unique %>% as.numeric, color = "blue", linetype = "dashed", size = 1) +
  labs(x = "Iteration", y = "Bernoulli deviance", title = paste(x,unique(GRID[[2]]$group)[i])) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1), panel.grid.major = element_line(color = "grey90"), title = element_text(size = 8))))
}
GRIDS
```

    ## [[1]]

![](impute_files/figure-gfm/dev%20best-1.png)<!-- -->

    ## 
    ## [[2]]

![](impute_files/figure-gfm/dev%20best-2.png)<!-- -->

    ## 
    ## [[3]]

![](impute_files/figure-gfm/dev%20best-3.png)<!-- -->

\#\#look at all deviance curves

``` r
####Deviance Curves
GRIDS <- foreach(x = NAMES) %do% {
load(paste0(folder, x, "_grid_results.Rdata"))
patchwork::wrap_plots(lapply(1:length(unique(GRID[[2]]$group)), function(i) GRID[[2]] %>% filter(group == unique(GRID[[2]]$group)[i]) %>% ggplot() +
  geom_line(aes(x = index, y = train), color = "black", size = 1) +
  geom_line(aes(x = index, y = valid), color = "green", size = 1) +
    geom_vline(xintercept = GRID[[2]] %>% filter(group == unique(GRID[[2]]$group)[i]) %>% dplyr::select(best.iter) %>% unique %>% as.numeric, color = "blue", linetype = "dashed", size = 1) +
  labs(x = "Iteration", y = "Bernoulli deviance", title = unique(GRID[[2]]$group)[i]) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = "transparent", color = "black", size = 1), panel.grid.major = element_line(color = "grey90"), title = element_text(size = 8))))
}
GRIDS
```

    ## [[1]]

![](impute_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

    ## 
    ## [[2]]

![](impute_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

    ## 
    ## [[3]]

![](impute_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

\#source
partial\_plot

``` r
source("https://raw.githubusercontent.com/fischhoff/fishbase2/master/partial_plotR.R")
```

\#\#get PDP results for variables w/ importance over 1, for models
including all variables

``` r
PDP <- foreach(x = NAMES, .errorhandling = "pass") %do% {
load(paste0(folder, x, "_bootstrap_results.Rdata"))
load(paste0(folder, "observedhist_", x, ".Rdata"))
VARS <- names(sort(apply(BOOT[[1]][, -1:-6], 2, mean)[apply(BOOT[[1]][, -1:-6], 2, mean) >= 1], decreasing = T))
tmp <- partial_plot(BOOT[[2]], hist.data = out_hist, vars = VARS, type = "mean", histogram = T)
ggsave(plot = tmp, filename = paste0(folder_graphs, "PDP", x,".jpg"), width = 10, height = 10)
#grab those variables over 1 importance and sort names by importance
partial_plot(BOOT[[2]], hist.data = out_hist, vars = VARS, type = "mean", histogram = T)
}
PDP
```

    ## [[1]]

![](impute_files/figure-gfm/pdp%20all-1.png)<!-- -->

    ## 
    ## [[2]]

![](impute_files/figure-gfm/pdp%20all-2.png)<!-- -->

    ## 
    ## [[3]]

![](impute_files/figure-gfm/pdp%20all-3.png)<!-- -->

\#\#get PDP results for variables w/ importance over 1, for models
including variables w/ importance over 1 in the first model

``` r
PDPM <- foreach(x = NAMES) %do% {
load(paste0(folder, x, "_bootstrap_results_min1.Rdata"))
load(paste0(folder, "observedhist_", x, "_min", ".Rdata"))
VARS <- names(sort(apply(BOOTR[[1]][, -1:-6], 2, mean)[apply(BOOTR[[1]][, -1:-6], 2, mean) >= 1], decreasing = T)) #grab those variables over 1 importance and sort names by importance
#save plot
tmp <- partial_plot(BOOTR[[2]], hist.data = out_hist, vars = VARS, type = "mean", histogram = T)
ggsave(plot = tmp, filename = paste0(folder_graphs, "PDP_min1.", x,".jpg"), width = 10, height = 10)

partial_plot(BOOTR[[2]], hist.data = out_hist, vars = VARS, type = "mean", histogram = T)
}
PDPM
```

    ## [[1]]

![](impute_files/figure-gfm/pdp%20min1-1.png)<!-- -->

    ## 
    ## [[2]]

![](impute_files/figure-gfm/pdp%20min1-2.png)<!-- -->

    ## 
    ## [[3]]

![](impute_files/figure-gfm/pdp%20min1-3.png)<!-- -->
