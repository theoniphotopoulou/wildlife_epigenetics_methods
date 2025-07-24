library(doParallel)
library(dplyr)
library(tidyr)

# make functions folder
source("fgecs.R")
source("probage.R")
source("healthmodels.R")
source("probage_fns.R")
source("age-transformations-zoller.R")

set.seed(47)

load("output/dolphin_data.Rdata")
health$VESOP_Flag = ifelse(health$VESOP_surv < 0.925, 1, 0)

#################
## ELASTIC NET ##
#################

# All combinations of parameters that we vary, each row is an elastic net model that is run
pars = expand.grid(use_weights = c(FALSE, TRUE), controls_only = c(FALSE, TRUE), reduce_size = c(FALSE, TRUE), 
                   loglin_age = c("lin", "loglin"), loglin_sexm = 15, 
                   en_alpha = c(0.01, seq(from = 0.1, to = 0.9, by = 0.1), 0.99),
                   n_cpg = c(500, 1000, 2000, 31157))
pars = pars |> filter(!(controls_only & reduce_size))
pars = pars |> mutate(par_set = 1:n())

# Run on bluewhale
cl <- makeCluster(20)
registerDoParallel(cl)
getDoParWorkers()

n_en = nrow(pars)
res <- foreach(kk = 1:n_en, 
               .packages = c("glmnet", "dplyr", "tidyr","stringr", 
                             "randomForest", "survival", "flexsurv",
                             "data.table")) %dopar% {
                               
                               source("age-transformations-zoller.R")
                               print(kk)
                               
                               m0 <- run_en(x = x, sample_info = sample_info, n_cpg = pars$n_cpg[kk], en_alpha = pars$en_alpha[kk], use_weights = pars$use_weights[kk], 
                                            controls_only = pars$controls_only[kk], reduce_size = pars$reduce_size[kk], loglin_age = pars$loglin_age[kk], 
                                            loglin_sexm = pars$loglin_sexm[kk])
                               
                               alltest <- m0$testpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               alltrain <- m0$trainpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               allval <- m0$valpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               
                               list(alltest = alltest, alltrain = alltrain, allval = allval)
                             }

stopCluster(cl)

res_en = res
pars_en = pars
save(res_en, pars_en, file = "output/en_model_output.Rdata")
rm(res, pars)

#################
## RANDOM FOREST ##
#################

# All combinations of parameters that we vary, each row is a RF model that is run
pars = expand.grid(use_weights = c(FALSE, TRUE), controls_only = c(FALSE, TRUE), reduce_size = c(FALSE, TRUE), 
                   loglin_age = c("lin", "loglin"), loglin_sexm = 15, 
                   rf_mtry = c(1/6, 1/3, 1/2),
                   n_cpg = c(500, 1000, 2000, 31157))
pars = pars |> filter(!(controls_only & reduce_size))
pars = pars |> mutate(par_set = 1:n())

# Run on bluewhale
cl <- makeCluster(20)
registerDoParallel(cl)
getDoParWorkers()

nn = nrow(pars)
res <- foreach(kk = 1:nn, 
               .packages = c("glmnet", "dplyr", "tidyr","stringr", 
                             "randomForest", "survival", "flexsurv",
                             "data.table")) %dopar% {
                               
                               source("age-transformations-zoller.R")
                               print(kk)
                               
                               m0 <- run_rf(x = x, sample_info = sample_info, rf_mtry = pars$rf_mtry[[kk]], use_weights = pars$use_weights[kk], 
                                            controls_only = pars$controls_only[kk], reduce_size = pars$reduce_size[kk], loglin_age = pars$loglin_age[kk], 
                                            loglin_sexm = pars$loglin_sexm[kk], n_cpg = pars$n_cpg[kk])
                               
                               alltest <- m0$testpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               alltrain <- m0$trainpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               allval <- m0$valpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               
                               list(alltest = alltest, alltrain = alltrain, allval = allval)
                             }

stopCluster(cl)

res_rf = res
pars_rf = pars
save(res_rf, pars_rf, file = "output/rf_model_output.Rdata")
rm(res, pars)

#################
## PROBAGE ##
#################

pars = expand.grid(use_weights = c(FALSE, TRUE), controls_only = c(FALSE, TRUE), reduce_size = c(FALSE, TRUE), 
                   loglin_age = c("lin", "loglin"), loglin_sexm = 15)
pars = pars |> filter(!(controls_only & reduce_size))
pars = pars |> mutate(par_set = 1:n())

# i = 1 is the easy one with everything set to FALSE

cl <- makeCluster(20)
registerDoParallel(cl)
getDoParWorkers()

n_en = nrow(pars)
res <- foreach(kk = 1:n_en, 
               .packages = c("glmnet", "dplyr", "tidyr","stringr", 
                             "randomForest", "survival", "flexsurv",
                             "data.table")) %dopar% {
                               
                               source("age-transformations-zoller.R")
                               print(kk)
                               
                               m0 <- run_probage(x = x, sample_info = sample_info, use_weights = pars$use_weights[kk], controls_only = pars$controls_only[kk],
                                                 reduce_size = pars$reduce_size[kk], loglin_age = pars$loglin_age[kk], loglin_sexm = pars$loglin_sexm[kk])
                               
                               alltest <- m0$testpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               alltrain <- m0$trainpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               
                               list(alltest = alltest, alltrain = alltrain)
                             }

stopCluster(cl)

res_pr = res
pars_pr = pars
save(res_pr, pars_pr, file = "output/pr_model_output.Rdata")
rm(res, pars)


#########################
## HEALTH ELASTIC NET ##
#########################

pars = expand.grid(use_weights = c(FALSE), controls_only = c(FALSE), reduce_size = c(FALSE), 
                   yvar = c("VESOP_Flag", "Lung_Flag", "Anemia_Flag", "Cortisol_lt30_Flag", "Neutrophils_Flag", 
                            "Globulin_Flag", "Glucose_Flag"),
                   y_transform = "lin",
                   en_family = "binomial",
                   en_alpha = c(0.01, seq(from = 0.1, to = 0.9, by = 0.1), 0.99),
                   n_cpg = c(500, 1000, 2000, 31157),
                   stringsAsFactors = FALSE) 
pars = pars |> filter(!(controls_only & reduce_size))
pars = pars |> mutate(par_set = 1:n())

# i = 1 is the easy one with everything set to FALSE

cl <- makeCluster(30)
registerDoParallel(cl)
getDoParWorkers()

nn = nrow(pars)
res <- foreach(kk = 1:nn, 
               .packages = c("glmnet", "dplyr", "tidyr","stringr", 
                             "randomForest", "survival", "flexsurv",
                             "data.table")) %dopar% {
                               
                               print(kk)
                               
                               m0 <- run_en_health(x = x, yvar = pars$yvar[kk], 
                                                   y_transform = pars$y_transform[kk], en_family = pars$en_family[kk], 
                                                   health = health, sample_info = sample_info, 
                                                   en_alpha = pars$en_alpha[kk], n_cpg = pars$n_cpg[kk], 
                                                   use_weights = pars$use_weights[kk], controls_only = pars$controls_only[kk], reduce_size = pars$reduce_size[kk])
                               
                               alltest <- m0$testpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               alltrain <- m0$trainpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               allval <- m0$valpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               
                               list(alltest = alltest, alltrain = alltrain, allval = allval)
                             }

stopCluster(cl)

res_en = res
pars_en = pars
save(res_en, pars_en, file = "output/health_en_model_output.Rdata")
rm(res, pars)

#########################
## HEALTH RANDOM FOREST ##
#########################

pars = expand.grid(use_weights = c(FALSE), controls_only = c(FALSE), reduce_size = c(FALSE), 
                   yvar = c("VESOP_Flag", "Lung_Flag", "Anemia_Flag", "Cortisol_lt30_Flag", "Neutrophils_Flag", 
                            "Globulin_Flag", "Glucose_Flag"),
                   y_transform = "factor",
                   rf_mtry = c(1/6, 1/3, 1/2),
                   n_cpg = c(500, 1000, 2000, 31157),
                   stringsAsFactors = FALSE)
pars = pars |> filter(!(controls_only & reduce_size))
pars = pars |> mutate(par_set = 1:n())

# i = 1 is the easy one with everything set to FALSE

cl <- makeCluster(30)
registerDoParallel(cl)
getDoParWorkers()

nn = nrow(pars)
res <- foreach(kk = 1:nn, 
               .packages = c("glmnet", "dplyr", "tidyr","stringr", 
                             "randomForest", "survival", "flexsurv",
                             "data.table")) %dopar% {
                               
                               print(kk)
                               
                               m0 <- run_rf_health(x = x, yvar = pars$yvar[kk], 
                                                   y_transform = pars$y_transform[kk], 
                                                   health = health, sample_info = sample_info, 
                                                   rf_mtry = pars$rf_mtry[kk], n_cpg = pars$n_cpg[kk], 
                                                   use_weights = pars$use_weights[kk], controls_only = pars$controls_only[kk], reduce_size = pars$reduce_size[kk])
                               
                               alltest <- m0$testpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               alltrain <- m0$trainpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               allval <- m0$valpreds |> left_join(sample_info |> select(DOLPHIN.ID, age), by = join_by(DOLPHIN.ID)) |> mutate(par_set = kk)
                               
                               list(alltest = alltest, alltrain = alltrain, allval = allval)
                             }

stopCluster(cl)

res_rf = res
pars_rf = pars
save(res_rf, pars_rf, file = "output/health_rf_model_output.Rdata")
rm(res, pars)

