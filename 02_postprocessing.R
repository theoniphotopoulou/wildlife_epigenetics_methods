library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

pars = expand.grid(use_weights = c(FALSE, TRUE), controls_only = c(FALSE, TRUE), reduce_size = c(FALSE, TRUE), 
                   loglin_age = c("lin", "loglin"), loglin_sexm = 15)
pars = pars |> filter(!(controls_only & reduce_size)) 
pars$par_set = 1:12

###############
# Elastic net #
###############

load("output/en_model_output.Rdata")
alltest = do.call(rbind, purrr::map(res_en, "alltest"))
alltrain = do.call(rbind, purrr::map(res_en, "alltrain"))
allval = do.call(rbind, purrr::map(res_en, "allval"))

# identify best combination of (alpha, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means lowest MAE on the validation set
valres = allval |> 
  pivot_longer(cols = starts_with("preds"), names_to = "fold", values_to = "predage") |> 
  mutate(ae = abs(predage - age)) |> 
  group_by(par_set) |> 
  summarize(mae = median(ae, na.rm = TRUE)) |> 
  ungroup() |> 
  left_join(pars_en, by = join_by(par_set)) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age)
valres2 = valres |> slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age))
valres2 = valres2 |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig)
best_parsets_en = unlist(valres2$par_set)

# identify best combination of (alpha = 0.5, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means lowest MAE on the validation set
best_parsets_en_n_cpg = allval |> 
  pivot_longer(cols = starts_with("preds"), names_to = "fold", values_to = "predage") |> 
  mutate(ae = abs(predage - age)) |> 
  group_by(par_set) |> 
  summarize(mae = median(ae, na.rm = TRUE)) |> 
  ungroup() |> 
  left_join(pars_en, by = join_by(par_set)) |> 
  filter(en_alpha == 0.5) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age) |> 
  slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age)) |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig) |> 
  pull(par_set)

# identify best combination of (alpha, n_cpg = 31192) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means lowest MAE on the validation set
best_parsets_en_alpha = allval |> 
  pivot_longer(cols = starts_with("preds"), names_to = "fold", values_to = "predage") |> 
  mutate(ae = abs(predage - age)) |> 
  group_by(par_set) |> 
  summarize(mae = median(ae, na.rm = TRUE)) |> 
  ungroup() |> 
  left_join(pars_en, by = join_by(par_set)) |> 
  filter(n_cpg == 31192) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age) |> 
  slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age)) |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig) |> 
  pull(par_set)

# Extracting different model runs at param combinations from global results object res_en
res_en_opt = res_en[best_parsets_en]                                            # optimise both number of preselected CpG and alpha 
res_en_opt_alpha = res_en[best_parsets_en_alpha]                                # optimise alpha with CpG sites set to max (all)
res_en_opt_n_cpg = res_en[best_parsets_en_n_cpg]                                # optimise both number of preselected CpG with alpha set to 0.5 
id_05 = pars_en |> filter(en_alpha == 0.5, n_cpg == 31192) |> pull(par_set)     # no optimisation, results for all CpG sites and alpha set to 0.5
res_en_05 = res_en[id_05]                                                       # no optimisation, results for all CpG sites and alpha set to 0.5

res_en_opt <- map(res_en_opt, ~ map(.x, ~ select(.x, -par_set)))
res_en_opt_alpha <- map(res_en_opt_alpha, ~ map(.x, ~ select(.x, -par_set)))
res_en_opt_n_cpg <- map(res_en_opt_n_cpg, ~ map(.x, ~ select(.x, -par_set)))
res_en_05 <- map(res_en_05, ~ map(.x, ~ select(.x, -par_set)))

###############
# Random forest #
###############

load("output/rf_model_output.Rdata")
alltest = do.call(rbind, purrr::map(res_rf, "alltest"))
alltrain = do.call(rbind, purrr::map(res_rf, "alltrain"))
allval = do.call(rbind, purrr::map(res_rf, "allval"))

# identify best combination of (rf_mtry, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means lowest MAE on the validation set
valres = allval |> 
  pivot_longer(cols = starts_with("preds"), names_to = "fold", values_to = "predage") |> 
  mutate(ae = abs(predage - age)) |> 
  group_by(par_set) |> 
  summarize(mae = median(ae, na.rm = TRUE)) |> 
  ungroup() |> 
  left_join(pars_rf, by = join_by(par_set)) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age)
valres2 = valres |> slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age)) |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig)
best_parsets_rf = unlist(valres2$par_set)

# identify best combination of (rf_mtry = 1/3, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means lowest MAE on the validation set
best_parsets_rf_n_cpg = allval |> 
  pivot_longer(cols = starts_with("preds"), names_to = "fold", values_to = "predage") |> 
  mutate(ae = abs(predage - age)) |> 
  group_by(par_set) |> 
  summarize(mae = median(ae, na.rm = TRUE)) |> 
  ungroup() |> 
  left_join(pars_rf, by = join_by(par_set)) |> 
  filter(rf_mtry == 1/3) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age) |> 
  slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age)) |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig) |> 
  pull(par_set)

# identify best combination of (rf_mtry, n_cpg = 31192) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means lowest MAE on the validation set
best_parsets_rf_mtry = allval |> 
  pivot_longer(cols = starts_with("preds"), names_to = "fold", values_to = "predage") |> 
  mutate(ae = abs(predage - age)) |> 
  group_by(par_set) |> 
  summarize(mae = median(ae, na.rm = TRUE)) |> 
  ungroup() |> 
  left_join(pars_rf, by = join_by(par_set)) |> 
  filter(n_cpg == 31192) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age) |> 
  slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age)) |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig) |> 
  pull(par_set)

# Extracting different model runs at param combinations from global results object res_rf
res_rf_opt = res_rf[best_parsets_rf]                                          # optimise both number of preselected CpG sites and splitting variables  
res_rf_opt_mtry = res_rf[best_parsets_rf_mtry]                                # optimise number of splitting variables 
res_rf_opt_n_cpg = res_rf[best_parsets_rf_n_cpg]                              # optimise number of preselected CpG sites 
id_05 = pars_rf |> filter(rf_mtry == 1/3, n_cpg == 31192) |> pull(par_set)    # no optimisation, results for all CpG sites and splitting variables = 1/3 of all available
res_rf_05 = res_rf[id_05]                                                     # no optimisation, results for all CpG sites and splitting variables = 1/3 of all available

res_rf_opt <- map(res_rf_opt, ~ map(.x, ~ select(.x, -par_set)))
res_rf_05 <- map(res_rf_05, ~ map(.x, ~ select(.x, -par_set)))
res_rf_opt_mtry <- map(res_rf_opt_mtry, ~ map(.x, ~ select(.x, -par_set)))
res_rf_opt_n_cpg <- map(res_rf_opt_n_cpg, ~ map(.x, ~ select(.x, -par_set)))


####### CHECK
# par_set is id for combination of parameters that goes into each model run
# this *might* make them comparable between all the model types so might not be needed now

ll = c(1:12)
for(i in 1:length(res_en_opt)){
  res_en_opt[[i]]$alltest$par_set = ll[i]
  res_en_opt[[i]]$alltrain$par_set = ll[i]
  res_en_opt_alpha[[i]]$alltest$par_set = ll[i]
  res_en_opt_alpha[[i]]$alltrain$par_set = ll[i]
  res_en_opt_n_cpg[[i]]$alltest$par_set = ll[i]
  res_en_opt_n_cpg[[i]]$alltrain$par_set = ll[i]
  res_en_05[[i]]$alltest$par_set = ll[i]
  res_en_05[[i]]$alltrain$par_set = ll[i]
  
  res_rf_opt[[i]]$alltest$par_set = ll[i]
  res_rf_opt[[i]]$alltrain$par_set = ll[i]
  res_rf_opt_mtry[[i]]$alltest$par_set = ll[i]
  res_rf_opt_mtry[[i]]$alltrain$par_set = ll[i]
  res_rf_opt_n_cpg[[i]]$alltest$par_set = ll[i]
  res_rf_opt_n_cpg[[i]]$alltrain$par_set = ll[i]
  res_rf_05[[i]]$alltest$par_set = ll[i]
  res_rf_05[[i]]$alltrain$par_set = ll[i]
}

load("output/pr_model_output.Rdata")

save(list = ls(pattern = "^res.*"), pars, file = "output/model_output.Rdata")

# 
# 
# 
# load("Worked_examples/full_model_output_par_1.Rdata")
# res1 = res
# load("Worked_examples/model_output_1.Rdata")
# res2 = res
# 
# metric_to_keep = names(table(res2[[1]]$alltest$metric))
# metric_to_keep = str_remove_all(metric_to_keep, "_0")
# for(i in 1:length(res1)){
#   res1[[i]]$alltest = res1[[i]]$alltest |> filter(metric %in% metric_to_keep)
#   res1[[i]]$alltrain = res1[[i]]$alltrain |> filter(metric %in% metric_to_keep)
#   res1[[i]]$alltest$par_set = i
#   res1[[i]]$alltrain$par_set = i
#   res1[[i]]$alltest = res1[[i]]$alltest |> mutate(metric = case_when(
#     str_detect(metric, "pheno") ~ paste0(metric, "_0"),
#     TRUE ~ metric))
#   res1[[i]]$alltrain = res1[[i]]$alltrain |> mutate(metric = case_when(
#     str_detect(metric, "pheno") ~ paste0(metric, "_0"),
#     TRUE ~ metric))
#   res1[[i]]$alltest  = res1[[i]]$alltest |> group_by(metric) |> mutate(nn = row_number()) |> ungroup() |> 
#     filter((str_detect(metric, "pheno") & nn > 426)|(!str_detect(metric, "pheno"))) |> select(-nn)
#   res1[[i]]$alltrain  = res1[[i]]$alltrain |> group_by(metric) |> mutate(nn = row_number()) |> ungroup() |> 
#     filter((str_detect(metric, "pheno") & nn > 426)|(!str_detect(metric, "pheno"))) |> select(-nn)
# }
# 
# for(i in 1:length(res2)){
#   res2[[i]]$alltest$par_set = i + 6
#   res2[[i]]$alltrain$par_set = i + 6
# }
# 
# xx = res1[[1]]$alltest
# res = c(res1,res2)
# rm(res1, res2)
