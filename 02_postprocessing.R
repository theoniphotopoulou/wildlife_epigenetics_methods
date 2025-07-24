library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

par_set_values = 1:12

pars = expand.grid(use_weights = c(FALSE, TRUE), controls_only = c(FALSE, TRUE), reduce_size = c(FALSE, TRUE), 
                   loglin_age = c("lin", "loglin"), loglin_sexm = 15)
pars = pars |> filter(!(controls_only & reduce_size)) 
pars$par_set = par_set_values

###############
# Elastic net #
###############

load("output/en_model_output.Rdata")
alltest = do.call(rbind, purrr::map(res_en, "alltest"))
alltrain = do.call(rbind, purrr::map(res_en, "alltrain"))
allval = do.call(rbind, purrr::map(res_en, "allval"))

# maximum number of n_cpg sites, used to select default models
max_ncpg = max(pars_en$n_cpg)

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
  filter(n_cpg == max_ncpg) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age) |> 
  slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age)) |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig) |> 
  pull(par_set)

res_en_opt = res_en[best_parsets_en]
res_en_opt_alpha = res_en[best_parsets_en_alpha]
res_en_opt_n_cpg = res_en[best_parsets_en_n_cpg]
id_05 = pars_en |> filter(en_alpha == 0.5, n_cpg == max_ncpg) |> pull(par_set)
res_en_05 = res_en[id_05]

# Element i in the list res_en_xxx was created with pars$par_set == i
# For later merging with parameter combos, set par_set = i
# Note: this replaces an existing par_set column in train, val, test, but
# this refers to the "pars_en" combination of parameters, which is used above
# but no longer needed

if(length(res_en_opt) != length(par_set_values)) stop("Check results objects, lengths don't match")

for(i in 1:length(par_set_values)){
  res_en_opt[[i]]$alltest$par_set = par_set_values[i]
  res_en_opt[[i]]$alltrain$par_set = par_set_values[i]
  
  res_en_opt_alpha[[i]]$alltest$par_set = par_set_values[i]
  res_en_opt_alpha[[i]]$alltrain$par_set = par_set_values[i]
  
  res_en_opt_n_cpg[[i]]$alltest$par_set = par_set_values[i]
  res_en_opt_n_cpg[[i]]$alltrain$par_set = par_set_values[i]
  
  res_en_05[[i]]$alltest$par_set = par_set_values[i]
  res_en_05[[i]]$alltrain$par_set = par_set_values[i]
}

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
  filter(n_cpg == max_ncpg) |> 
  arrange(use_weights, controls_only, reduce_size, loglin_age) |> 
  slice_min(mae, by = c(use_weights, controls_only, reduce_size, loglin_age)) |> 
  left_join(pars |> rename(par_set_orig = par_set), by = join_by(use_weights, controls_only, reduce_size, loglin_age)) |> 
  arrange(par_set_orig) |> 
  pull(par_set)

res_rf_opt = res_rf[best_parsets_rf]
res_rf_opt_mtry = res_rf[best_parsets_rf_mtry]
res_rf_opt_n_cpg = res_rf[best_parsets_rf_n_cpg]
id_05 = pars_rf |> filter(rf_mtry == 1/3, n_cpg == max_ncpg) |> pull(par_set)
res_rf_05 = res_rf[id_05]

if(length(res_rf_opt) != length(par_set_values)) stop("Check results objects, lengths don't match")

for(i in 1:length(res_en_opt)){
  res_rf_opt[[i]]$alltest$par_set = par_set_values[i]
  res_rf_opt[[i]]$alltrain$par_set = par_set_values[i]
  
  res_rf_opt_mtry[[i]]$alltest$par_set = par_set_values[i]
  res_rf_opt_mtry[[i]]$alltrain$par_set = par_set_values[i]
  
  res_rf_opt_n_cpg[[i]]$alltest$par_set = par_set_values[i]
  res_rf_opt_n_cpg[[i]]$alltrain$par_set = par_set_values[i]
  
  res_rf_05[[i]]$alltest$par_set = par_set_values[i]
  res_rf_05[[i]]$alltrain$par_set = par_set_values[i]
}


###############
# Probage #
###############

load("output/pr_model_output.Rdata")

save(list = ls(pattern = "^res.*"), pars, file = "output/model_output.Rdata")
