library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(mgcv)
library(scico)
library(scales)
library(patchwork)
library(glue)
library(kableExtra)
source("age-transformations-zoller.R")

load("output/model_output.Rdata")
load("output/dolphin_data.Rdata")

change_metric = function(res, ensemble = FALSE){
  
  if(ensemble){
    
    res = res |> 
      filter(metric %in% c("Unpub20xx", "Universal", "Robeck2021", 
                           "fgec_rf05_31192", "fgec_en05_31192", 
                           "fgec_rfopt_31192", "fgec_enopt_31192", 
                           "ens_all", "ens_lm",
      #                   "pheno_gompertz_survival_3650_31192_0", "dnampheno_gompertz_survival_3650_31192_0",
                           "probage_mle_bias_1000", "probage_mle_acc_1000"))
    
    res = res |> 
      mutate(metric = factor(metric, 
                             levels = c("Robeck2021", "Unpub20xx", "Universal",
                                        "fgec_rf05_31192", "fgec_en05_31192", 
                                        "fgec_rfopt_31192", "fgec_enopt_31192", 
                                        "ens_all", "ens_lm",
       #                                 "pheno_gompertz_survival_3650_31192_0", "dnampheno_gompertz_survival_3650_31192_0",
                                        "probage_mle_acc_1000", "probage_mle_bias_1000"),
                             labels = c("BND", "Cetaceans", "Mammals", 
                                        "Random forest", "Elastic net", 
                                        "Random forest opt", "Elastic net opt",
                                        "Ensemble", "Ensemble opt", 
      #                                  "PhenoAge", "DNAm-PhenoAge",
                                        "Probage Acc", "Probage Bias")))
    
  } else {
    
    res = res |> 
      filter(metric %in% c("Unpub20xx", "Universal", "Robeck2021", 
                           "fgec_rf05_31192", "fgec_en05_31192", 
                           "fgec_rfopt_31192", "fgec_enopt_31192", 
    #                       "pheno_gompertz_survival_3650_31192_0", "dnampheno_gompertz_survival_3650_31192_0",
                           "probage_mle_bias_1000", "probage_mle_acc_1000"))
    
    res = res |> 
      mutate(metric = factor(metric, 
                             levels = c("Robeck2021", "Unpub20xx", "Universal",
                                        "fgec_rf05_31192", "fgec_en05_31192", 
                                        "fgec_rfopt_31192", "fgec_enopt_31192", 
    #                                    "pheno_gompertz_survival_3650_31192_0", "dnampheno_gompertz_survival_3650_31192_0",
                                        "probage_mle_acc_1000", "probage_mle_bias_1000"),
                             labels = c("BND", "Cetaceans", "Mammals", 
                                        "Random forest", "Elastic net", 
                                        "Random forest opt", "Elastic net opt",
    #                                    "PhenoAge", "DNAm-PhenoAge",
                                        "Probage Acc", "Probage Bias")))
  }
  
  res
}


calculate_r <- function(p, df) {
  t <- qt(1 - p / 2, df = df)  # Convert p-value to t-statistic
  if (df + t^2 <= 0) {
    stop("Invalid input: denominator must be positive.")
  }
  r <- t / sqrt(df + t^2)
  return(r)  # Return both positive and negative roots
}

select = dplyr::select

rm(health_all, health_survmod, x, vars_to_keep, y)
health = health |> select(-control) |> 
  left_join(sample_info |> select(DOLPHIN.ID, control), by = join_by(DOLPHIN.ID))
health <- health |> 
  mutate(BMI = weight / (length/100)^2,
         left_lung_disease = (Left.Lung.Score %in% c("Moderate", "Severe")),
         lung_disease = (Left.Lung.Score %in% c("Moderate", "Severe")) | (Right.Lung.Score %in% c("Moderate", "Severe")))
health = health |> select(DOLPHIN.ID, Site_Field.ID, sex, control, pregnant, BMI, left_lung_disease, lung_disease, ALT, AST, GGT,
                          VESOP_surv, contains("_flag"))

# pars

pars = pars |> 
  mutate(uw = ifelse(use_weights, "Wtd", "Un"),
         con_rs = case_when(
           !controls_only & !reduce_size ~ "All", 
           controls_only ~ "Ctl",
           !controls_only & reduce_size ~ "Sub"),
         lla = case_when(
           loglin_age == "lin" ~ "Age",
           loglin_age == "log" ~ "L",
           loglin_age == "loglin" ~ "LL")) |> 
  mutate(lbl = paste0(lla, "\n", con_rs, "\n", uw, "\n")) |> 
  select(-uw, -con_rs, -lla)

# comb_data_dnam_trim |> filter(!is.na(BMI), !is.na(left_lung_disease)) |> 
#   select(BMI, left_lung_disease, contains("residual")) |> cor() |> round(2)

# alltest = do.call(rbind, purrr::map(res, "alltest"))
# missing_cols = setdiff(names(alltest_en_05), names(alltest))
# missing_cols = setdiff(names(alltest_rf_05), names(alltest))
# alltest[missing_cols] = NA
# alltest = alltest %>% select(all_of(names(alltest_en_05)))
alltest_en_05 = do.call(rbind, purrr::map(res_en_05, "alltest")) |> mutate(metric = "fgec_en05_31192")
alltest_en_opt = do.call(rbind, purrr::map(res_en_opt, "alltest")) |> mutate(metric = "fgec_enopt_31192")
alltest_rf_05 = do.call(rbind, purrr::map(res_rf_05, "alltest")) |> mutate(metric = "fgec_rf05_31192")
alltest_rf_opt = do.call(rbind, purrr::map(res_rf_opt, "alltest")) |> mutate(metric = "fgec_rfopt_31192")
alltest_pheno = do.call(rbind, purrr::map(res_ph, "alltest"))
alltest_probage = do.call(rbind, purrr::map(res_pr, "alltest"))
alltest = alltest_en_05 |> 
  rbind.data.frame(alltest_en_opt) |> 
  rbind.data.frame(alltest_rf_05) |> 
  rbind.data.frame(alltest_rf_opt) |> 
  rbind.data.frame(alltest_pheno) |> 
  rbind.data.frame(alltest_probage)

alltrain_en_05 = do.call(rbind, purrr::map(res_en_05, "alltrain")) |> mutate(metric = "fgec_en05_31192")
alltrain_en_opt = do.call(rbind, purrr::map(res_en_opt, "alltrain")) |> mutate(metric = "fgec_enopt_31192")
alltrain_rf_05 = do.call(rbind, purrr::map(res_rf_05, "alltrain")) |> mutate(metric = "fgec_rf05_31192")
alltrain_rf_opt = do.call(rbind, purrr::map(res_rf_opt, "alltrain")) |> mutate(metric = "fgec_rfopt_31192")
alltrain_pheno = do.call(rbind, purrr::map(res_ph, "alltrain"))
alltrain_probage = do.call(rbind, purrr::map(res_pr, "alltrain"))
alltrain = alltrain_en_05 |> 
  rbind.data.frame(alltrain_en_opt) |> 
  rbind.data.frame(alltrain_rf_05) |> 
  rbind.data.frame(alltrain_rf_opt) |> 
  rbind.data.frame(alltrain_pheno) |> 
  rbind.data.frame(alltrain_probage)

# alltest$par_set = rep(1:12, each = nrow(alltest)/12)
# alltrain$par_set = rep(1:12, each = nrow(alltrain)/12)

# m0, m1 output 463 preds 
# m2-m13 output 926, 463 for each of phenoage and dnamphenoage
# m14 out 926, 463 for each of acc and bias

# 2 * 463 + 13 * 926 = 12964 rows per parameter set
# 12 * 12964 = 155568 rows in total

alltrain = alltrain |> rename_with(.fn = ~ paste0("train", .), .cols = starts_with("preds"))
alltest = alltest |> rename_with(.fn = ~ paste0("test", .), .cols = starts_with("preds"))

allpreds = alltrain |> left_join(alltest |> select(-age), by = join_by(DOLPHIN.ID, par_set, metric))

allpreds = left_join(allpreds, pars, by = join_by(par_set))

# calculate predicted age from folds 
allpreds = allpreds |> 
  mutate(testpredage = rowMeans(pick(starts_with("testpreds")), na.rm = TRUE)) |> 
  mutate(testpredage = ifelse(is.nan(testpredage), NA, testpredage)) |> 
  mutate(trainpredage = rowMeans(pick(starts_with("trainpreds")), na.rm = TRUE)) |> 
  mutate(trainpredage = ifelse(is.nan(trainpredage), NA, trainpredage)) 

# AARs are differences (positive for "unhealthy" [epi age > chron age])
# for probage acc and bias can be directly interpreted as AARs
allpreds = allpreds |> 
  mutate(aar_diff_test = case_when(
    str_detect(metric, "probage") ~ testpredage, 
    TRUE ~ testpredage - age),
    aar_diff_train = case_when(
      str_detect(metric, "probage") ~ trainpredage, 
      TRUE ~ trainpredage - age))

# checks
xx = allpreds |> filter(metric == "dnampheno_gompertz_survival_3650_31192_0")
table(xx$loglin_age)

# drop log
allpreds = allpreds |> filter(loglin_age != "log")

# AARs are residuals from linear model 
allpreds = allpreds |> 
  group_by(metric, par_set) |> 
  mutate(aar_lm_test = case_when(
    str_detect(metric, "probage") ~ testpredage, 
    TRUE ~ as.numeric(predict(lm(age ~ testpredage, data = cur_data(), na.action = na.exclude)) - age))) |> 
  # don't end up using these so no adjustment for probage
  mutate(aar_lm_train_v1 = as.numeric(predict(lm(age ~ trainpredage, data = cur_data(), na.action = na.exclude), 
                                              newdata = data.frame(trainpredage = cur_data() |> pull(testpredage))) - age)) |> 
  mutate(aar1 = as.numeric(predict(lm(age ~ trainpreds1, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds1 = cur_data() |> pull(testpreds1))) - age),
         aar2 = as.numeric(predict(lm(age ~ trainpreds2, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds2 = cur_data() |> pull(testpreds2))) - age),
         aar3 = as.numeric(predict(lm(age ~ trainpreds3, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds3 = cur_data() |> pull(testpreds3))) - age),
         aar4 = as.numeric(predict(lm(age ~ trainpreds4, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds4 = cur_data() |> pull(testpreds4))) - age),
         aar5 = as.numeric(predict(lm(age ~ trainpreds5, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds5 = cur_data() |> pull(testpreds5))) - age),
         aar6 = as.numeric(predict(lm(age ~ trainpreds6, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds6 = cur_data() |> pull(testpreds6))) - age),
         aar7 = as.numeric(predict(lm(age ~ trainpreds7, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds7 = cur_data() |> pull(testpreds7))) - age),
         aar8 = as.numeric(predict(lm(age ~ trainpreds8, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds8 = cur_data() |> pull(testpreds8))) - age),
         aar9 = as.numeric(predict(lm(age ~ trainpreds9, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds9 = cur_data() |> pull(testpreds9))) - age),
         aar10 = as.numeric(predict(lm(age ~ trainpreds10, data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds10 = cur_data() |> pull(testpreds10))) - age)) |>
  ungroup() |> 
  mutate(aar_lm_train_v2 = rowMeans(pick(starts_with("aar[0-9]")), na.rm = TRUE)) |> 
  mutate(aar_lm_train_v2 = ifelse(is.nan(aar_lm_train_v2), NA, aar_lm_train_v2)) 

# in one case all trainpreds are constants and gam's fail - jitter these.
allpreds = allpreds |> 
  mutate(rn = runif(n(),-0.0001,0.0001), across(trainpreds1:trainpreds5, ~ .x + rn)) |> 
  mutate(rn = runif(n(),-0.0001,0.0001), across(testpreds1:testpreds5, ~ .x + rn)) |> 
  select(-rn)

# AARs are residuals from gam 
allpreds = allpreds |> 
  group_by(metric, par_set) |> 
  mutate(aar_gam_test = case_when(
    str_detect(metric, "probage") ~ testpredage, 
    TRUE ~ as.numeric(predict(gam(age ~ s(testpredage, k = 4), data = cur_data(), na.action = na.exclude)) - age))) |> 
  # don't end up using these so no adjustment for probage
  mutate(aar_gam_train_v1 = as.numeric(predict(gam(age ~ s(trainpredage,  k = 4), data = cur_data(), na.action = na.exclude), 
                                               newdata = data.frame(trainpredage = cur_data() |> pull(testpredage))) - age)) |> 
  mutate(aar1 = as.numeric(predict(gam(age ~ s(trainpreds1, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds1 = cur_data() |> pull(testpreds1))) - age),
         aar2 = as.numeric(predict(gam(age ~ s(trainpreds2, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds2 = cur_data() |> pull(testpreds2))) - age),
         aar3 = as.numeric(predict(gam(age ~ s(trainpreds3, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds3 = cur_data() |> pull(testpreds3))) - age),
         aar4 = as.numeric(predict(gam(age ~ s(trainpreds4, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds4 = cur_data() |> pull(testpreds4))) - age),
         aar5 = as.numeric(predict(gam(age ~ s(trainpreds5, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds5 = cur_data() |> pull(testpreds5))) - age),
         # aar6 = as.numeric(predict(gam(age ~ s(trainpreds6, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds6 = cur_data() |> pull(testpreds6))) - age),
         # aar7 = as.numeric(predict(gam(age ~ s(trainpreds7, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds7 = cur_data() |> pull(testpreds7))) - age),
         # aar8 = as.numeric(predict(gam(age ~ s(trainpreds8, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds8 = cur_data() |> pull(testpreds8))) - age),
         # aar9 = as.numeric(predict(gam(age ~ s(trainpreds9, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds9 = cur_data() |> pull(testpreds9))) - age),
         aar10 = as.numeric(predict(gam(age ~ s(trainpreds10, k = 4), data = cur_data(), na.action = na.exclude), newdata = data.frame(trainpreds10 = cur_data() |> pull(testpreds10))) - age)) |>
  ungroup() |> 
  mutate(aar_gam_train_v2 = rowMeans(pick(starts_with("aar[0-9]")), na.rm = TRUE)) |> 
  mutate(aar_gam_train_v2 = ifelse(is.nan(aar_gam_train_v2), NA, aar_gam_train_v2)) 

allpreds = allpreds |> select(-starts_with("aar[0-9]"))
allpreds = allpreds |> select(-starts_with("testpreds[0-9]")) |> select(-starts_with("trainpreds[0-9]"))

# add in existing clocks

load("output/existing_clocks.Rdata")
existing_clocks = ys.output |> mutate(DOLPHIN.ID = as.character(DOLPHIN.ID))
colnames(existing_clocks)
existing_clocks = existing_clocks |> select(DOLPHIN.ID, age, 
                                            Robeck2021 = DNAmAge.BottlenoseSkinAge.LogLinear2,
                                            Unpub20xx = DNAmAge.CetaceanTursiopsTruncatusSkinAge.LogLinear2,
                                            Universal = DNAmAge.Universal3_Age.LogLinearRelAdult)
dd = expand.grid(par_set = c(1:12), DOLPHIN.ID = existing_clocks$DOLPHIN.ID, stringsAsFactors = FALSE)
dd = left_join(dd, pars, by = join_by(par_set))
dd = left_join(dd, existing_clocks, by = join_by(DOLPHIN.ID))

existing_clocks = dd |> 
  pivot_longer(c("Robeck2021", "Unpub20xx", "Universal"),
               names_to = "metric", values_to = "testpredage")

# AARs (none for training cos no training)
existing_clocks = existing_clocks |> 
  mutate(aar_diff_test = (testpredage - age)) |> 
  group_by(metric, par_set) |> 
  mutate(aar_lm_test = as.numeric(predict(lm(age ~ testpredage, data = cur_data(), na.action = na.exclude)) - age)) |> 
  mutate(aar_gam_test = as.numeric(predict(gam(age ~ s(testpredage, k = 4), data = cur_data(), na.action = na.exclude)) - age)) |> 
  ungroup() |> 
  mutate(trainpredage = NA) |> 
  mutate(aar_diff_train = NA, aar_lm_train_v1 = NA, aar_lm_train_v2 = NA, aar_gam_train_v1 = NA, aar_gam_train_v2 = NA) 

colnames(existing_clocks)
colnames(allpreds)

allpreds = bind_rows(allpreds, existing_clocks) 

# ensembles

table(allpreds$metric)

# ensemble names for unoptimised version
ens_names = c("fgec_enopt_31192", "fgec_rfopt_31192", "Robeck2021", "Unpub20xx", "Universal") 
ensembles = allpreds |> 
  filter(metric %in% ens_names) |> 
  select(DOLPHIN.ID, age, metric, lbl, testpredage) |> 
  pivot_wider(
    names_from = metric,
    values_from = c(testpredage)
  ) |> 
  group_by(lbl) |>
  rowwise() |> 
  mutate(ens_all = mean(c_across(fgec_enopt_31192:Universal), na.rm = TRUE)) |> 
  mutate(ens_existing = mean(c_across(Robeck2021:Universal), na.rm = TRUE)) |> 
  ungroup() |> 
  select(DOLPHIN.ID, age, lbl, ens_existing, ens_all)

# estimate best coefficients
ens_names = c("fgec_enopt_31192", "fgec_rfopt_31192", "Robeck2021", "Unpub20xx", "Universal") 
fml_ens = as.formula(paste0("age ~ ", paste0("trainpredage_", ens_names, collapse = ' + ')))
ensembles_trained = allpreds |> 
  filter(metric %in% ens_names) |> 
  select(DOLPHIN.ID, age, metric, lbl, testpredage, trainpredage) |> 
  pivot_wider(
    names_from = metric,
    values_from = c(trainpredage, testpredage)
  ) |> 
  mutate(trainpredage_Robeck2021 = testpredage_Robeck2021, 
         trainpredage_Unpub20xx = testpredage_Unpub20xx,
         trainpredage_Universal = testpredage_Universal) |> 
  group_by(lbl) |>
  mutate(ens_lm = as.numeric(predict(lm(fml_ens, data = cur_data(), na.action = na.exclude), 
                                     newdata = cur_data() |> select(starts_with("test")) |> 
                                       rename_with(~ gsub("test", "train", .), starts_with("test"))))) |> 
  ungroup() |> 
  select(DOLPHIN.ID, lbl, ens_lm)

ensembles$ens_lm = ensembles_trained$ens_lm

ensembles = ensembles |> 
  pivot_longer(starts_with("ens"),
               names_to = "metric", values_to = "testpredage")

ensembles = ensembles |> left_join(pars, by = join_by(lbl))

# AARs (none for training cos no training)
ensembles = ensembles |> 
  mutate(aar_diff_test = (testpredage - age)) |> 
  group_by(metric, par_set) |> 
  mutate(aar_lm_test = as.numeric(predict(lm(age ~ testpredage, data = cur_data(), na.action = na.exclude)) - age)) |> 
  mutate(aar_gam_test = as.numeric(predict(gam(age ~ s(testpredage, k = 4), data = cur_data(), na.action = na.exclude)) - age)) |> 
  ungroup() |> 
  mutate(trainpredage = NA) |> 
  mutate(aar_diff_train = NA, aar_lm_train_v1 = NA, aar_lm_train_v2 = NA, aar_gam_train_v1 = NA, aar_gam_train_v2 = NA) 

colnames(ensembles)
colnames(allpreds)

allpreds = bind_rows(allpreds, ensembles)         

#

allpreds = allpreds |> 
  mutate(rb_diff_test = aar_diff_test / age,
         rb_lm_test = aar_lm_test / age,
         rb_gam_test = aar_gam_test / age,
         rb_lm_test_v1 = aar_lm_train_v1 / age,
         rb_gam_test_v1 = aar_gam_train_v1 / age)

allpreds = allpreds |> 
  mutate(age_class = case_when(
    age <= 2 ~ "Calf",
    age <= 15 ~ "Subadult",
    age <= 25 ~ "Adult",
    age > 25 ~ "Older"
  )) |> 
  left_join(health, by = join_by(DOLPHIN.ID))

# model results
tableres = allpreds |> 
  mutate(age_class = factor(age_class, levels = c("Calf", "Subadult", "Adult", "Older"))) |> 
  group_by(metric, par_set, lbl, age_class) |> 
  summarize(mae_test = list(round(as.numeric(quantile(abs(aar_diff_test), prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 1)),
            mrb_test = list(round(100 * as.numeric(quantile(rb_diff_test, prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 0))) |> 
  mutate(mae_test = purrr::map_chr(mae_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})")),
         mrb_test = purrr::map_chr(mrb_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})"))) |> 
  ungroup()

tableres_all = allpreds |> 
  mutate(age_class = "All") |> 
  group_by(metric, par_set, lbl, age_class) |> 
  summarize(mae_test = list(round(as.numeric(quantile(abs(aar_diff_test), prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 1)),
            mrb_test = list(round(100 * as.numeric(quantile(rb_diff_test, prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 0))) |> 
  mutate(mae_test = purrr::map_chr(mae_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})")),
         mrb_test = purrr::map_chr(mrb_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})"))) |> 
  ungroup()

tableres = rbind(tableres, tableres_all)

# 10 or 16
tablepar = 10
tableres1 = tableres |> filter(par_set == tablepar) 
tableres1 = tableres1 |> pivot_wider(id_cols = metric, names_from = age_class, values_from = mae_test, names_sort = FALSE)
tableres1 = tableres1 |> change_metric(ensemble = TRUE) |> arrange(metric) #|> filter(!str_detect(metric, "Pheno"))

tableres2 = tableres |> filter(par_set == tablepar) 
tableres2 = tableres2 |> pivot_wider(id_cols = metric, names_from = age_class, values_from = mrb_test, names_sort = FALSE)
tableres2 = tableres2 |> change_metric(ensemble = TRUE) |> arrange(metric) #|> filter(!str_detect(metric, "Pheno"))

tableres12 = rbind(tableres1, tableres2)
tableres12
kable(tableres12, format = "latex", booktabs = TRUE)

# checks
xx = allpreds |> change_metric(ensemble = TRUE)
xx |> filter(par_set == 10) |> 
  group_by(metric) |> 
  summarize(min= min(testpredage,na.rm=TRUE),
            max = max(testpredage, na.rm = TRUE))

# removing log transformation
allpreds = allpreds |> filter(loglin_age != "log")

res = allpreds |> 
  group_by(metric, lbl) |> 
  summarize(cor_testpredage = cor.test(age, testpredage, use = "na.or.complete")$estimate,
            cor_trainpredage = cor(age, trainpredage, use = "na.or.complete"),
            df_testpredage = cor.test(age, testpredage, use = "na.or.complete")$parameter,
            p_testpredage = cor.test(age, testpredage, use = "na.or.complete")$p.value,
            me_test = median(aar_diff_test, na.rm = TRUE),
            mlme_test = median(aar_lm_test, na.rm = TRUE),
            mgame_test = median(aar_gam_test, na.rm = TRUE), 
            mae_test = median(abs(aar_diff_test), na.rm = TRUE),
            malme_test = median(abs(aar_lm_test), na.rm = TRUE),
            magame_test = median(abs(aar_gam_test), na.rm = TRUE),
            mrb_diff_test = 100 * median(rb_diff_test, na.rm = TRUE),
            mrb_lm_test = 100 * median(rb_lm_test, na.rm = TRUE),
            mrb_gam_test = 100 * median(rb_gam_test, na.rm = TRUE),
            malme_test_v1 = median(abs(aar_lm_train_v1), na.rm = TRUE),
            magame_test_v1 = median(abs(aar_gam_train_v1), na.rm = TRUE),
            mrb_lm_test_v1 = 100 * median(rb_lm_test_v1, na.rm = TRUE),
            mrb_gam_test_v1 = 100 * median(rb_gam_test_v1, na.rm = TRUE)) |> 
  ungroup()

res_2g = res |> filter(str_detect(metric, "pheno"))

res = res |> change_metric(ensemble = TRUE)

# color clock names by the kind of clocks they are
## for figures where all clocks are shown
# with pheno
# ycols = c(rep("purple", 2), rep("red", 2), rep("lightblue", 2), rep("blue", 4), rep("black", 3))
# no pheno
ycols = c(rep("purple", 2), rep("lightblue", 2), rep("blue", 4), rep("black", 3))
## for figures where pheno and probage clocks not shown
ycols23 = c(rep("lightblue", 2), rep("blue", 4), rep("black", 4))

# Create a summarized dataset where we merge tiles for a specific row
# Identify the row where all values are the same
row_to_merge <- c("BND", "Cetaceans", "Mammals")  # Change this to the desired y level

# Extract the unique fill value for that row
merged_row <- res %>%
  filter(metric %in% row_to_merge) %>%
  group_by(metric) |> 
  summarize(metric = unique(metric), lbl = first(lbl), 
            cor_testpredage = first(cor_testpredage),
            cor_testpredage = first(cor_testpredage)) 

xt = res |> group_by(metric) |> summarize(minx = min(cor_testpredage, na.rm = TRUE),
                                               maxx = max(cor_testpredage, na.rm = TRUE))

xt = res |> group_by(metric, lbl) |> summarize(minx = min(cor_testpredage, na.rm = TRUE),
                                          maxx = max(cor_testpredage, na.rm = TRUE))

xt = res |> filter(metric %in% c("Elastic net", "Elastic net opt", "Random forest", "Random forest opt")) |>
  mutate(optim = ifelse(str_detect(metric, "opt"), "Optimized", "Unoptimized")) |> 
  mutate(model = str_extract(metric, "\\w+")) |> 
  select(lbl, model, optim, value = mae_test) |> 
  pivot_wider(names_from = optim, values_from = value) |> 
  mutate(opt_imp = (Optimized - Unoptimized)) |> 
  mutate(opt_imp_rel = 100 * (Optimized - Unoptimized)/Unoptimized) |> 
  arrange(model, desc(opt_imp_rel))

xt |> group_by(model) |> summarize(mdf1 = mean(opt_imp), mf1 = mean(Optimized))

f0 = res |> 
  mutate(cor_testpredage = ifelse(str_detect(metric, "[Pp]robage"), -cor_testpredage, cor_testpredage)) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_testpredage)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = merged_row, aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_testpredage), width = Inf, height = 1, color = "black") +
  scale_fill_distiller(palette = "YlOrBr", direction = 1, limits = c(0, 1)) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(a) Correlation between observed and DNAm-predicted age") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols)) 

f0

# better for controls than not
res_with_control = allpreds |> group_by(control, metric, lbl) |> 
  summarize(cor_testpredage = cor.test(age, testpredage, use = "na.or.complete")$estimate,
            cor_trainpredage = cor(age, trainpredage, use = "na.or.complete"),
            df_testpredage = cor.test(age, testpredage, use = "na.or.complete")$parameter,
            p_testpredage = cor.test(age, testpredage, use = "na.or.complete")$p.value,
            me_test = median(aar_diff_test, na.rm = TRUE),
            mlme_test = median(aar_lm_test, na.rm = TRUE),
            mgame_test = median(aar_gam_test, na.rm = TRUE), 
            mae_test = median(abs(aar_diff_test), na.rm = TRUE),
            malme_test = median(abs(aar_lm_test), na.rm = TRUE),
            magame_test = median(abs(aar_gam_test), na.rm = TRUE),
            mrb_diff_test = 100 * median(rb_diff_test, na.rm = TRUE),
            mrb_lm_test = 100 * median(rb_lm_test, na.rm = TRUE),
            mrb_gam_test = 100 * median(rb_gam_test, na.rm = TRUE)) |> 
  ungroup() |> 
  mutate(cor_testpredage = ifelse(str_detect(metric, "[Pp]robage"), -cor_testpredage, cor_testpredage)) 

res_with_control = res_with_control |> change_metric(ensemble = TRUE)

f1 = res_with_control |> 
  pivot_wider(id_cols = c("metric", "lbl"), names_from = control, 
              values_from = cor_testpredage, names_prefix = "control") |> 
  mutate(cordiff = control1 - control0) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cordiff)) + 
  geom_tile(colour = "black") + 
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5, 0.5)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) Change in correlation in control vs exposed group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f1

f1control = res_with_control |> 
  pivot_wider(id_cols = c("metric", "lbl"), names_from = control, 
              values_from = cor_testpredage, names_prefix = "control") |> 
  mutate(cordiff = control1 - control0) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = control1)) + 
  geom_tile(colour = "black") + 
  scale_fill_distiller(palette = "YlOrBr", direction = 1, limits = c(0, 1)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) Correlation in control group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f1exposed = res_with_control |> 
  pivot_wider(id_cols = c("metric", "lbl"), names_from = control, 
              values_from = cor_testpredage, names_prefix = "control") |> 
  mutate(cordiff = control1 - control0) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = control0)) + 
  geom_tile(colour = "black") + 
  scale_fill_distiller(palette = "YlOrBr", direction = 1, limits = c(0, 1)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) Correlation in exposed group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f1control
f1exposed

#

# scatters

fscat = allpreds |> change_metric(ensemble = TRUE) |> 
  mutate(control = factor(control, levels = c(0, 1), labels = c("Exposed population", "Healthy population"))) |> 
  filter(par_set == 10 & !str_detect(metric, "Probage") & metric != "PhenoAge") |> 
  #filter(par_set == 16 & !str_detect(metric, "Probage")) |> 
  ggplot(aes(x = age, y = testpredage)) + 
  geom_point(aes(fill = control), colour = "black", shape=21) +
  scale_fill_manual(values = c("white", "gray")) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = FALSE, colour = "blue", linewidth = 0.7) +
  coord_fixed(xlim = c(0, 65), ylim = c(-1, 65)) +
  labs(subtitle = "(b) Observed vs DNAm-predicted ages") +
  facet_wrap(~metric, ncol = 5) +
  labs(x = "Chronological age", y = "Epigenetic age") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.margin = margin(-20, 0, 0, 0),
        strip.background = element_rect(fill = NA), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

fscat

fscat_tr_bb = allpreds |> change_metric(ensemble = FALSE) |> 
  filter((metric %in% c("Elastic net opt", "Random forest opt"))) |> 
  filter(par_set == 8 & !str_detect(metric, "Probage") & metric != "PhenoAge") |> 
  select(DOLPHIN.ID, par_set, metric, age, control, trainpreds1:trainpreds6,testpredage) |> 
  pivot_longer(cols = c(trainpreds1:trainpreds6,testpredage), values_to = "predage") |> 
  ggplot(aes(x = age, y = predage)) + 
  geom_point(colour = "black", shape=21) +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = FALSE, colour = "blue", linewidth = 0.7) +
  coord_fixed(xlim = c(0, 65), ylim = c(-1, 65)) +
  facet_grid(metric ~ name) +
  labs(x = "Chronological age", y = "Epigenetic age") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

fscat_tr_bb
ggsave("output/figs/test_train_bb.png", fscat_tr_bb, width = 9, height = 2.5, dpi = 300)


fscat_tr = allpreds |> change_metric(ensemble = TRUE) |> 
  filter(!(metric %in% c("BND", "Cetaceans", "Mammals"))) |> 
  mutate(control = factor(control, levels = c(0, 1), labels = c("Exposed population", "Healthy population"))) |> 
  filter(par_set == 8 & !str_detect(metric, "Probage") & metric != "PhenoAge") |> 
  ggplot(aes(x = age, y = trainpredage)) + 
  geom_point(aes(fill = control), colour = "black", shape=21) +
  scale_fill_manual(values = c("white", "gray")) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4), se = FALSE, colour = "blue", linewidth = 0.7) +
  coord_fixed(xlim = c(0, 65), ylim = c(-1, 65)) +
  labs(subtitle = "(c) Chronological age versus DNAm-predicted age") +
  facet_wrap(~metric, ncol = 4) +
  labs(x = "Chronological age", y = "Epigenetic age") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

fscat_tr


# RB

f2 = res |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_diff_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_diff_test = first(mrb_diff_test),
                        mrb_diff_test = first(mrb_diff_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_diff_test), width = Inf, height = 1, color = "black") +
  scico::scale_fill_scico(palette = "vik", limits = c(-25, 25)) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(a) Relative bias of difference-based AARs") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2

# better for controls than not
res_with_control_bias = res_with_control |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  pivot_wider(id_cols = c("metric", "lbl"), names_from = control, 
              #        values_from = mrb_diff_test, names_prefix = "mrb_diff_test",
              values_from = c(mrb_diff_test, mrb_lm_test),
              names_glue = "{.value}{control}") |> 
  mutate(mrbdiff = mrb_diff_test0 - mrb_diff_test1) |> 
  left_join(res |> select(metric, lbl, mrb_diff_test), by = join_by(metric, lbl)) |> 
  mutate(prop_mrbdiff1 = 100 * (mrb_diff_test - mrb_diff_test1) / mrb_diff_test,
         prop_mrbdiff0 = 100 * (mrb_diff_test - mrb_diff_test0) / mrb_diff_test)

f2control = res_with_control_bias |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_diff_test1)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_bias %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_diff_test1 = first(mrb_diff_test1),
                        mrb_diff_test1 = first(mrb_diff_test1)),
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_diff_test1), width = Inf, height = 1, color = "black") +
  scico::scale_fill_scico(palette = "vik", limits = c(-50, 50)) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(a) Relative bias of AARs in control group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2control

f2exposed = res_with_control_bias |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_diff_test0)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_bias %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_diff_test0 = first(mrb_diff_test0),
                        mrb_diff_test0 = first(mrb_diff_test0)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_diff_test0), width = Inf, height = 1, color = "black") +
  scico::scale_fill_scico(palette = "vik", limits = c(-50, 50)) +
  scale_x_discrete(expand = expansion(add = 0)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) Relative bias of AARs in exposed group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2exposed

f2control_lm = res_with_control_bias |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test1)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_bias %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_lm_test1 = first(mrb_lm_test1),
                        mrb_lm_test1 = first(mrb_lm_test1)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test1), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-25, 25)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(a) Relative bias of AARs in control group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2control_lm

f2exposed_lm = res_with_control_bias |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test0)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_bias %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_lm_test0 = first(mrb_lm_test0),
                        mrb_lm_test0 = first(mrb_lm_test0)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test0), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-25, 25)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) Relative bias of AARs in exposed group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2exposed_lm

# better for lm than diff
# res_lm_mae = allpreds |>
#   mutate(mae_lm_diff = abs(aar_diff_test) - abs(aar_lm_test),
#          lm_gam_diff = abs(aar_lm_test) - abs(aar_gam_test),
#          mae_gam_diff = abs(aar_diff_test) - abs(aar_gam_test)) |>
#   group_by(metric, lbl) |>
#   summarize(mae_lm_diff = median(mae_lm_diff, na.rm = TRUE),
#             lm_gam_diff = median(lm_gam_diff, na.rm = TRUE),
#             mae_gam_diff = median(mae_gam_diff, na.rm = TRUE)) |>
#   ungroup() |>
#   change_metric(ensemble = TRUE)

res_lm_mae = res |>
  mutate(mae_lm_diff = 100 * (mrb_diff_test - mrb_lm_test) / mrb_diff_test,
         lm_gam_diff = 100 * (mrb_lm_test - mrb_gam_test) / malme_test,
         mae_gam_diff = 100 * (mrb_diff_test - mrb_gam_test) / mae_test,
         lm_te_tr_diff = 100 * (mrb_lm_test_v1 - mrb_lm_test) / mrb_lm_test,
         gam_te_tr_diff = 100 * (mrb_gam_test_v1 - mrb_gam_test) / mrb_gam_test) |>
  ungroup() 

f2lm = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_lm_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_lm_test = first(mrb_lm_test),
                        mrb_lm_test = first(mrb_lm_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-25, 25)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) Relative bias of AARs from linear model") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2lm

f2lm_tr = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test_v1)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_lm_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_lm_test_v1 = first(mrb_lm_test_v1),
                        mrb_lm_test_v1 = first(mrb_lm_test_v1)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test_v1), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-25, 25)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) Relative bias of AARs from linear model fit to training data") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2lm_tr

f2gam = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_gam_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_lm_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_gam_test = first(mrb_gam_test),
                        mrb_gam_test = first(mrb_gam_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_gam_test), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-25, 25)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(c) Relative bias of AARs from a GAM") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2gam

f2gam_tr = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test_v1)) + 
  geom_tile(colour = "black") + 
  scico::scale_fill_scico(palette = "vik", limits = c(-25, 25)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(f) Relative bias of AARs from a GAM fit to training data") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f2gam_tr

# MAE

f3 = res |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mae_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_lm_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mae_test = first(mae_test),
                        mae_test = first(mae_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mae_test), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "acton", limits = c(1, 3.25), direction = -1) +
  # scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(1, 3.25)) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(d) MAE of difference-based AARs") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3

# better for controls than not
res_with_control_mae = res_with_control |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  pivot_wider(id_cols = c("metric", "lbl"), names_from = control, 
              #  values_from = mae_test, names_prefix = "mae_test",
              values_from = c(mae_test, malme_test),
              names_glue = "{.value}{control}") |> 
  mutate(maediff = mae_test0 - mae_test1) |> 
  left_join(res |> select(metric, lbl, mae_test), by = join_by(metric, lbl)) |> 
  mutate(prop_maediff1 = 100 * (mae_test - mae_test1) / mae_test,
         prop_maediff0 = 100 * (mae_test - mae_test0) / mae_test) 

f3control = res_with_control_mae |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mae_test1)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mae_test1 = first(mae_test1),
                        mae_test1 = first(mae_test1)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mae_test1), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "acton", limits = c(1, 4), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(c) MAE of AARs in control group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3control

f3exposed = res_with_control_mae |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mae_test0)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mae_test0 = first(mae_test0),
                        mae_test0 = first(mae_test0)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mae_test0), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "acton", limits = c(1, 4), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(d) MAE of AARs in exposed group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3exposed

f3control_lm = res_with_control_mae |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test1)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        malme_test1 = first(malme_test1),
                        malme_test1 = first(malme_test1)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test1), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "acton", limits = c(1, 3.25), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(c) MAE of AARs in control group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3control_lm

f3exposed_lm = res_with_control_mae |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test0)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        malme_test0 = first(malme_test0),
                        malme_test0 = first(malme_test0)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test0), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "acton", limits = c(1, 3.25), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(d) MAE of AARs in exposed group") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3exposed_lm

# better for lm than diff
# res_lm_mae = allpreds |>
#   mutate(mae_lm_diff = abs(aar_diff_test) - abs(aar_lm_test),
#          lm_gam_diff = abs(aar_lm_test) - abs(aar_gam_test),
#          mae_gam_diff = abs(aar_diff_test) - abs(aar_gam_test)) |>
#   group_by(metric, lbl) |>
#   summarize(mae_lm_diff = median(mae_lm_diff, na.rm = TRUE),
#             lm_gam_diff = median(lm_gam_diff, na.rm = TRUE),
#             mae_gam_diff = median(mae_gam_diff, na.rm = TRUE)) |>
#   ungroup() |>
#   change_metric(ensemble = TRUE)

res_lm_mae = res |>
  mutate(mae_lm_diff = 100 * (mae_test - malme_test) / mae_test,
         lm_gam_diff = 100 * (malme_test - magame_test) / malme_test,
         mae_gam_diff = 100 * (mae_test - magame_test) / mae_test,
         lm_te_tr_diff = 100 * (malme_test_v1 - malme_test) / malme_test,
         gam_te_tr_diff = 100 * (magame_test_v1 - magame_test) / magame_test) |>
  ungroup() 

f3lm = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_lm_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        malme_test = first(malme_test),
                        malme_test = first(malme_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "acton", limits = c(1, 3.25), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(e) MAE of AARs from linear model fit") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3lm

f3lm_tr = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test_v1)) + 
  geom_tile(colour = "black") + 
  scico::scale_fill_scico(palette = "acton", limits = c(1, 3.25), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(d) MAE of AARs from linear model fit to training data") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3lm_tr

f3gam = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = magame_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_lm_mae %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        magame_test = first(magame_test),
                        magame_test = first(magame_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = magame_test), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "acton", limits = c(1, 3.25), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(f) MAE of AARs from GAM") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3gam

f3gam_tr = res_lm_mae |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = magame_test_v1)) + 
  geom_tile(colour = "black") + 
  scico::scale_fill_scico(palette = "acton", limits = c(1, 3.25), direction = -1) +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(f) MAE of AARs from GAM fit to training data") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols23))

f3gam_tr

f2t = f2 + theme(axis.text.x = element_blank()) 
f2lmt = f2lm + theme(axis.text.x = element_blank()) 
f2gamt = f2gam  
f2gamtr = f2gam_tr 

f3t = f3 + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 
f3lmt = f3lm + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 
f3gamt = f3gam + theme(axis.text.y = element_blank()) 
f3gamtr = f3gam_tr + theme(axis.text.y = element_blank()) 

# f3lm_trt = f3lm_tr + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
# f3gam_trt = f3gam_tr + theme(axis.text.y = element_blank())

f2all = f2t + f2lmt + f2gamt & theme(legend.position = "bottom", legend.box.margin = margin(-20, 0, 0, 0),
                                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
f2all = f2all + plot_layout(ncol = 1, guides = "collect", byrow = TRUE) 
f2all

f3all = f3t + f3lmt + f3gamt & theme(legend.position = "bottom", legend.box.margin = margin(-20, 0, 0, 0),
                                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
f3all = f3all + plot_layout(ncol = 1, guides = "collect", byrow = TRUE) 
f3all

f23all = f2all | f3all
f23all

f2ct = f2control_lm + theme(axis.text.x = element_blank())
f2exp = f2exposed_lm 

f3ct = f3control_lm + theme(axis.text.x = element_blank(), axis.text.y = element_blank())
f3exp = f3exposed_lm + theme(axis.text.y = element_blank())

f4all = (f2ct + f2exp) &  theme(legend.position = "bottom", legend.box.margin = margin(-20, 0, 0, 0),
                                axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
f4all = f4all + plot_layout(ncol = 1, guides = "collect", byrow = TRUE) 
f4all

f5all = (f3ct + f3exp) &  theme(legend.position = "bottom", legend.box.margin = margin(-20, 0, 0, 0),
                                axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
f5all = f5all + plot_layout(ncol = 1, guides = "collect", byrow = TRUE) 
f5all

f45all = f4all | f5all
f45all

###

f0t = f0 + theme(axis.ticks.y = element_blank(), legend.position = "bottom") 
f1t = f1 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom")
fcor = f0t + fscat + plot_layout(ncol = 2, widths = c(3,5))
fcor

#fcor2 = (f0/f1) | fscat 
fcorp2 = fscat + fscat_tr + plot_layout(ncol = 1, widths = c(2,1))
fcorp2 <- (wrap_elements(fscat) / wrap_elements(fscat_tr)) + 
  plot_layout(heights = c(1.55, 1))
fcorp2
fcor2 = f0 | fcorp2
fcor2

ggsave("output/figs/fcorxy.png", fcor, width = 12, height = 5, dpi = 300)
ggsave("output/figs/faartype.png", f23all, width = 10, height = 7, dpi = 300)
ggsave("output/figs/faarbygroup.png", f45all, width = 10, height = 5, dpi = 300)

# ggsave("output/figs/fcor.png", fcor, width = 10, height = 5, dpi = 300)
# ggsave("output/figs/fcor2.png", fcor2, width = 10, height = 7, dpi = 300)
# ggsave("output/figs/fbias.png", f2all, width = 10, height = 5, dpi = 300)
# ggsave("output/figs/fmae.png", f3all, width = 10, height = 5, dpi = 300)
# ggsave("output/figs/fscatter.png", fscat, width = 7, height = 6, dpi = 300)

# ggsave("output/figs/fcor2.png", fcor2, width = 10, height = 7, dpi = 300)


# HEALTH

allpreds = allpreds |> 
  mutate(pregnant = case_when(
    pregnant == "No" ~ 0,
    pregnant == "Probable" ~ 1,
    pregnant == "Suspect" ~ 1,
    pregnant == "Yes" ~ 1,
    TRUE ~ NA
  ))

res_health = allpreds |> group_by(metric, lbl) |> 
  mutate(resids_pregnant_agesex = as.numeric(residuals(glm(pregnant ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_VESOP_agesex = as.numeric(residuals(glm(VESOP_surv ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_BMI05_agesex = as.numeric(residuals(glm(BMI_05_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_Anemia_agesex = as.numeric(residuals(glm(Anemia_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_Cortisol_agesex = as.numeric(residuals(glm(Cortisol_lt30_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_Neutrophils_agesex = as.numeric(residuals(glm(Neutrophils_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_TotalT3_agesex = as.numeric(residuals(glm(TotalT3_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_BMI_25_agesex = as.numeric(residuals(glm(BMI_25_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude))))  |> 
  mutate(resids_Alkaline.Phosphatase_agesex = as.numeric(residuals(glm(Alkaline.Phosphatase_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude))))  |> 
  mutate(resids_Lung_agesex = as.numeric(residuals(glm(Lung_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_Globulin_agesex = as.numeric(residuals(glm(Globulin_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  mutate(resids_Glucose_agesex = as.numeric(residuals(glm(Glucose_Flag ~ age + sex, family = "binomial", data = cur_data(), na.action = na.exclude)))) |> 
  ungroup()

res_health = res_health |> group_by(metric, lbl) |> 
  summarize(cor_diff_pregnant = cor.test(aar_diff_test, resids_pregnant_agesex, use = "na.or.complete")$estimate,
            cor_lm_pregnant = cor.test(aar_lm_test, resids_pregnant_agesex, use = "na.or.complete")$estimate,
            cor_diff_VESOP = cor.test(aar_diff_test, resids_VESOP_agesex, use = "na.or.complete")$estimate,
            cor_lm_VESOP = cor.test(aar_lm_test, resids_VESOP_agesex, use = "na.or.complete")$estimate,
            cor_diff_BMI05 = cor.test(aar_diff_test, resids_BMI05_agesex, use = "na.or.complete")$estimate,
            cor_lm_BMI05 = cor.test(aar_lm_test, resids_BMI05_agesex, use = "na.or.complete")$estimate,
            cor_diff_Lung = cor.test(aar_diff_test, resids_Lung_agesex, use = "na.or.complete")$estimate,
            cor_lm_Lung = cor.test(aar_lm_test, resids_Lung_agesex, use = "na.or.complete")$estimate,
            cor_diff_Anemia = cor.test(aar_diff_test, resids_Anemia_agesex, use = "na.or.complete")$estimate,
            cor_lm_Anemia = cor.test(aar_lm_test, resids_Anemia_agesex, use = "na.or.complete")$estimate,
            cor_diff_Cortisol = cor.test(aar_diff_test, resids_Cortisol_agesex, use = "na.or.complete")$estimate,
            cor_lm_Cortisol = cor.test(aar_lm_test, resids_Cortisol_agesex, use = "na.or.complete")$estimate,
            cor_diff_BMI25 = cor.test(aar_diff_test, resids_BMI_25_agesex, use = "na.or.complete")$estimate,
            cor_lm_BMI25 = cor.test(aar_lm_test, resids_BMI_25_agesex, use = "na.or.complete")$estimate,
            cor_diff_Neutrophils = cor.test(aar_diff_test, resids_Neutrophils_agesex, use = "na.or.complete")$estimate,
            cor_lm_Neutrophils = cor.test(aar_lm_test, resids_Neutrophils_agesex, use = "na.or.complete")$estimate,
            cor_diff_TotalT3 = cor.test(aar_diff_test, resids_TotalT3_agesex, use = "na.or.complete")$estimate,
            cor_lm_TotalT3 = cor.test(aar_lm_test, resids_TotalT3_agesex, use = "na.or.complete")$estimate,
            cor_diff_AlkPhos = cor.test(aar_diff_test, resids_Alkaline.Phosphatase_agesex, use = "na.or.complete")$estimate,
            cor_lm_AlkPhos = cor.test(aar_lm_test, resids_Alkaline.Phosphatase_agesex, use = "na.or.complete")$estimate,
            cor_diff_Globulin = cor.test(aar_diff_test, resids_Globulin_agesex, use = "na.or.complete")$estimate,
            cor_lm_Globulin = cor.test(aar_lm_test, resids_Globulin_agesex, use = "na.or.complete")$estimate,
            cor_diff_Glucose = cor.test(aar_diff_test, resids_Glucose_agesex, use = "na.or.complete")$estimate,
            cor_lm_Glucose = cor.test(aar_lm_test, resids_Glucose_agesex, use = "na.or.complete")$estimate,
            #
            df_diff_pregnant = cor.test(aar_diff_test, resids_pregnant_agesex, use = "na.or.complete")$parameter,
            df_lm_pregnant = cor.test(aar_lm_test, resids_pregnant_agesex, use = "na.or.complete")$parameter,
            df_diff_VESOP = cor.test(aar_diff_test, resids_VESOP_agesex, use = "na.or.complete")$parameter,
            df_lm_VESOP = cor.test(aar_lm_test, resids_VESOP_agesex, use = "na.or.complete")$parameter,
            df_diff_BMI05 = cor.test(aar_diff_test, resids_BMI05_agesex, use = "na.or.complete")$parameter,
            df_lm_BMI05 = cor.test(aar_lm_test, resids_BMI05_agesex, use = "na.or.complete")$parameter,
            df_diff_Lung = cor.test(aar_diff_test, resids_Lung_agesex, use = "na.or.complete")$parameter,
            df_lm_Lung = cor.test(aar_lm_test, resids_Lung_agesex, use = "na.or.complete")$parameter,
            df_diff_Anemia = cor.test(aar_diff_test, resids_Anemia_agesex, use = "na.or.complete")$parameter,
            df_lm_Anemia = cor.test(aar_lm_test, resids_Anemia_agesex, use = "na.or.complete")$parameter,
            df_diff_Cortisol = cor.test(aar_diff_test, resids_Cortisol_agesex, use = "na.or.complete")$parameter,
            df_lm_Cortisol = cor.test(aar_lm_test, resids_Cortisol_agesex, use = "na.or.complete")$parameter,
            df_diff_BMI25 = cor.test(aar_diff_test, resids_BMI_25_agesex, use = "na.or.complete")$parameter,
            df_lm_BMI25 = cor.test(aar_lm_test, resids_BMI_25_agesex, use = "na.or.complete")$parameter,
            df_diff_Neutrophils = cor.test(aar_diff_test, resids_Neutrophils_agesex, use = "na.or.complete")$parameter,
            df_lm_Neutrophils = cor.test(aar_lm_test, resids_Neutrophils_agesex, use = "na.or.complete")$parameter,
            df_diff_TotalT3 = cor.test(aar_diff_test, resids_TotalT3_agesex, use = "na.or.complete")$parameter,
            df_lm_TotalT3 = cor.test(aar_lm_test, resids_TotalT3_agesex, use = "na.or.complete")$parameter,
            df_diff_AlkPhos = cor.test(aar_diff_test, resids_Alkaline.Phosphatase_agesex, use = "na.or.complete")$parameter,
            df_lm_AlkPhos = cor.test(aar_lm_test, resids_Alkaline.Phosphatase_agesex, use = "na.or.complete")$parameter,
            df_diff_Globulin = cor.test(aar_diff_test, resids_Globulin_agesex, use = "na.or.complete")$parameter,
            df_lm_Globulin = cor.test(aar_lm_test, resids_Globulin_agesex, use = "na.or.complete")$parameter,
            df_diff_Glucose = cor.test(aar_diff_test, resids_Glucose_agesex, use = "na.or.complete")$parameter,
            df_lm_Glucose = cor.test(aar_lm_test, resids_Glucose_agesex, use = "na.or.complete")$parameter,
            #
            p_diff_pregnant = cor.test(aar_diff_test, resids_pregnant_agesex, use = "na.or.complete")$p.value,
            p_lm_pregnant = cor.test(aar_lm_test, resids_pregnant_agesex, use = "na.or.complete")$p.value,
            p_diff_VESOP = cor.test(aar_diff_test, resids_VESOP_agesex, use = "na.or.complete")$p.value,
            p_lm_VESOP = cor.test(aar_lm_test, resids_VESOP_agesex, use = "na.or.complete")$p.value,
            p_diff_BMI05 = cor.test(aar_diff_test, resids_BMI05_agesex, use = "na.or.complete")$p.value,
            p_lm_BMI05 = cor.test(aar_lm_test, resids_BMI05_agesex, use = "na.or.complete")$p.value,
            p_diff_Lung = cor.test(aar_diff_test, resids_Lung_agesex, use = "na.or.complete")$p.value,
            p_lm_Lung = cor.test(aar_lm_test, resids_Lung_agesex, use = "na.or.complete")$p.value,
            p_diff_Anemia = cor.test(aar_diff_test, resids_Anemia_agesex, use = "na.or.complete")$p.value,
            p_lm_Anemia = cor.test(aar_lm_test, resids_Anemia_agesex, use = "na.or.complete")$p.value,
            p_diff_Cortisol = cor.test(aar_diff_test, resids_Cortisol_agesex, use = "na.or.complete")$p.value,
            p_lm_Cortisol = cor.test(aar_lm_test, resids_Cortisol_agesex, use = "na.or.complete")$p.value,
            p_diff_BMI25 = cor.test(aar_diff_test, resids_BMI_25_agesex, use = "na.or.complete")$p.value,
            p_lm_BMI25 = cor.test(aar_lm_test, resids_BMI_25_agesex, use = "na.or.complete")$p.value,
            p_diff_Neutrophils = cor.test(aar_diff_test, resids_Neutrophils_agesex, use = "na.or.complete")$p.value,
            p_lm_Neutrophils = cor.test(aar_lm_test, resids_Neutrophils_agesex, use = "na.or.complete")$p.value,
            p_diff_TotalT3 = cor.test(aar_diff_test, resids_TotalT3_agesex, use = "na.or.complete")$p.value,
            p_lm_TotalT3 = cor.test(aar_lm_test, resids_TotalT3_agesex, use = "na.or.complete")$p.value,
            p_diff_AlkPhos = cor.test(aar_diff_test, resids_Alkaline.Phosphatase_agesex, use = "na.or.complete")$p.value,
            p_lm_AlkPhos = cor.test(aar_lm_test, resids_Alkaline.Phosphatase_agesex, use = "na.or.complete")$p.value,
            p_diff_Globulin = cor.test(aar_diff_test, resids_Globulin_agesex, use = "na.or.complete")$p.value,
            p_lm_Globulin = cor.test(aar_lm_test, resids_Globulin_agesex, use = "na.or.complete")$p.value,
            p_diff_Glucose = cor.test(aar_diff_test, resids_Glucose_agesex, use = "na.or.complete")$p.value,
            p_lm_Glucose = cor.test(aar_lm_test, resids_Glucose_agesex, use = "na.or.complete")$p.value) |> 
  ungroup() |> 
  change_metric(ensemble = TRUE)

res_health_long1 = res_health %>% pivot_longer(
  cols = cor_diff_pregnant:p_lm_Glucose,
  names_to = c("stat", "resid_type", "outcome"),
  names_pattern = "(.*)_(.*)_(.*)",
  values_to = "value") |> 
  filter(outcome %in% c("Lung", "Cortisol", "Neutrophils", "Anemia", "AlkPhos", "VESOP")) 

res_health_long2 = res_health_long1 |>
  filter(stat == "p") |> 
  mutate(stat = "p.adj") |> 
  mutate(value = p.adjust(value, method = "fdr")) 

res_health_long1 = rbind.data.frame(res_health_long1, res_health_long2)
rm(res_health_long2)

xt = res_health_long1 |> filter(stat == "cor") |> 
  group_by(metric, outcome, resid_type) |>
  summarize(meancor = mean((value), na.rm = TRUE)) |> ungroup() |> 
  pivot_wider(names_from = resid_type, values_from = meancor) |> 
  mutate(cordiff = (lm - diff)/ diff) |> 
  mutate(lmbiggercor = ifelse((abs(lm) > abs(diff)), 1, 0)) 

xt |> group_by(metric) |> summarize(meanlmbigger = mean(lmbiggercor))
  
res_health_wide = res_health_long1 |> 
  pivot_wider(
    names_from = c(stat, resid_type),
    names_glue = "{stat}_{resid_type}",
    values_from = value
  ) 

res_health_wide = res_health_wide |> 
  mutate(cor_diff = ifelse(p.adj_diff > 0.05, NA, cor_diff)) |> 
  mutate(cor_lm = ifelse(p.adj_lm > 0.05, NA, cor_lm)) 

res_health_wide = res_health_wide %>%
  pivot_longer(
    cols = c(starts_with("cor_"), starts_with("df_"), starts_with("p_"), starts_with("p.adj_")),
    names_to = "metric_type",
    values_to = "value"
  ) %>%
  unite("metric_outcome", metric_type, outcome, sep = "_") %>%
  pivot_wider(
    names_from = metric_outcome,
    values_from = value
  ) 

res_health_NAcor = res_health_wide

## Lung
f6 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Lung)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_diff_Lung = first(cor_diff_Lung),
                        cor_diff_Lung = first(cor_diff_Lung)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Lung), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between AAR (diff) and Lung Flag") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))


f7 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Lung)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_lm_Lung = first(cor_lm_Lung),
                        cor_lm_Lung = first(cor_lm_Lung)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Lung), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between lm residual and Lung Flag") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f6
f7

## VESOP
f10 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_VESOP)) + 
  geom_tile(colour = "black") +  
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_diff_VESOP = first(cor_diff_VESOP),
                        cor_diff_VESOP = first(cor_diff_VESOP)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_VESOP), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between diff and VESOP") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f11 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_VESOP)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_lm_VESOP = first(cor_lm_VESOP),
                        cor_lm_VESOP = first(cor_lm_VESOP)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_VESOP), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between lm residual and VESOP") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f10
f11

## Cortisol
f12 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Cortisol)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_diff_Cortisol = first(cor_diff_Cortisol),
                        cor_diff_Cortisol = first(cor_diff_Cortisol)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Cortisol), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between diff and Cortisol") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f13 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Cortisol)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_lm_Cortisol = first(cor_lm_Cortisol),
                        cor_lm_Cortisol = first(cor_lm_Cortisol)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Cortisol), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between lm residual and Cortisol") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f12
f13

## Anemia
f14 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Anemia)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_diff_Anemia = first(cor_diff_Anemia),
                        cor_diff_Anemia = first(cor_diff_Anemia)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Anemia), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between diff and Anemia") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f15 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Anemia)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_lm_Anemia = first(cor_lm_Anemia),
                        cor_lm_Anemia = first(cor_lm_Anemia)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Anemia), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between lm residual and Anemia") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f14
f15

## Neutrophils
f16 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Neutrophils)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_diff_Neutrophils = first(cor_diff_Neutrophils),
                        cor_diff_Neutrophils = first(cor_diff_Neutrophils)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_diff_Neutrophils), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between diff and Neutrophils") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f17 = res_health_NAcor |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Neutrophils)) + 
  geom_tile(colour = "black") +
  geom_tile(data = res_health_NAcor %>%
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        cor_lm_Neutrophils = first(cor_lm_Neutrophils),
                        cor_lm_Neutrophils = first(cor_lm_Neutrophils)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = cor_lm_Neutrophils), width = Inf, height = 1, color = "black") +
  scale_x_discrete(expand = expansion(add = 0)) +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "Correlation between lm residual and Neutrophils") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

f16
f17

f10t = f10 + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), legend.position = "none") + 
  labs(subtitle = "(a) Correlation between difference-based AAR and VESOP")
f6t = f6 + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), legend.position = "none") + 
  labs(subtitle = "(b) Correlation between difference-based AAR and lung disease")
f12t = f12 + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") + 
  labs(subtitle = "(c) Correlation between difference-based AAR and low cortisol")
f11t = f11 + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none") + 
  labs(subtitle = "(d) Correlation between model-based AAR and VESOP")
f7t = f7 + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none") + 
  labs(subtitle = "(e) Correlation between model-based AAR and lung disease")
f13t = f13 + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none") + 
  labs(subtitle = "(f) Correlation between model-based AAR and low cortisol")


fhealthall = (f10t + f6t + f12t + f11t + f7t + f13t) &  
  theme(legend.position = "bottom", legend.box.margin = margin(-20, 0, 0, 0),
                                axis.ticks.x = element_blank(), axis.ticks.y = element_blank()) 
fhealthall = fhealthall + plot_layout(ncol = 2, guides = "collect", byrow = FALSE) 
fhealthall

ggsave("output/figs/fhealth.png", fhealthall, width = 10.5, height = 9, dpi = 300)

# plot of all health outcomes

res_health_long1 = res_health %>% pivot_longer(
  cols = cor_diff_pregnant:p_lm_Glucose,
  names_to = c("stat", "resid_type", "outcome"),
  names_pattern = "(.*)_(.*)_(.*)",
  values_to = "value") |> 
  filter(outcome %in% c("VESOP", "Lung", "Anemia", "Cortisol", "Neutrophils", "Globulin", "Glucose")) 

res_health_long2 = res_health_long1 |>
  filter(stat == "p") |> 
  mutate(stat = "p.adj") |> 
  mutate(value = p.adjust(value, method = "fdr")) 

res_health_long1 = rbind.data.frame(res_health_long1, res_health_long2)
rm(res_health_long2)

res_health_wide = res_health_long1 |> 
  pivot_wider(
    names_from = c(stat, resid_type),
    names_glue = "{stat}_{resid_type}",
    values_from = value
  ) 

res_health_wide = res_health_wide |> 
  mutate(cor_diff = ifelse(p.adj_diff > 0.05, NA, cor_diff)) |> 
  mutate(cor_lm = ifelse(p.adj_lm > 0.05, NA, cor_lm)) 

res_health_NAcor = res_health_wide %>%
  pivot_longer(
    cols = c(starts_with("cor_"), starts_with("df_"), starts_with("p_"), starts_with("p.adj_")),
    names_to = "metric_type",
    values_to = "value"
  ) 


t1 = res_health_NAcor  |> 
  filter(lbl == "Age\nCtl\nWtd\n") |> 
  filter(metric_type == "cor_lm") |> 
  ggplot(aes(x = outcome, y = forcats::fct_rev(metric), fill = value)) + 
  geom_tile(colour = "black") +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(b) With untransformed age") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-20, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

t1

t2 = res_health_NAcor  |> 
  filter(lbl == "LL\nCtl\nWtd\n") |> 
  filter(metric_type == "cor_lm") |> 
  ggplot(aes(x = outcome, y = forcats::fct_rev(metric), fill = value)) + 
  geom_tile(colour = "black") +
  scico::scale_fill_scico(palette = "vik", limits = c(-0.5,0.5), na.value="white") +
  labs(x = element_blank(), y = element_blank(),
       subtitle = "(a) Correlations between log-linear age acceleration and health") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-10, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(colour = ycols))

t2

t3 = t1 + theme(axis.text.y = element_blank()) 

t4 = t2 + t3 & theme(legend.position = "bottom", legend.box.margin = margin(-10, 0, 0, 0),
                     axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
t4 = t4 + plot_layout(ncol = 2, guides = "collect", byrow = TRUE) 
t4

ggsave("output/figs/all_health_cors.png", t4, width = 12, height = 4, dpi = 300)

# 
# 
# 
# ########
# res_health_long = res_health %>% 
#   select(metric, lbl, cor_diff_pregnant:p_lm_AlkPhos) |> 
#   pivot_longer(
#     cols = cor_diff_pregnant:p_lm_AlkPhos,
#     names_to = c("stat", "resid_type", "outcome"),
#     names_pattern = "(.*)_(.*)_(.*)",
#     values_to = "value") |> 
#   filter(stat == "cor", outcome %in% c("VESOP", "Lung", "Anemia", "Cortisol", "Neutrophils")) |> 
#   pivot_wider(names_from = resid_type, values_from = value) |> 
#   mutate(lm_better_cor = abs(lm) > abs(diff),
#          lm_same_cor = abs(lm) == abs(diff)) |> 
#   filter(!str_detect(metric, "Probage")) |> 
#   filter(!str_detect(metric, "Pheno")) |> 
#   filter(!(metric == "BND" & lbl != "Age\nAll\nUn\n"))|> 
#   filter(!(metric == "Cetaceans" & lbl != "Age\nAll\nUn\n"))|> 
#   filter(!(metric == "Mammals" & lbl != "Age\nAll\nUn\n"))
# 
# table(res_health_long$metric)  
# 
# res_health_long |> summarize(proplm = sum(lm_better_cor) / n(), propsame = sum(lm_same_cor) / n())
# res_health_long |> group_by(metric) |> summarize(proplm = sum(lm_better_cor) / n(), propsame = sum(lm_same_cor) / n())
# res_health_long |> group_by(lbl) |> summarize(proplm = sum(lm_better_cor) / n(), propsame = sum(lm_same_cor) / n())
# res_health_long |> group_by(outcome) |> summarize(proplm = sum(lm_better_cor) / n(), propsame = sum(lm_same_cor) / n())
# 
# res_health_long |> filter(!str_detect(metric, "Pheno")) |> summarize(proplm = sum(lm_better_cor) / n(), propsame = sum(lm_same_cor) / n())
