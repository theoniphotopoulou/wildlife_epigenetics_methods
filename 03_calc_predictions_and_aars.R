library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(mgcv)
# library(scico)
# library(scales)
# library(patchwork)
# library(glue)
# library(kableExtra)
source("age-transformations-zoller.R")

load("output/model_output.Rdata")
load("output/dolphin_data.Rdata")

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

alltest_en_05 = do.call(rbind, purrr::map(res_en_05, "alltest")) |> mutate(metric = "fgec_en05")
alltest_en_opt = do.call(rbind, purrr::map(res_en_opt, "alltest")) |> mutate(metric = "fgec_enopt")
alltest_rf_05 = do.call(rbind, purrr::map(res_rf_05, "alltest")) |> mutate(metric = "fgec_rf05")
alltest_rf_opt = do.call(rbind, purrr::map(res_rf_opt, "alltest")) |> mutate(metric = "fgec_rfopt")
alltest_probage = do.call(rbind, purrr::map(res_pr, "alltest"))
alltest = alltest_en_05 |> 
  rbind.data.frame(alltest_en_opt) |> 
  rbind.data.frame(alltest_rf_05) |> 
  rbind.data.frame(alltest_rf_opt) |> 
  rbind.data.frame(alltest_probage)

alltrain_en_05 = do.call(rbind, purrr::map(res_en_05, "alltrain")) |> mutate(metric = "fgec_en05")
alltrain_en_opt = do.call(rbind, purrr::map(res_en_opt, "alltrain")) |> mutate(metric = "fgec_enopt")
alltrain_rf_05 = do.call(rbind, purrr::map(res_rf_05, "alltrain")) |> mutate(metric = "fgec_rf05")
alltrain_rf_opt = do.call(rbind, purrr::map(res_rf_opt, "alltrain")) |> mutate(metric = "fgec_rfopt")
alltrain_probage = do.call(rbind, purrr::map(res_pr, "alltrain"))
alltrain = alltrain_en_05 |> 
  rbind.data.frame(alltrain_en_opt) |> 
  rbind.data.frame(alltrain_rf_05) |> 
  rbind.data.frame(alltrain_rf_opt) |> 
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

# AARs are residuals from linear model 
allpreds = allpreds |> 
  group_by(metric, par_set) |> 
  mutate(aar_lm_test = case_when(
    str_detect(metric, "probage") ~ testpredage, 
    TRUE ~ as.numeric(predict(lm(age ~ testpredage, data = cur_data(), na.action = na.exclude)) - age))) |> 
  ungroup()

# AARs are residuals from gam 
allpreds = allpreds |> 
  group_by(metric, par_set) |> 
  mutate(aar_gam_test = case_when(
    str_detect(metric, "probage") ~ testpredage, 
    TRUE ~ as.numeric(predict(gam(age ~ s(testpredage, k = 4), data = cur_data(), na.action = na.exclude)) - age))) |> 
  ungroup()

allpreds = allpreds |> select(-starts_with("aar[0-9]"))
allpreds = allpreds |> select(-starts_with("testpreds[0-9]")) 

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
  ungroup() 

colnames(existing_clocks)
colnames(allpreds)

allpreds = bind_rows(allpreds, existing_clocks) 

# ensembles

table(allpreds$metric)

# ensemble names for unoptimised version
ens_names = c("fgec_enopt", "fgec_rfopt", "Robeck2021", "Unpub20xx", "Universal") 
ensembles = allpreds |> 
  filter(metric %in% ens_names) |> 
  select(DOLPHIN.ID, age, metric, lbl, testpredage) |> 
  pivot_wider(
    names_from = metric,
    values_from = c(testpredage)
  ) |> 
  group_by(lbl) |>
  rowwise() |> 
  mutate(ens_all = mean(c_across(fgec_enopt:Universal), na.rm = TRUE)) |> 
  mutate(ens_existing = mean(c_across(Robeck2021:Universal), na.rm = TRUE)) |> 
  ungroup() |> 
  select(DOLPHIN.ID, age, lbl, ens_existing, ens_all)

# estimate best coefficients
ens_names = c("fgec_enopt", "fgec_rfopt", "Robeck2021", "Unpub20xx", "Universal") 
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
  ungroup() 

colnames(ensembles)
colnames(allpreds)

allpreds = bind_rows(allpreds, ensembles)         

#

allpreds = allpreds |> 
  mutate(rb_diff_test = aar_diff_test / age,
         rb_lm_test = aar_lm_test / age,
         rb_gam_test = aar_gam_test / age)

allpreds = allpreds |> 
  mutate(age_class = case_when(
    age <= 2 ~ "Calf",
    age <= 15 ~ "Subadult",
    age <= 25 ~ "Adult",
    age > 25 ~ "Older"
  )) |> 
  left_join(health, by = join_by(DOLPHIN.ID))

save(allpreds, file = "output/predicted_ages_and_aars.Rdata")
