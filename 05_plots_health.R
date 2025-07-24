library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(pROC)
library(tidyr)
library(glue)
library(patchwork)

source("get_performance_stats.R")

load("output/dolphin_data.Rdata")
health$VESOP_Flag = ifelse(health$VESOP_surv < 0.925, 1, 0)

###############
# Elastic net #
###############

load("output/health_en_model_output.Rdata")

yy = c()
ps = c()
for(i in 1:length(pars_en$yvar)){
  yy = c(yy, health[, pars_en$yvar[i]])
}
df = data.frame(obs = yy)

alltest = do.call(rbind, purrr::map(res_en, "alltest")) |> cbind(df)
alltrain = do.call(rbind, purrr::map(res_en, "alltrain")) |> cbind(df)
allval = do.call(rbind, purrr::map(res_en, "allval")) |> cbind(df)

allval = allval |> 
  mutate(predy = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
  mutate(predy = ifelse(is.nan(predy), NA, predy)) |> 
  left_join(pars_en, by = join_by(par_set))

alltest = alltest |> 
  mutate(predy = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
  mutate(predy = ifelse(is.nan(predy), NA, predy)) |> 
  left_join(pars_en, by = join_by(par_set))

# identify best combination of (alpha, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means highest AUC on the validation set
val_aucs <- allval |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      auc = as.numeric(auc(roc_obj))
    )
  })
val_aucs = val_aucs |> left_join(pars_en, by = join_by(par_set))
val_aucs_best = val_aucs |> slice_max(auc, by = yvar)
best_parsets_en = unlist(val_aucs_best$par_set)

# identify best combination of (alpha = 0.5, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means highest AUC on the validation set
val_aucs <- allval |> 
  filter(en_alpha == 0.5) |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      auc = as.numeric(auc(roc_obj))
    )
  })
val_aucs = val_aucs |> left_join(pars_en, by = join_by(par_set))
val_aucs_best = val_aucs |> slice_max(auc, by = yvar)
best_parsets_en_n_cpg = unlist(val_aucs_best$par_set)

# identify best combination of (alpha, n_cpg = all) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means highest AUC on the validation set
val_aucs <- allval |> 
  filter(n_cpg == ncol(x)) |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      auc = as.numeric(auc(roc_obj))
    )
  })
val_aucs = val_aucs |> left_join(pars_en, by = join_by(par_set))
val_aucs_best = val_aucs |> slice_max(auc, by = yvar)
best_parsets_en_alpha = unlist(val_aucs_best$par_set)

# Test set ROC curves for all par_sets 
roc_data_all <- alltest |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      specificity = rev(roc_obj$specificities),
      sensitivity = rev(roc_obj$sensitivities),
      auc = as.numeric(auc(roc_obj))
    )
  }) |> 
  mutate(model = "Elastic net") |> 
  left_join(pars_en, by = "par_set") |> 
  mutate(yvar_auc = glue("{yvar} (AUC = {round(auc, 3)})")) 

# Performance statistics for all elastic net 
f1_points_en = get_performance_stats(valdata = allval, testdata = alltest, opt_threshold = TRUE) |> 
  left_join(pars_en, by = join_by(par_set))

# Test set ROC curves for unoptimized elastic net 
roc_data_en_05 <- roc_data_all |> 
  filter(en_alpha == 0.5, n_cpg == ncol(x)) |> 
  mutate(optim = "Unoptimized")

f1_en_05 = f1_points_en |> 
  filter(en_alpha == 0.5, n_cpg == ncol(x)) |> 
  mutate(optim = "Unoptimized")

# Test set ROC curves for optimized elastic net (both alpha and n_cpg)
roc_data_en_opt <- roc_data_all |> 
  filter(par_set %in% best_parsets_en) |> 
  mutate(optim = "Optimized")

f1_en_opt <- f1_points_en |> 
  filter(par_set %in% best_parsets_en) |> 
  mutate(optim = "Optimized")

# Test set ROC curves for optimized elastic net (n_cpg only)
roc_data_en_opt_n_cpg <- roc_data_all |> 
  filter(par_set %in% best_parsets_en_n_cpg) |> 
  mutate(optim = "Optimized nCpG")

f1_en_opt_n_cpg <- f1_points_en |> 
  filter(par_set %in% best_parsets_en_n_cpg)  |> 
  mutate(optim = "Optimized nCpG")

# Test set ROC curves for optimized elastic net (alpha only)
roc_data_en_opt_alpha <- roc_data_all |> 
  filter(par_set %in% best_parsets_en_alpha)  |> 
  mutate(optim = "Optimized alpha")

f1_en_opt_alpha <- f1_points_en |> 
  filter(par_set %in% best_parsets_en_alpha) |> 
  mutate(optim = "Optimized alpha")

roc_data_en = rbind.data.frame(roc_data_en_05, roc_data_en_opt, roc_data_en_opt_alpha, roc_data_en_opt_n_cpg)
f1_en = rbind.data.frame(f1_en_05, f1_en_opt, f1_en_opt_alpha, f1_en_opt_n_cpg)

roc_data_en = roc_data_en |> mutate(model = "Elastic net")
f1_en = f1_en |> mutate(model = "Elastic net")

#################
# Random forest #
#################

load("output/health_rf_model_output.Rdata")

yy = c()
ps = c()
for(i in 1:length(pars_rf$yvar)){
  yy = c(yy, health[, pars_rf$yvar[i]])
}
df = data.frame(obs = yy)

alltest = do.call(rbind, purrr::map(res_rf, "alltest")) |> cbind(df)
alltrain = do.call(rbind, purrr::map(res_rf, "alltrain")) |> cbind(df)
allval = do.call(rbind, purrr::map(res_rf, "allval")) |> cbind(df)

allval = allval |> 
  mutate(predy = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
  mutate(predy = ifelse(is.nan(predy), NA, predy)) |> 
  left_join(pars_rf, by = join_by(par_set))

alltest = alltest |> 
  mutate(predy = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
  mutate(predy = ifelse(is.nan(predy), NA, predy)) |> 
  left_join(pars_rf, by = join_by(par_set))

# identify best combination of (mtry, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means highest AUC on the validation set
val_aucs <- allval |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      auc = as.numeric(auc(roc_obj))
    )
  })
val_aucs = val_aucs |> left_join(pars_rf, by = join_by(par_set))
val_aucs_best = val_aucs |> slice_max(auc, by = yvar)
best_parsets_rf = unlist(val_aucs_best$par_set)

# identify best combination of (mtry = 1/3, n_cpg) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means highest AUC on the validation set
val_aucs <- allval |> 
  filter(rf_mtry == 1/3) |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      auc = as.numeric(auc(roc_obj))
    )
  })
val_aucs = val_aucs |> left_join(pars_rf, by = join_by(par_set))
val_aucs_best = val_aucs |> slice_max(auc, by = yvar)
best_parsets_rf_n_cpg = unlist(val_aucs_best$par_set)

# identify best combination of (mtry, n_cpg = all) for each modelling condition (c(use_weights, controls_only, reduce_size, loglin_age))
# "best" means highest AUC on the validation set
val_aucs <- allval |> 
  filter(n_cpg == ncol(x)) |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      auc = as.numeric(auc(roc_obj))
    )
  })
val_aucs = val_aucs |> left_join(pars_rf, by = join_by(par_set))
val_aucs_best = val_aucs |> slice_max(auc, by = yvar)
best_parsets_rf_mtry = unlist(val_aucs_best$par_set)

# Test set ROC curves for all par_sets 
roc_data_all <- alltest |> 
  group_split(par_set) %>%
  map_dfr(~{
    par_set <- .x$par_set[1]
    roc_obj <- roc(response = .x$obs, predictor = .x$predy, quiet = TRUE, smooth = FALSE)
    tibble(
      par_set = par_set,
      specificity = rev(roc_obj$specificities),
      sensitivity = rev(roc_obj$sensitivities),
      auc = as.numeric(auc(roc_obj))
    )
  }) |> 
  mutate(model = "Elastic net") |> 
  left_join(pars_rf, by = "par_set") |> 
  mutate(yvar_auc = glue("{yvar} (AUC = {round(auc, 3)})")) 

# Performance statistics for all rf 
f1_points_rf = get_performance_stats(valdata = allval, testdata = alltest, opt_threshold = TRUE) |> 
  left_join(pars_rf, by = join_by(par_set))

# Performance statistics for all rf *optimizing on the test set* - to demo wrong way
f1_points_rf_wrong = get_performance_stats(valdata = alltest, testdata = alltest, opt_threshold = TRUE) |> 
  left_join(pars_rf, by = join_by(par_set))

# Test set ROC curves for unoptimized elastic net 
roc_data_rf_05 <- roc_data_all |> 
  filter(rf_mtry == 1/3, n_cpg == ncol(x)) |> 
  mutate(optim = "Unoptimized")

f1_rf_05 = f1_points_rf |> 
  filter(rf_mtry == 1/3, n_cpg == ncol(x)) |> 
  mutate(optim = "Unoptimized")

# Test set ROC curves for optimized elastic net (both mtry and n_cpg)
roc_data_rf_opt <- roc_data_all |> 
  filter(par_set %in% best_parsets_rf) |> 
  mutate(optim = "Optimized")

f1_rf_opt <- f1_points_rf |> 
  filter(par_set %in% best_parsets_rf) |> 
  mutate(optim = "Optimized")

# Test set ROC curves for optimized elastic net (n_cpg only)
roc_data_rf_opt_n_cpg <- roc_data_all |> 
  filter(par_set %in% best_parsets_rf_n_cpg) |> 
  mutate(optim = "Optimized nCpG")

f1_rf_opt_n_cpg <- f1_points_rf |> 
  filter(par_set %in% best_parsets_rf_n_cpg)  |> 
  mutate(optim = "Optimized nCpG")

# Test set ROC curves for optimized elastic net (mtry only)
roc_data_rf_opt_mtry <- roc_data_all |> 
  filter(par_set %in% best_parsets_rf_mtry)  |> 
  mutate(optim = "Optimized mtry")

f1_rf_opt_mtry <- f1_points_rf |> 
  filter(par_set %in% best_parsets_rf_mtry) |> 
  mutate(optim = "Optimized mtry")

roc_data_rf = rbind.data.frame(roc_data_rf_05, roc_data_rf_opt, roc_data_rf_opt_mtry, roc_data_rf_opt_n_cpg)
f1_rf = rbind.data.frame(f1_rf_05, f1_rf_opt, f1_rf_opt_mtry, f1_rf_opt_n_cpg)

roc_data_rf = roc_data_rf |> mutate(model = "Random forest")
f1_rf = f1_rf |> mutate(model = "Random forest")

##########
# Tables #
##########

outcomes_to_keep = unique(roc_data_en$yvar)
roc_vars = intersect(names(roc_data_en), names(roc_data_rf))
f1_vars = intersect(names(f1_en), names(f1_rf))

roc_data = rbind.data.frame(roc_data_en |> select(all_of(roc_vars)), 
                            roc_data_rf |> select(all_of(roc_vars))) |> 
  filter(yvar %in% outcomes_to_keep) |> 
  mutate(yvar = str_extract(yvar, "^[^_]+")) 

f1_data = rbind.data.frame(f1_en |> select(all_of(f1_vars)), 
                           f1_rf |> select(all_of(f1_vars))) |> 
  filter(yvar %in% outcomes_to_keep) |> 
  mutate(yvar = str_extract(yvar, "^[^_]+")) |> 
  mutate(yvar = str_replace(yvar, "Neutrophils", "N'phils."))

xwide = f1_data |> filter(optim %in% c("Optimized", "Unoptimized")) |>  mutate(across(auc_val:recall, ~ . * 100))
x <- xwide %>% pivot_longer(cols = auc_val:recall, names_to = "metric")
x <- x %>% mutate(pt_size = ifelse(metric == "Accuracy", 2, 1))
x <- x %>% mutate(metric = factor(metric, levels = c("threshold", "acc", "sensitivity", 
                                                     "specificity", "precision", "recall", "f1", "auc_val", "auc_test"), 
                                  labels = c("Threshold", "Accuracy", "Sensitivity", 
                                             "Specificity", "Precision", "Recall", "F1", "AUC (val)", "AUC")))

base_acc = f1_data |> mutate(base_acc = 100 - 100 * (TP + FN) / (TP + TN + FP + FN)) |> 
  filter(optim %in% c("Optimized", "Unoptimized")) 

t6 <- x %>% filter(metric %in% c("Precision", "Recall", "F1")) %>% 
  ggplot(aes(x = yvar, y = value, col = optim, shape = metric, size = pt_size, group = optim)) +
  geom_point(position = position_dodge(0.5)) +
  geom_bar(data = xwide, width = 0.4, aes(x = yvar, y = acc, fill = optim, group = optim), inherit.aes = FALSE, alpha = 0.3, stat = "identity", position = position_dodge(0.5)) +
  geom_point(data = base_acc, aes(x = yvar, y = base_acc, colour = optim, group = optim), inherit.aes = FALSE, shape = "-", size = 9, position = position_dodge(0.5)) +
  geom_linerange(data = xwide, linewidth = 1, aes(x = yvar, col = optim, ymin = recall, ymax = precision, group = optim), inherit.aes = FALSE, position = position_dodge(0.5)) +
  facet_grid(. ~ model) +
  scale_radius(range = c(1.5,3)) +
  scale_colour_manual(values = c("#CC3300", "#99C9FF")) +
  guides(size = "none", shape = guide_legend(title = "Metric"), colour = guide_legend(title = "Tuning"), fill = guide_legend(title = "Tuning")) +
  xlab("Outcome") + ylab("Performance") + 
  ylim(0,100) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box.margin = margin(-10, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA)) 

t6

ggsave("output/figs/health_model_prf1.png", t6, width = 8, height = 3.5, dpi = 300)

# values for text
x |> filter(metric == "F1") |> group_by(optim, model) |> summarize(mv = mean(value))
x |> filter(metric == "F1", optim == "Optimized", model == "Elastic net") |> 
  arrange(desc(value)) |> select(yvar, value)
x |> filter(metric == "F1", optim == "Optimized", model == "Random forest") |> 
  arrange(desc(value)) |> select(yvar, value)
x |> filter(metric == "Precision", optim == "Optimized") |> 
  select(model, yvar, value)
x |> filter(metric == "F1") |> group_by(optim, model) |> summarize(mv = mean(value))

xt = x |> filter(metric == "F1") |> select(yvar, model, optim, metric, value) |> 
  pivot_wider(names_from = optim, values_from = value) |> 
  mutate(opt_imp = (Optimized - Unoptimized)) |> 
  mutate(opt_imp_rel = (Optimized - Unoptimized)/Unoptimized) |> 
  arrange(model, desc(opt_imp))
xt |> group_by(model) |> summarize(mdf1 = mean(opt_imp), mf1 = mean(Optimized)) 
xt 

x |> filter(metric == "AUC") |> group_by(optim, model) |> summarize(mv = mean(value))
x |> filter(metric == "AUC", optim == "Optimized", model == "Elastic net") |> 
  arrange(desc(value)) |> select(yvar, value)
x |> filter(metric == "AUC", optim == "Optimized", model == "Random forest") |> 
  arrange(desc(value)) |> select(yvar, value)

xt = x |> filter(metric == "AUC") |> select(yvar, model, optim, metric, value) |> 
  pivot_wider(names_from = optim, values_from = value) |> 
  mutate(opt_imp = (Optimized - Unoptimized)) |> 
  mutate(opt_imp_rel = (Optimized - Unoptimized)/Unoptimized) |> 
  arrange(model, desc(Optimized))
xt |> group_by(model) |> summarize(mdf1 = mean(opt_imp), mf1 = mean(Optimized)) 


#########
# Plots #
#########

outcomes_to_keep = c("VESOP", "Lung", "Cortisol", "Neutrophils")
roc_vars = intersect(names(roc_data_en), names(roc_data_rf))
f1_vars = intersect(names(f1_en), names(f1_rf))

roc_data = rbind.data.frame(roc_data_en |> select(all_of(roc_vars)), 
                            roc_data_rf |> select(all_of(roc_vars))) |> 
  mutate(yvar = str_extract(yvar, "^[^_]+")) |> 
  filter(optim == "Optimized") |> 
  filter(yvar %in% outcomes_to_keep)  
#  mutate(model = paste0(model, " opt"))

f1_data = rbind.data.frame(f1_en |> select(all_of(f1_vars)), 
                           f1_rf |> select(all_of(f1_vars))) |> 
  mutate(yvar = str_extract(yvar, "^[^_]+")) |> 
  filter(optim == "Optimized") |> 
  filter(yvar %in% outcomes_to_keep)  
 # mutate(model = paste0(model, " opt"))

xwide = f1_data |> filter(optim %in% c("Optimized", "Unoptimized")) |>  mutate(across(auc_val:recall, ~ . * 100))
x <- xwide %>% pivot_longer(cols = auc_val:recall, names_to = "metric")
x <- x %>% mutate(pt_size = ifelse(metric == "Accuracy", 2, 1))
x <- x %>% mutate(metric = factor(metric, levels = c("threshold", "acc", "sensitivity", 
                                                     "specificity", "precision", "recall", "f1", "auc_val", "auc_test"), 
                                  labels = c("Threshold", "Accuracy", "Sensitivity", 
                                             "Specificity", "Precision", "Recall", "F1", "AUC (val)", "AUC")))

p0_en_data = f1_en |> select(all_of(f1_vars)) |> 
  mutate(yvar = str_extract(yvar, "^[^_]+")) |> 
  filter(yvar %in% outcomes_to_keep) |> 
  filter(optim %in% c("Optimized", "Unoptimized")) |>  
  mutate(optim = factor(optim, levels = c("Optimized", "Unoptimized"))) |> 
  mutate(yvar = str_sub(yvar, 1, 1))

p0_rf_data = f1_rf |> select(all_of(f1_vars)) |> 
  mutate(yvar = str_extract(yvar, "^[^_]+")) |> 
  filter(yvar %in% outcomes_to_keep) |> 
  filter(optim %in% c("Optimized", "Unoptimized")) |>  
  mutate(optim = factor(optim, levels = c("Optimized", "Unoptimized"))) |> 
  mutate(yvar = str_sub(yvar, 1, 1))

p0_inset_EN <- p0_en_data |> 
  ggplot(aes(x = yvar, y = 100 * auc_test, fill = yvar)) +
  geom_bar(data = p0_en_data |> filter(optim == "Optimized"), stat = "identity", width = 0.4) +
  geom_point(data = p0_en_data |> filter(optim == "Unoptimized"), shape = "-", size = 6, color = "black") +
  annotate("text", x = 2.5, y = 95, label = "AUCs") +
  scale_fill_brewer(palette = "Set2") +
  guides(size = "none", fill = "none") +
  labs(x = element_blank(), y = element_blank()) +
  ylim(0, 100) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white"))

p0_inset_RF <- p0_rf_data |> 
  ggplot(aes(x = yvar, y = 100 * auc_test, fill = yvar)) +
  geom_bar(data = p0_rf_data |> filter(optim == "Optimized"), stat = "identity", width = 0.4) +
  geom_point(data = p0_rf_data |> filter(optim == "Unoptimized"), shape = "-", size = 6, color = "black") +
  annotate("text", x = 2.5, y = 95, label = "AUCs") +
  scale_fill_brewer(palette = "Set2") +
  guides(size = "none", fill = "none") +
  labs(x = element_blank(), y = element_blank()) +
  ylim(0, 100) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_rect(fill = "white"))

p0_inset_EN
p0_inset_RF
grob_inset_EN <- ggplotGrob(p0_inset_EN)
grob_inset_RF <- ggplotGrob(p0_inset_RF)

p0_EN = ggplot(roc_data |> filter(model == "Elastic net"), 
            aes(x = 1 - specificity, y = sensitivity, color = yvar, group = yvar)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, colour = "gray", linetype = 2) +
  geom_point(data = f1_data|> filter(model == "Elastic net"), aes(x = 1 - specificity, y = sensitivity, color = yvar), 
             size = 3, shape = 21, fill = "white", stroke = 1) +
  facet_grid(. ~ model) +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Parameter Set"
  ) +
  scale_colour_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-10, 0, 0, 0),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA)) 

p0_RF = ggplot(roc_data |> filter(model == "Random forest"), 
            aes(x = 1 - specificity, y = sensitivity, color = yvar, group = yvar)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, colour = "gray", linetype = 2) +
  geom_point(data = f1_data |> filter(model == "Random forest"), aes(x = 1 - specificity, y = sensitivity, color = yvar), 
             size = 3, shape = 21, fill = "white", stroke = 1) +
  facet_grid(. ~ model) +
  labs(
    x = "False Positive Rate (1 - Specificity)",
    y = element_blank(),
    color = "Parameter Set"
  ) +
  scale_colour_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-10, 0, 0, 0),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = NA)) 


p0_EN <- p0_EN +
  annotation_custom(grob = grob_inset_EN, 
                    xmin = 0.50, xmax = 1.05, ymin = -0.1, ymax = 0.50) 

p0_RF <- p0_RF +
  annotation_custom(grob = grob_inset_RF, 
                    xmin = 0.50, xmax = 1.05, ymin = -0.1, ymax = 0.50)

p0 = p0_EN + p0_RF &  
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.box.margin = margin(-10, 0, 0, 0),
        plot.margin = unit(c(0, 0, 0, 0), "pt")) 

p0 = p0 + plot_layout(ncol = 2, guides = "collect") 

p0

ggsave("output/figs/health_model_ROCs.png", p0, width = 8, height = 3.5, dpi = 300)

# ggplot(roc_data_rf, 
#        aes(x = 1 - specificity, y = sensitivity, color = yvar, group = yvar)) +
#   geom_line() +
#   geom_point(data = f1_rf, aes(x = 1 - specificity, y = sensitivity, color = yvar), 
#              size = 3, shape = 21, fill = "white", stroke = 1) +
#   facet_grid(. ~ optim) +
#   labs(
#     title = "ROC Curves by Parameter Set",
#     x = "False Positive Rate (1 - Specificity)",
#     y = "True Positive Rate (Sensitivity)",
#     color = "Parameter Set"
#   ) +
#   theme_minimal()
# 
# ggplot(roc_data_en, 
#        aes(x = 1 - specificity, y = sensitivity, color = optim)) +
#   geom_line() + 
#   geom_point(data = f1_en, aes(x = 1 - specificity, y = sensitivity, color = optim), 
#              size = 3, shape = 21, fill = "white", stroke = 1) +
#   facet_wrap(~yvar, ncol = 4) +
#   labs(
#     title = "ROC Curves by Parameter Set",
#     x = "False Positive Rate (1 - Specificity)",
#     y = "True Positive Rate (Sensitivity)",
#     color = "Parameter Set"
#   ) +
#   theme_minimal()
# 
# ggplot(roc_data_rf, 
#        aes(x = 1 - specificity, y = sensitivity, color = optim)) +
#   geom_line() + 
#   geom_point(data = f1_rf, aes(x = 1 - specificity, y = sensitivity, color = optim), 
#              size = 3, shape = 21, fill = "white", stroke = 1) +
#   facet_wrap(~yvar, ncol = 4) +
#   labs(
#     title = "ROC Curves by Parameter Set",
#     x = "False Positive Rate (1 - Specificity)",
#     y = "True Positive Rate (Sensitivity)",
#     color = "Parameter Set"
#   ) +
#   theme_minimal()
# 
# 
# 
# 
# 
