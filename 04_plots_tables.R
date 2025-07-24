library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scico)
library(scales)
library(patchwork)
library(glue)
library(kableExtra)
source("age-transformations-zoller.R")

load("output/dolphin_data.Rdata")
load("output/model_output.Rdata")
rm(list = setdiff(ls(), c("sample_info", "pars"))) # just keep these
load("output/predicted_ages_and_aars.Rdata")

# helper function for filtering and labelling clocks 
change_metric = function(res, ensemble = FALSE){
  
  if(ensemble){
    
    res = res |> 
      filter(metric %in% c("Unpub20xx", "Universal", "Robeck2021", 
                           "fgec_rf05", "fgec_en05", 
                           "fgec_rfopt", "fgec_enopt", 
                           "ens_all", "ens_lm",
                           "probage_mle_bias_1000", "probage_mle_acc_1000"))
    
    res = res |> 
      mutate(metric = factor(metric, 
                             levels = c("Robeck2021", "Unpub20xx", "Universal",
                                        "fgec_rf05", "fgec_en05", 
                                        "fgec_rfopt", "fgec_enopt", 
                                        "ens_all", "ens_lm",
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
                           "fgec_rf05", "fgec_en05", 
                           "fgec_rfopt", "fgec_enopt", 
                           "probage_mle_bias_1000", "probage_mle_acc_1000"))
    
    res = res |> 
      mutate(metric = factor(metric, 
                             levels = c("Robeck2021", "Unpub20xx", "Universal",
                                        "fgec_rf05", "fgec_en05", 
                                        "fgec_rfopt", "fgec_enopt", 
                                        "probage_mle_acc_1000", "probage_mle_bias_1000"),
                             labels = c("BND", "Cetaceans", "Mammals", 
                                        "Random forest", "Elastic net", 
                                        "Random forest opt", "Elastic net opt",
                                        "Probage Acc", "Probage Bias")))
  }
  
  res
}

################################################################################
## Table 1 : Median absolute error in years (lower quartile, upper quartile) 
## and relative bias (proportional error) of DNAm age..
################################################################################

tableres = allpreds |> 
  mutate(age_class = factor(age_class, levels = c("Calf", "Subadult", "Adult", "Older"))) |> 
  group_by(metric, par_set, lbl, age_class) |> 
  summarize(mae_test = list(round(as.numeric(quantile(abs(aar_diff_test), prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 1)),
            mrb_test = list(round(100 * as.numeric(quantile(rb_diff_test, prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 0)),
            n = n()) |> 
  mutate(mae_test = purrr::map_chr(mae_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})")),
         mrb_test = purrr::map_chr(mrb_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})"))) |> 
  ungroup()

tableres_all = allpreds |> 
  mutate(age_class = "All") |> 
  group_by(metric, par_set, lbl, age_class) |> 
  summarize(mae_test = list(round(as.numeric(quantile(abs(aar_diff_test), prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 1)),
            mrb_test = list(round(100 * as.numeric(quantile(rb_diff_test, prob = c(0.25, 0.5, 0.75), na.rm = TRUE)), 0)),
            n = n()) |> 
  mutate(mae_test = purrr::map_chr(mae_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})")),
         mrb_test = purrr::map_chr(mrb_test, ~ glue::glue("{.x[2]} ({.x[1]}; {.x[3]})"))) |> 
  ungroup()

tableres = rbind(tableres, tableres_all)

# use par_set == 10 to be as close as possible to Barratclough et al (2024)
# (loglinear age transformation used here rather than log)
tablepar = 10
tableres1 = tableres |> filter(par_set == tablepar) |> mutate(age_n = paste0(age_class, " (n = ", n, ")")) 
tableres1 = tableres1 |> pivot_wider(id_cols = metric, names_from = age_n, values_from = mae_test, names_sort = FALSE)
tableres1 = tableres1 |> change_metric(ensemble = TRUE) |> arrange(metric) #|> filter(!str_detect(metric, "Pheno"))

tableres2 = tableres |> filter(par_set == tablepar) |> mutate(age_n = paste0(age_class, " (n = ", n, ")")) 
tableres2 = tableres2 |> pivot_wider(id_cols = metric, names_from = age_n, values_from = mrb_test, names_sort = FALSE)
tableres2 = tableres2 |> change_metric(ensemble = TRUE) |> arrange(metric) #|> filter(!str_detect(metric, "Pheno"))

tableres12 = rbind(tableres1, tableres2)
tableres12
kable(tableres12, format = "latex", booktabs = TRUE)

# Preprocessing for figures

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
            mrb_gam_test = 100 * median(rb_gam_test, na.rm = TRUE)) |> 
  ungroup()

res = res |> change_metric(ensemble = TRUE)

# min-max correls per clock
res |> group_by(metric) |> summarize(minc = min(cor_testpredage, na.rm = TRUE),
                                     maxc = max(cor_testpredage, na.rm = TRUE))


# color clock names by the kind of clocks they are
## for figures where all clocks are shown
ycols = c(rep("purple", 2), rep("lightblue", 2), rep("blue", 4), rep("black", 3))
## for figures where pheno and probage clocks not shown
ycols23 = c(rep("lightblue", 2), rep("blue", 4), rep("black", 4))

# Predictions from existing clocks don't change with parameter settings so
# figures should merge tiles across the row for these clocks
row_to_merge <- c("BND", "Cetaceans", "Mammals") 

# Extract the unique fill value for that row
merged_row <- res %>%
  filter(metric %in% row_to_merge) %>%
  group_by(metric) |> 
  summarize(metric = unique(metric), lbl = first(lbl), 
            cor_testpredage = first(cor_testpredage),
            cor_testpredage = first(cor_testpredage)) 

################################################################################
## Figure 1: Associations between chronological and epigenetic age across 
## different modelling choices...
################################################################################

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

f0t = f0 + theme(axis.ticks.y = element_blank(), legend.position = "bottom") 
fcor = f0t + fscat + plot_layout(ncol = 2, widths = c(3,5))
fcor
ggsave("output/figs/fcorxy.png", fcor, width = 12, height = 5, dpi = 300)

################################################################################
## Figure 2: Accuracy of age predictions from epigenetic clocks (vertical axis) 
## across different modelling choices (horizontal axis), as measured by relative
## bias and MAE
################################################################################

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

f2lm = res |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res %>%
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

f2gam = res |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_gam_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res %>%
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

# MAE

f3 = res |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mae_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res %>%
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

f3lm = res |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res %>%
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

f3gam = res |> 
  filter(!str_detect(metric, "[Pp]robage") & !str_detect(metric, "[Pp]heno")) |> 
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = magame_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res %>%
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

f2t = f2 + theme(axis.text.x = element_blank()) 
f2lmt = f2lm + theme(axis.text.x = element_blank()) 
f2gamt = f2gam  

f3t = f3 + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 
f3lmt = f3lm + theme(axis.text.x = element_blank(), axis.text.y = element_blank()) 
f3gamt = f3gam + theme(axis.text.y = element_blank()) 

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
ggsave("output/figs/faartype.png", f23all, width = 10, height = 7, dpi = 300)

################################################################################
## Figure 3: Differences in the accuracy (relative bias, MAE) of age predictions 
## from epigenetic clocks (vertical axis) across different modelling choices 
## (horizontal axis) across different cohorts of the test sample.
################################################################################

# results within control and exposed groups
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

f2control_lm = res_with_control |> filter(control == 1) |> filter(!str_detect(metric, "[Pp]robage")) |>  
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control |> filter(control == 1) |> 
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_lm_test = first(mrb_lm_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test), width = Inf, height = 1, color = "black") +
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

f2exposed_lm = res_with_control |> filter(control == 0) |> filter(!str_detect(metric, "[Pp]robage")) |>  
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control |> filter(control == 1) |>  
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        mrb_lm_test = first(mrb_lm_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = mrb_lm_test), width = Inf, height = 1, color = "black") +
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

f3control_lm = res_with_control |> filter(control == 1) |> filter(!str_detect(metric, "[Pp]robage")) |>
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control |> filter(control == 0) |> 
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        malme_test = first(malme_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test), width = Inf, height = 1, color = "black") +
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

f3exposed_lm = res_with_control |> filter(control == 0) |> filter(!str_detect(metric, "[Pp]robage")) |>
  ggplot(aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test)) + 
  geom_tile(colour = "black") + 
  geom_tile(data = res_with_control |> filter(control == 0) |> 
              filter(metric %in% row_to_merge) %>%
              group_by(metric) |> 
              summarize(metric = unique(metric), lbl = first(lbl), 
                        malme_test = first(malme_test)), 
            aes(x = lbl, y = forcats::fct_rev(metric), fill = malme_test), width = Inf, height = 1, color = "black") +
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
ggsave("output/figs/faarbygroup.png", f45all, width = 10, height = 5, dpi = 300)

################################################################################
################################################################################

# HEALTH figures

## preprocessing

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

# min-max correlations for text
xt = res_health_long1 |>
  filter(stat == "cor") |> filter(!str_detect(lbl, "LL")) |>  
  group_by(metric, outcome, resid_type) |>
  summarize(minc = min(value, na.rm = TRUE),
            maxc = max(value, na.rm = TRUE)) |> ungroup()

# bigger correlations with lm than diff aars
xt = res_health_long1 |> 
  filter(stat == "cor") |> 
  filter(!str_detect(metric, "Probage")) |> 
  filter(!str_detect(metric, "Ensemble")) |> 
  group_by(metric, outcome, resid_type) |>
  summarize(meancor = mean((value), na.rm = TRUE)) |> ungroup() |> 
  pivot_wider(names_from = resid_type, values_from = meancor) |> 
  mutate(cordiff = (lm - diff)/ diff) |> 
  mutate(lmbiggercor = ifelse((abs(lm) > abs(diff)), 1, 0)) 

xt |> group_by(metric) |> summarize(meanlmbigger = mean(lmbiggercor))
mean(xt$lmbiggercor)

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

################################################################################
## Figure 4: Associations between two health outcomes and age acceleration 
## residuals (AARs) from epigenetic clocks...
################################################################################

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

################################################################################
## Figure 5: Associations between an expanded set of health outcomes (horizontal axis) and age
## acceleration from epigenetic clocks
################################################################################

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

