# data consists of:
## x: betas
## y: age
## weights: vector of sample weights ~ uncertainty in actual age
## control: vector, 1 = healthy population, 0 = other

# Load packages
library(glmnet)
library(dplyr)
library(tidyr)
library(stringr)
library(here)
library(readxl)
library(lubridate)
library(survival)
library(flexsurv)
library(randomForestSRC)

set.seed(47)
nfolds = 10 

betas = readRDS("data/Barratclough2024_betas476.rds")
betas = betas |> mutate(year = str_sub(Date, 1, 4))
# check if any missing age
which(is.na(betas$age)) 

# going to match health and betas so read in health here
health <- read_xlsx(path = here("data/wild_BND_health_filtered_sheet.xlsx"), sheet = "health_sheet")

# # drop any animals here
# betas <-  betas |> filter(Site != "USN")

# DOLPHIN.ID as combination of Site + year + Field.ID 
betas = betas |> mutate(DOLPHIN.ID = paste0(Site, "_", year, "_", Field.ID),
                        Site_Field.ID = paste0(Site, "_", Field.ID))

health = health |> mutate(year = year(Date)) |> 
  mutate(DOLPHIN.ID = paste0(Site, "_", year, "_", Field.ID),
         Site_Field.ID = paste0(Site, "_", Field.ID))

table(table(betas$Site_Field.ID)) # number of samples per individual (Site_Field.ID)
table(table(health$Site_Field.ID)) # number of samples per individual (Site_Field.ID)

# add vesop and flags to health data
load("data/VESOP_w_flags.Rdata")
flags = flags |> mutate(year = str_sub(Date, 1, 4)) |> 
  mutate(DOLPHIN.ID = paste0(Site, "_", year, "_", Field.ID)) |> 
  select(DOLPHIN.ID, VESOP_surv:Lung_Flag)
health = health |> left_join(flags, by = join_by(DOLPHIN.ID))

# bb = betas |> add_count(Site_Field.ID) |> select(DOLPHIN.ID, Site_Field.ID, year, Site, Date, Age, n)
# hh = health |> add_count(Site_Field.ID) |> select(DOLPHIN.ID, Site_Field.ID, year, Site, Date, Age.at.sampling, n) |> mutate(inhh = 1)
# 
# bb0 = bb |> filter(n > 1)
# hh0 = hh |> filter(n > 1) 
# 
# bb = left_join(bb, hh |> select(DOLPHIN.ID, Date, inhh), by = join_by(DOLPHIN.ID))

# # for animals with >1 sample, take the first
# betas = betas |> group_by(Site_Field.ID) |> slice_min(year, n = 1) |> ungroup()

betas = betas |> arrange(DOLPHIN.ID)

#############################################
### DNAm betas
#############################################

betas.ind <- str_detect(colnames(betas), "^cg")
colnames(betas)[!betas.ind]
x = as.matrix(betas[, betas.ind])
y = as.numeric(betas$Age)
y[y == 0] <- 0

#############################################
### Sample metadata
#############################################

sample_info = betas |> 
  select(DOLPHIN.ID, Site_Field.ID, Field.ID, year, site = Site, age = Age, weight = Sample.weight, sex = Sex) |> 
  mutate(control = ifelse(site %in% c("USN", "CHS", "SAR", "IRL"), 1, 0)) |> 
  mutate(year = as.numeric(year)) |> 
  mutate(age = pmax(0, age))

# add fold id

# ensure all individual's observations go into the same fold
fold_alloc = sample_info |> group_by(Site_Field.ID) |> 
  summarize(n = n(), control = first(control)) |> ungroup() 

fold_alloc = fold_alloc |> 
  group_by(control) |> 
  mutate(kf.foldid = rep(1:nfolds,times=ceiling(n()/nfolds),length.out=n())) |> 
  ungroup() 

# Create a list of indicators per fold
for (i in 1:nfolds) {
  fold_alloc <- fold_alloc |>
    group_by(control) |>
    mutate(
      is_test = (kf.foldid == i)
    ) |>
    group_modify(~ {
      n_rows <- nrow(.x)
      is_test_vec <- .x$is_test
      
      # Number of non-test items in this group
      n_non_test <- sum(!is_test_vec)
      
      # Create val allocation vector for non-test items
      val_inds <- sample(c(rep(TRUE, round(0.15 * n_non_test)), 
                           rep(FALSE, n_non_test - round(0.15 * n_non_test))))
      
      # Full-length vector: insert TRUE only in non-test positions
      temp_val_vec <- replace(rep(FALSE, n_rows), which(!is_test_vec), val_inds)
      
      # Return with new train/val columns for this fold
      .x |>
        mutate(
          !!paste0("kf.train", i) := (!is_test_vec & !temp_val_vec),
          !!paste0("kf.val", i) := (!is_test_vec & temp_val_vec)
        )
    }) |>
    ungroup() 
}

fold_alloc <- fold_alloc |>
  relocate(starts_with("kf.train"), .after = last_col()) |>
  relocate(starts_with("kf.val"), .after = last_col()) |> 
  rowwise() |> 
  mutate(ntrain = sum(c_across(starts_with("kf.train"))),
         nval = sum(c_across(starts_with("kf.val")))) |> 
  ungroup()


# add another fold id so that the N(dolphins in kf.redfoldid i) = N(control dolphins in kf.foldid i),
# drawn from both control and non-control individiduals
tt <- sum(fold_alloc$control) / nfolds
fold_alloc <- fold_alloc |> 
  slice_sample(n = nrow(fold_alloc)) |> 
  group_by(kf.foldid) |> 
  mutate(kf.redfoldid = ifelse(row_number() <= tt, kf.foldid, 0)) |> 
  ungroup() 

# Create a list of indicators per fold
for (i in 1:nfolds) {
  fold_alloc <- fold_alloc |>
    group_by(control) |>
    mutate(
      is_test = (kf.redfoldid == i)
    ) |>
    group_modify(~ {
      n_rows <- nrow(.x)
      is_test_vec <- .x$is_test
      
      # Number of non-test items in this group
      n_non_test <- sum(!is_test_vec)
      
      # Create val allocation vector for non-test items
      val_inds <- sample(c(rep(TRUE, round(0.15 * n_non_test)), 
                           rep(FALSE, n_non_test - round(0.15 * n_non_test))))
      
      # Full-length vector: insert TRUE only in non-test positions
      temp_val_vec <- replace(rep(FALSE, n_rows), which(!is_test_vec), val_inds)
      
      # Return with new train/val columns for this fold
      .x |>
        mutate(
          !!paste0("kf.redtrain", i) := (!is_test_vec & !temp_val_vec),
          !!paste0("kf.redval", i) := (!is_test_vec & temp_val_vec)
        )
    }) |>
    ungroup()
}

fold_alloc <- fold_alloc |>
  relocate(starts_with("kf.redtrain"), .after = last_col()) |>
  relocate(starts_with("kf.redval"), .after = last_col())

# fold_alloc = fold_alloc |> select(Site_Field.ID, kf.foldid, kf.redfoldid)
fold_alloc = fold_alloc |> select(Site_Field.ID, starts_with("kf."), starts_with("kf.red"))

sample_info <- sample_info |> left_join(fold_alloc, by = join_by(Site_Field.ID))

sample_info = sample_info |> arrange(DOLPHIN.ID)

table(sample_info$control, sample_info$kf.foldid)
table(sample_info$control, sample_info$kf.redfoldid)

table(sample_info$site, sample_info$kf.foldid)
table(sample_info$site, sample_info$kf.redfoldid)

#############################################
### Health data matched to DNAm
#############################################

# add indicator variable = 1 if have DNAm data for this DOLPHIN.ID (id_site_year combo)
dnam.ids = unique(paste0(sample_info$DOLPHIN.ID))
health = health |> mutate(betas.ind = ifelse(DOLPHIN.ID %in% dnam.ids, 1, 0)) |> 
  mutate(control = ifelse(Site %in% c("USN", "CHS", "SAR", "IRL"), 1, 0))

# TRUE if have DNAm for all but Navy animals
sum(health$betas.ind) == nrow(betas |> filter(Site != "USN"))

# add Navy animals to health dataframe (will be NAs for all columns except DOLPHIN.ID)
health_all = health
health = data.frame(DOLPHIN.ID = dnam.ids) |> left_join(health, by = join_by(DOLPHIN.ID))

health <- health |> select(DOLPHIN.ID, Site_Field.ID, site = Site, sex, 
                           first_sighting = Initial.sighting, last_sighting = Status_date, 
                           status = Last_known_status, sampling_date = Date, age_at_sampling = Age.at.sampling, 
                           everything())

# check order of rows is same in x/betas, health, and sample_info df's
ordbetas <- betas$DOLPHIN.ID
ordhealth = health$DOLPHIN.ID
ordsinfo = sample_info$DOLPHIN.ID
which(ordbetas != ordsinfo)  # should be empty
which(ordhealth != ordsinfo) # should be empty

save(x, y, sample_info, health, health_all, 
     file = "output/dolphin_data.Rdata")


# not used from here on - PhenoAge type analysis

# #############################################
# ### Health data for survival model
# #############################################
# 
# # for animals with >1 health measurement two ways to do things:
# # (1) take one observation, entry_time = that sampling time, exit_time = last_sighting, status = status at last_sighting
# #     - in this case probably makes sense to take the first sample, but sort of arbitrary
# # (2) use all observations, with entry_time = sampling time and exit_time either next sampling time or last sighting
# #.    depending on some rules (see below). Need to ensure correct status is used.
# # I put both of these into one data frame so one can filter out later on.
# 
# health_all <- health_all |> select(DOLPHIN.ID, Site_Field.ID, site = Site, sex, 
#                                    first_sighting = Initial.sighting, last_sighting = Status_date, 
#                                    status = Last_known_status, sampling_date = Date, age_at_sampling = Age.at.sampling, 
#                                    everything())
# 
# # adding relevant time variables
# health_all <- health_all |> group_by(Site_Field.ID) |> 
#   arrange(sampling_date) |> 
#   mutate(
#     age_at_sampling = as.numeric(age_at_sampling),
#     gt1_sample = n() > 1,
#     first_sample = sampling_date == min(sampling_date),
#     last_sample = sampling_date == max(sampling_date),
#     next_sampling_date = lead(sampling_date)) |> 
#   ungroup() |> 
#   mutate(birth_date = sampling_date - 365 * 24 * 60 * 60 * age_at_sampling) 
# 
# study_start = min(health_all$sampling_date, na.rm = TRUE)
# 
# # NB: condition-dependent entry and exit times
# # a) one sample and last_sighting == sampling_date: drop
# # b) >1 sample, last_sample == TRUE, last_sighting == sampling_date: drop
# # c) one sample and last_sighting > sampling_date: entry_time = sampling date, exit_time = last_sighting, status = status
# # d) >1 sample, last_sample == FALSE, next_sampling_date > sampling_date, next_sampling_date < last_sighting: entry_time = sampling_date, exit_time = next_sampling_date, status = "Alive"
# # e) >1 sample, last_sample == TRUE, last_sighting > sampling_date: entry_time = sampling_date, exit_time = last_sighting, status = status.
# 
# # conditions a, b
# health_all = health_all |> 
#   filter(gt1_sample | (last_sighting != sampling_date)) |> 
#   filter(!last_sample | (last_sighting != sampling_date))
# 
# # conditions c-e
# # times / 100 else get convergence problems with flexsurvreg gompertz model
# health_all = health_all |> 
#   mutate(
#     entry_time = case_when(
#       !gt1_sample & (last_sighting > sampling_date) ~ as.numeric(sampling_date - study_start, "days") / 100,
#       gt1_sample & !last_sample & (next_sampling_date > sampling_date) & (last_sighting > next_sampling_date) ~ as.numeric(sampling_date - study_start, "days") / 100,
#       gt1_sample & last_sample & (last_sighting > sampling_date) ~ as.numeric(sampling_date - study_start, "days") / 100,
#       TRUE ~ NA
#     ),
#     exit_time = case_when(
#       !gt1_sample & (last_sighting > sampling_date) ~ as.numeric(last_sighting - study_start, "days") / 100,
#       gt1_sample & !last_sample & (next_sampling_date > sampling_date) & (last_sighting > next_sampling_date) ~ as.numeric(next_sampling_date - study_start, "days") / 100,
#       gt1_sample & last_sample & (last_sighting > sampling_date) ~ as.numeric(last_sighting - study_start, "days") / 100,
#       TRUE ~ NA
#     ),
#     exit_status = case_when(
#       !gt1_sample & (last_sighting > sampling_date) ~ status,
#       gt1_sample & !last_sample & (next_sampling_date > sampling_date) & (last_sighting > next_sampling_date) ~ "Alive",
#       gt1_sample & last_sample & (last_sighting > sampling_date) ~ status,
#       TRUE ~ NA
#     )
#   )
# 
# # if you're going to use one observation per animal then the entry and exit times should be, respectively, the 
# # earliest entry_time and latest exit_time for that animal
# health_all = health_all |> 
#   group_by(Site_Field.ID) |>
#   mutate(entry_time_1obs = min(entry_time),
#          exit_time_1obs = max(exit_time),
#          exit_status_1obs = last(status, order_by = exit_time)) |> 
#   ungroup() |> 
#   mutate(time_from_entry_1obs = exit_time_1obs - entry_time_1obs)
# 
# # add indicator TRUE if this observation is the one to use if one obs per animal 
# health_all = health_all |> group_by(Site_Field.ID) |>
#   mutate(ind_1obs = entry_time_1obs == min(entry_time)) |> 
#   ungroup()
# 
# # add indicator for dead
# health = health |> mutate(dead = ifelse(status == "Dead", 1, 0)) 
# health_all = health_all |> 
#   mutate(dead = ifelse(exit_status == "Dead", 1, 0),
#          dead_1obs = ifelse(exit_status_1obs == "Dead", 1, 0)) 
# 
# health = health |> rename(age = age_at_sampling)
# health_all = health_all |> rename(age = age_at_sampling)
# 
# # previously (a) removed NA age's, (b) cases where entry_time > exit_time,
# # (c) pre-screened health indicators (variables), (b) scaled health indicators
# # now left to do later in analysis
# 
# # deciding which indicators to include
# 
# # https://www.randomforestsrc.org/articles/survival.html#illustration
# health_tmp <- health_all |> select(time_from_entry_1obs, dead_1obs, A.G.Ratio:WBC)
# rf0 <- rfsrc(Surv(time = time_from_entry_1obs, event = dead_1obs) ~., data = health_tmp, ntree = 1000, nodesize = 5, nsplit = 50, importance = TRUE)
# # health_tmp <- health_all |> select(entry_time, exit_time, dead, A.G.Ratio:WBC)
# # rf0 <- rfsrc(Surv(time = entry_time, time2 = exit_time, event = dead) ~., data = health_tmp, ntree = 1000, nodesize = 5, nsplit = 50, importance = TRUE)
# jk.obj <- subsample(rf0)
# sort(jk.obj$rf$importance, decreasing = TRUE)[1:15]
# # plot(jk.obj, xlab = "Variable Importance (x 100)")
# vars_to_keep = names(sort(jk.obj$rf$importance, decreasing = TRUE)[1:10])
# # number of NAs per good variable
# prop_NA = apply(health |> select(all_of(vars_to_keep)), 2, function(x) sum(is.na(x)) / length(x)) 
# prop_NA
# # drop any variables with more than 25% missing
# vars_to_keep = vars_to_keep[prop_NA < 0.25]
# vars_to_keep
# 
# health_survmod <- health_all |> filter(!is.na(age)) |> filter(exit_time > entry_time)
# # health_all <- health_all |> filter_at(vars_to_keep, all_vars(!is.na(.)))
# health_survmod <- health_survmod |> select(DOLPHIN.ID, Site_Field.ID, control, entry_time, exit_time, dead, age, all_of(vars_to_keep)) 
# health_survmod <- health_survmod[complete.cases(health_survmod), ]
# 
# # if obs appears in sample_info then it gets the same fold ids, and randomly allocate the rest
# health_survmod1 = health_survmod |> left_join(fold_alloc, by = join_by(Site_Field.ID)) 
# health_survmod2 = health_survmod1 |> filter(is.na(kf.foldid)) |> select(-kf.foldid, -kf.redfoldid)
# 
# # these have folds
# health_survmod1 = health_survmod1 |> filter(!is.na(kf.foldid))
# 
# # to add folds to those obs that need them, use same fold allocation procedure as above
# fold_alloc = health_survmod2 |> group_by(Site_Field.ID) |> 
#   summarize(n = n(), control = first(control)) |> ungroup() 
# 
# fold_alloc = fold_alloc |> 
#   group_by(control) |> 
#   mutate(kf.foldid = rep(1:nfolds,times=ceiling(n()/nfolds),length.out=n())) |> 
#   ungroup() 
# 
# tt <- sum(fold_alloc$control) / nfolds
# fold_alloc <- fold_alloc |> 
#   slice_sample(n = nrow(fold_alloc)) |> 
#   group_by(kf.foldid) |> 
#   mutate(kf.redfoldid = ifelse(row_number() <= tt, kf.foldid, 0)) |> 
#   ungroup() 
# 
# fold_alloc = fold_alloc |> select(Site_Field.ID, kf.foldid, kf.redfoldid)
# 
# health_survmod2 <- health_survmod2 |> left_join(fold_alloc, by = join_by(Site_Field.ID))
# 
# health_survmod <- rbind(health_survmod1, health_survmod2) |> arrange(DOLPHIN.ID)
# rm(health_survmod1, health_survmod2)
# 
# sample_info = sample_info |> mutate(have_health = DOLPHIN.ID %in% health_survmod$DOLPHIN.ID)
# 
# # check order of rows is same in x/betas, health, and sample_info df's
# ordbetas <- betas$DOLPHIN.ID
# ordhealth = health$DOLPHIN.ID
# ordsinfo = sample_info$DOLPHIN.ID
# which(ordbetas != ordsinfo)  # should be empty
# which(ordhealth != ordsinfo) # should be empty
# 
# save(x, y, sample_info, health, health_all, health_survmod, vars_to_keep, 
#      file = "output/dolphin_data.Rdata")
# 
