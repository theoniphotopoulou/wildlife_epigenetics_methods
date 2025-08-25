library(dplyr)
library(tidyr)

set.seed(47)

# define number of CPG sites
n_cpg <- 31139
n_samples <- 476

# create health data
synth_health <- data.frame(DOLPHIN.ID = c(1:n_samples),
                           site = c(rep("A",111), 
                                    rep("B",14),
                                    rep("C",65),
                                    rep("D",86),
                                    rep("E",11),
                                    rep("F",7), 
                                    rep("G",112),
                                    rep("H",20),
                                    rep("I",50)
                                    ),
                           sex = sample(c("Male","Female"), size=n_samples, replace=TRUE), 
                           age = round(runif(n=n_samples, min=0, max=50),digits=1),
                           VESOP_Flag = rbinom(n=n_samples, size=1, prob=0.925),
                           Lung_Flag = rbinom(n=n_samples, size=1, prob=0.2),
                           Anemia_Flag = rbinom(n=n_samples, size=1, prob=0.2),
                           Cortisol_lt30_Flag = rbinom(n=n_samples, size=1, prob=0.2),
                           Neutrophils_Flag = rbinom(n=n_samples, size=1, prob=0.2),
                           Globulin_Flag = rbinom(n=n_samples, size=1, prob=0.2),
                           Glucose_Flag = rbinom(n=n_samples, size=1, prob=0.2),
                           control = rbinom(n=n_samples, size=1, prob=0.5),
                           weight = runif(n=n_samples, min=0.27, max=1),
                           year = sample(x=(1999:2018), size=n_samples, replace=TRUE))
head(synth_health)


# Create beta values (proportion of methylation) 
random_digits <- vector("numeric")
random_string <- vector("character", length=n_cpg)

for (i in 1:n_cpg){
  random_digits <- sample(c(0:9), size = 8, replace = TRUE) # sample 8 random digits
  random_string[i] <- paste0("cg", paste(random_digits, collapse = "")) # stick them together in a string
}

synth_x <- matrix(data = rbeta(n_cpg*n_samples, 0.2, 0.2), 
             nrow = n_samples, ncol = n_cpg,
             dimnames = list(NULL,random_string))
dimnames(synth_x)
dim(synth_x)
mean(synth_x)

synth_dat <- cbind(synth_health, synth_x)
dim(synth_dat)

# create sample_info including fold info
synth_sample_info = synth_health |> 
  select(DOLPHIN.ID, year, site, age, weight, sex) |> 
  mutate(control = ifelse(site %in% c("A", "C", "D", "G"), 1, 0)) |> 
  mutate(year = as.numeric(year))

# add fold id

nfolds = 10 

fold_alloc = synth_sample_info |> 
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

fold_alloc = fold_alloc |> select(DOLPHIN.ID, starts_with("kf."), starts_with("kf.red"))

synth_sample_info <- synth_sample_info |> left_join(fold_alloc, by = join_by(DOLPHIN.ID))

synth_sample_info = synth_sample_info |> arrange(DOLPHIN.ID)

table(synth_sample_info$control, synth_sample_info$kf.foldid)
table(synth_sample_info$control, synth_sample_info$kf.redfoldid)

table(synth_sample_info$site, synth_sample_info$kf.foldid)
table(synth_sample_info$site, sample_info$kf.redfoldid)

x = synth_x
y = synth_sample_info$age
health = synth_health
sample_info = synth_sample_info

# optionally save the data
save(x, y, sample_info, health, 
     file = "output/synth_dolphin_data.Rdata")




