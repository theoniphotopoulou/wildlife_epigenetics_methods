# # remove comment if want to run in standalone mode rather than through run_all.R and go to line 240
# library(glmnet)
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(randomForest)
# source("age-transformations-zoller.R")
# set.seed(47)
# load("output/dolphin_data.Rdata")

# function to run elastic net
run_en_nocatch = function(x, sample_info, en_alpha = 0.5, n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE, loglin_age = "lin", loglin_sexm = NULL){
  # downsample matches sample size to what is returned if controls_only == TRUE. Only used if controls_only == FALSE
  
  nfolds = max(sample_info$kf.foldid)
  
  test.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  train.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  val.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  
  for(i in 1:nfolds){
    
    if (loglin_age == "loglin") { 
      y = fun_llin3.trans(sample_info$age, loglin_sexm)
    } else if (loglin_age == "log"){
      y = log(pmax(sample_info$age, 0.01))
    } else if (loglin_age == "lin") {
      y = sample_info$age
    } else { stop("Unknown age transformation, must be one of log, loglin, lin")}
    
    # inclusion indicator for training model (all except fold i, taking filters in account)
    train_col <- paste0("kf.train", i)
    train_include <- sample_info[[train_col]] 
    redtrain_col <- paste0("kf.redtrain", i)
    redtrain_include <- sample_info[[redtrain_col]] 
    if (controls_only) {
      inc <- train_include & sample_info$control == 1
    } else if (reduce_size) {
      inc <- redtrain_include & sample_info$kf.redfoldid > 0
    } else {
      inc <- train_include
    }
    
    xt <- x[inc, ]
    yt <- y[inc]
    
    if(!is.null(n_cpg)){
      cors = cor(xt, yt)
      inds = order(abs(cors), decreasing = TRUE)[1:n_cpg]
      xt = xt[, inds] 
      x = x[, inds]
    }
    
    # inclusion indicator for validation
    val_col <- paste0("kf.val", i)
    val_include <- sample_info[[val_col]] 
    redval_col <- paste0("kf.redval", i)
    redval_include <- sample_info[[redval_col]] 
    if (controls_only) {
      inc_val <- val_include & sample_info$control == 1
    } else if (reduce_size) {
      inc_val <- redval_include & sample_info$kf.redfoldid > 0
    } else {
      inc_val <- val_include
    }
    # inclusion indicator for testing (all fold i observations)
    inc_pred <- sample_info$kf.foldid == i
    
    if(use_weights){ 
      wts <- sample_info$weights[inc]
      fit <- cv.glmnet(xt, yt, foldid = 1:nrow(xt), alpha = en_alpha, weights = wts)
    } else {
      fit <- cv.glmnet(xt, yt, foldid = 1:nrow(xt), alpha = en_alpha)
    } 
    
    pred <- as.numeric(predict(fit, s = fit$lambda.min, newx = x))
    
    if (loglin_age == "loglin") { 
      pred = fun_llin3.inv(pred, loglin_sexm)
    } else if (loglin_age == "log"){
      pred = exp(pred)
    } 
    
    # predict
    test.pred <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc_pred], preds = pred[inc_pred])
    train.pred <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc], preds = pred[inc])
    val.pred <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc_val], preds = pred[inc_val])
    
    test.pred.df <- left_join(test.pred.df, test.pred, by = join_by(DOLPHIN.ID))
    train.pred.df <- left_join(train.pred.df, train.pred, by = join_by(DOLPHIN.ID))
    val.pred.df <- left_join(val.pred.df, val.pred, by = join_by(DOLPHIN.ID))
    
  }
  
  #tmp = ifelse(is.null(n_cpg), ncol(x), n_cpg)
  test.pred.df$metric = paste0("fgec_en")
  train.pred.df$metric = paste0("fgec_en")
  val.pred.df$metric = paste0("fgec_en")
  colnames(test.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(train.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(val.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  
  return(list(testpreds = test.pred.df, trainpreds = train.pred.df, valpreds = val.pred.df))
  
}

# function to run random forest
run_rf_nocatch = function(x, sample_info, n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE, 
                          loglin_age = "lin", loglin_sexm = NULL, rf_mtry = 1/3){
  # downsample matches sample size to what is returned if controls_only == TRUE. Only used if controls_only == FALSE
  
  nfolds = max(sample_info$kf.foldid)
  
  test.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  train.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  val.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  
  for(i in 1:nfolds){
    
    if (loglin_age == "loglin") { 
      y = fun_llin3.trans(sample_info$age, loglin_sexm)
    } else if (loglin_age == "log"){
      y = log(pmax(sample_info$age, 0.01))
    } else if (loglin_age == "lin") {
      y = sample_info$age
    } else { stop("Unknown age transformation, must be one of log, loglin, lin")}
    
    # inclusion indicator for training model (all except fold i, taking filters in account)
    train_col <- paste0("kf.train", i)
    train_include <- sample_info[[train_col]] 
    redtrain_col <- paste0("kf.redtrain", i)
    redtrain_include <- sample_info[[redtrain_col]] 
    if (controls_only) {
      inc <- train_include & sample_info$control == 1
    } else if (reduce_size) {
      inc <- redtrain_include & sample_info$kf.redfoldid > 0
    } else {
      inc <- train_include
    }
    
    xt <- x[inc, ]
    yt <- y[inc]
    
    if(!is.null(n_cpg)){
      cors = cor(xt, yt)
      inds = order(abs(cors), decreasing = TRUE)[1:n_cpg]
      xt = xt[, inds] 
      x = x[, inds]
    }
    
    # inclusion indicator for validation
    val_col <- paste0("kf.val", i)
    val_include <- sample_info[[val_col]] 
    redval_col <- paste0("kf.redval", i)
    redval_include <- sample_info[[redval_col]] 
    if (controls_only) {
      inc_val <- val_include & sample_info$control == 1
    } else if (reduce_size) {
      inc_val <- redval_include & sample_info$kf.redfoldid > 0
    } else {
      inc_val <- val_include
    }
    # inclusion indicator for testing (all fold i observations)
    inc_pred <- sample_info$kf.foldid == i
    
    if(use_weights){ 
      wts <- sample_info$weights[inc]
      fit <- randomForest(xt, yt, mtry = max(500, floor(ncol(xt) * rf_mtry)), importance = FALSE, ntree = 1000, weights = wts)
    } else {
      fit <- randomForest(xt, yt, mtry = max(500, floor(ncol(xt) * rf_mtry)), importance = FALSE, ntree = 1000)
    } 
    
    pred = as.numeric(predict(fit, newdata = x))
    
    if (loglin_age == "loglin") { 
      pred = fun_llin3.inv(pred, loglin_sexm)
    } else if (loglin_age == "log"){
      pred = exp(pred)
    } 
    
    # predict
    test.pred <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc_pred], preds = pred[inc_pred])
    train.pred <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc], preds = pred[inc])
    val.pred <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc_val], preds = pred[inc_val])
    
    test.pred.df <- left_join(test.pred.df, test.pred, by = join_by(DOLPHIN.ID))
    train.pred.df <- left_join(train.pred.df, train.pred, by = join_by(DOLPHIN.ID))
    val.pred.df <- left_join(val.pred.df, val.pred, by = join_by(DOLPHIN.ID))
    
  }
  
  #tmp = ifelse(is.null(n_cpg), ncol(x), n_cpg)
  test.pred.df$metric = paste0("fgec_rf")
  train.pred.df$metric = paste0("fgec_rf")
  val.pred.df$metric = paste0("fgec_rf")
  colnames(test.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(train.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(val.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  
  return(list(testpreds = test.pred.df, trainpreds = train.pred.df, valpreds = val.pred.df))
  
}

run_en = function(x, sample_info, en_alpha = 0.5, n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE, loglin_age = "lin", loglin_sexm = NULL){
  
  tmp <- try(
    run_en_nocatch(x = x, sample_info = sample_info, en_alpha = en_alpha, n_cpg = n_cpg, use_weights = use_weights, controls_only = controls_only, 
                   reduce_size = reduce_size, loglin_age = loglin_age, loglin_sexm = loglin_sexm), 
    silent = TRUE)
  if(class(tmp) != "try-error"){
    tmp = tmp
  } else {
    metric_name = paste0("fgec_en")
    empty_df <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
    empty_df[paste0("preds", 1:nfolds)] <- NA
    tmp = list(testpreds = empty_df |> mutate(metric = metric_name),
               trainpreds = empty_df |> mutate(metric = metric_name))
  }
  tmp
}

run_rf = function(x, sample_info, rf_mtry = 1/3, n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE, loglin_age = "lin", loglin_sexm = NULL){
  
  tmp <- try(
    run_rf_nocatch(x = x, sample_info = sample_info, rf_mtry = rf_mtry, n_cpg = n_cpg, use_weights = use_weights, controls_only = controls_only, 
                   reduce_size = reduce_size, loglin_age = loglin_age, loglin_sexm = loglin_sexm), 
    silent = TRUE)
  if(class(tmp) != "try-error"){
    tmp = tmp
  } else {
    metric_name = paste0("fgec_rf")
    empty_df <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
    empty_df[paste0("preds", 1:nfolds)] <- NA
    tmp = list(testpreds = empty_df |> mutate(metric = metric_name),
               trainpreds = empty_df |> mutate(metric = metric_name))
  }
  tmp
}

# # run once to check
# m0 <- run_en_nocatch(x = x, n_cpg = NULL, sample_info = sample_info, en_alpha = 0.5, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)
# m1 <- run_rf_nocatch(x = x, n_cpg = 500, sample_info = sample_info, rf_mtry = 1/3, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)

# m0 <- run_en(x = x, n_cpg = NULL, sample_info = sample_info, en_alpha = 0.5, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)
# m1 <- run_rf(x = x, n_cpg = 500, sample_info = sample_info, rf_mtry = 1/3, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)

# # leaving checks on alpha and age transform in here; also sensitive to fold allocations 
# m0 <- run_en(x = x, n_cpg = NULL, sample_info = sample_info, en_alpha = 0.12, use_weights = TRUE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "log", loglin_sexm = 15)
# m1 <- run_en(x = x, n_cpg = NULL, sample_info = sample_info, en_alpha = 0.5, use_weights = TRUE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "log", loglin_sexm = 15)
# m2 <- run_en(x = x, n_cpg = NULL, sample_info = sample_info, en_alpha = 0.12, use_weights = TRUE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)
# m3 <- run_en(x = x, n_cpg = NULL, sample_info = sample_info, en_alpha = 0.5, use_weights = TRUE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)
# 
# test0 = m0$testpreds |>
#   filter(if_any(starts_with("preds"), ~ !is.na(.x))) |>
#   mutate(predage = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
#   mutate(age = sample_info$age,
#          control = sample_info$control) |> 
#   mutate(age_class = case_when(
#     age <= 2 ~ "Calf",
#     age <= 15 ~ "Subadult",
#     age <= 25 ~ "Adult",
#     age > 25 ~ "Older")) |> 
#   mutate(age_class = factor(age_class, levels = c("Calf", "Subadult", "Adult", "Older")))
# 
# test1 = m1$testpreds |>
#   filter(if_any(starts_with("preds"), ~ !is.na(.x))) |>
#   mutate(predage = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
#   mutate(age = sample_info$age,
#          control = sample_info$control) |> 
#   mutate(age_class = case_when(
#     age <= 2 ~ "Calf",
#     age <= 15 ~ "Subadult",
#     age <= 25 ~ "Adult",
#     age > 25 ~ "Older")) |> 
#   mutate(age_class = factor(age_class, levels = c("Calf", "Subadult", "Adult", "Older")))
# 
# test2 = m2$testpreds |>
#   filter(if_any(starts_with("preds"), ~ !is.na(.x))) |>
#   mutate(predage = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
#   mutate(age = sample_info$age,
#          control = sample_info$control) |> 
#   mutate(age_class = case_when(
#     age <= 2 ~ "Calf",
#     age <= 15 ~ "Subadult",
#     age <= 25 ~ "Adult",
#     age > 25 ~ "Older")) |> 
#   mutate(age_class = factor(age_class, levels = c("Calf", "Subadult", "Adult", "Older")))
# 
# test3 = m3$testpreds |>
#   filter(if_any(starts_with("preds"), ~ !is.na(.x))) |>
#   mutate(predage = rowMeans(pick(starts_with("preds")), na.rm = TRUE)) |> 
#   mutate(age = sample_info$age,
#          control = sample_info$control) |> 
#   mutate(age_class = case_when(
#     age <= 2 ~ "Calf",
#     age <= 15 ~ "Subadult",
#     age <= 25 ~ "Adult",
#     age > 25 ~ "Older")) |> 
#   mutate(age_class = factor(age_class, levels = c("Calf", "Subadult", "Adult", "Older")))
# 
# t1 = test0 |> group_by(age_class) |> summarize(mae = median(abs(predage - age))) |> 
#   mutate(alpha = 0.12, age_transform = "log", run = 1)
# t2 = test1 |> group_by(age_class) |> summarize(mae = median(abs(predage - age))) |> 
#   mutate(alpha = 0.5, age_transform = "log", run = 1)
# t3 = test2 |> group_by(age_class) |> summarize(mae = median(abs(predage - age))) |> 
#   mutate(alpha = 0.12, age_transform = "loglin", run = 1)
# t4 = test3 |> group_by(age_class) |> summarize(mae = median(abs(predage - age))) |> 
#   mutate(alpha = 0.5, age_transform = "loglin", run = 1)
# 
# rbind(t1, t2, t3, t4)