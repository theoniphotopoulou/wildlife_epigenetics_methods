# # remove comment if want to run in standalone mode rather than through run_all.R and go to line 265

# library(glmnet)
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(randomForest)
# source("age-transformations-zoller.R")
# set.seed(47)
# load("output/dolphin_data.Rdata")
# health$VESOP_Flag = ifelse(health$VESOP_surv < 0.925, 1, 0)

# function to run elastic net
run_en_health_nocatch = function(x, yvar, y_transform = "lin", health, sample_info, en_alpha = 0.5, en_family = "gaussian", n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE){
  # downsample matches sample size to what is returned if controls_only == TRUE. Only used if controls_only == FALSE
  
  nfolds = max(sample_info$kf.foldid)
  
  invlogit <- function(x){exp(x)/(exp(x) + 1)}
  logit <- function(x){log(x/(1-x))}
  
  test.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  train.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  val.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  
  for(i in 1:nfolds){
    
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
    
    y = health[, yvar]
    
    if (y_transform == "logit") { 
      y = logit(y)
    } else if (y_transform == "log"){
      y = log(pmax(y, 0.001))
    } else if (y_transform == "lin") {
      y = y
    } else { stop("Unknown transformation, must be one of log, logit, lin")}
    
    inc = inc & !is.na(y)
    
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
      fit <- cv.glmnet(xt, yt, foldid = 1:nrow(xt), alpha = en_alpha, weights = wts, family = en_family)
    } else {
      fit <- cv.glmnet(xt, yt, foldid = 1:nrow(xt), alpha = en_alpha, family = en_family)
    } 
    
    pred <- as.numeric(predict(fit, s = fit$lambda.min, newx = x, type = "response"))
    
    if (y_transform == "logit") { 
      pred = invlogit(pred)
    } else if (y_transform == "log"){
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
  
  test.pred.df$metric = paste0("fgec_en")
  train.pred.df$metric = paste0("fgec_en")
  val.pred.df$metric = paste0("fgec_en")
  colnames(test.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(train.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(val.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  
  return(list(testpreds = test.pred.df, trainpreds = train.pred.df, valpreds = val.pred.df))
  
}

# function to run random forest
run_rf_health_nocatch = function(x, yvar, y_transform = "lin", health, sample_info, rf_mtry = 1/3, n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE){
  # downsample matches sample size to what is returned if controls_only == TRUE. Only used if controls_only == FALSE
  
  nfolds = max(sample_info$kf.foldid)
  
  test.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  train.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  val.pred.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  
  for(i in 1:nfolds){
    
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
    
    y = health[, yvar]
    
    if (y_transform == "logit") { 
      y = logit(y)
    } else if (y_transform == "log"){
      y = log(pmax(y, 0.001))
    } else if (y_transform == "factor"){
      y_unt = y
      y = factor(y)
    } else if (y_transform == "lin") {
      y = y
    } else { stop("Unknown transformation, must be one of log, logit, factor, lin")}
    
    inc = inc & !is.na(y)
    
    xt <- x[inc, ]
    yt <- y[inc]
    
    if(!is.null(n_cpg)){
      if(y_transform == "factor"){
        cors = cor(xt, y_unt[inc])
      } else {
        cors = cor(xt, yt)
      }
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
    
    pred = as.numeric(predict(fit, newdata = x, type = "prob")[, 2])
    
    if (y_transform == "logit") { 
      pred = invlogit(pred)
    } else if (y_transform == "log"){
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
  
  test.pred.df$metric = paste0("fgec_rf")
  train.pred.df$metric = paste0("fgec_rf")
  val.pred.df$metric = paste0("fgec_rf")
  colnames(test.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(train.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(val.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  
  return(list(testpreds = test.pred.df, trainpreds = train.pred.df, valpreds = val.pred.df))
  
}

run_en_health = function(x, yvar, y_transform = "lin", health, sample_info, en_alpha = 0.5, en_family = "gaussian", n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE){
  
  tmp <- try(
    run_en_health_nocatch(x = x, yvar = yvar, y_transform = y_transform, health = health, sample_info = sample_info, en_alpha = en_alpha, en_family = en_family, n_cpg = n_cpg, use_weights = use_weights, controls_only = controls_only, 
                          reduce_size = reduce_size), 
    silent = TRUE)
  if(class(tmp) != "try-error"){
    tmp = tmp
  } else {
    metric_name = paste0("fgec_en")
    empty_df <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
    empty_df[paste0("preds", 1:nfolds)] <- NA
    tmp = list(testpreds = empty_df |> mutate(metric = metric_name),
               trainpreds = empty_df |> mutate(metric = metric_name),
               valpreds = empty_df |> mutate(metric = metric_name))
  }
  tmp
}

run_rf_health = function(x, yvar, y_transform = "lin", health, sample_info, rf_mtry = 1/3, n_cpg = NULL, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE){
  
  tmp <- try(
    run_rf_health_nocatch(x = x, yvar = yvar, y_transform = y_transform, 
                          health = health, sample_info = sample_info, 
                          rf_mtry = rf_mtry, n_cpg = n_cpg, 
                          use_weights = use_weights, controls_only = controls_only, reduce_size = reduce_size), 
    silent = TRUE)
  if(class(tmp) != "try-error"){
    tmp = tmp
  } else {
    metric_name = paste0("fgec_rf")
    empty_df <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
    empty_df[paste0("preds", 1:nfolds)] <- NA
    tmp = list(testpreds = empty_df |> mutate(metric = metric_name),
               trainpreds = empty_df |> mutate(metric = metric_name),
               valpreds = empty_df |> mutate(metric = metric_name))
  }
  tmp
}

# m0 <- run_en_health(x = x, yvar = "VESOP_Flag", y_transform = "lin", health = health, sample_info = sample_info, en_alpha = 0.5, n_cpg = 1000, en_family = "binomial", use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE)
# m1 <- run_rf_health(x = x, yvar = "VESOP_Flag", y_transform = "factor", health = health, sample_info = sample_info, rf_mtry = 1/3, n_cpg = 500, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE)
# 
# m0 <- run_en_nocatch(x = x, n_cpg = 500, sample_info = sample_info, en_alpha = 0.5, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)
# m1 <- run_rf_nocatch(x = x, n_cpg = 500, sample_info = sample_info, rf_mtry = 1/3, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)

# m0 <- run_en(x = x, n_cpg = 500, sample_info = sample_info, en_alpha = 0.5, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)
# m1 <- run_rf(x = x, n_cpg = 500, sample_info = sample_info, rf_mtry = 1/3, use_weights = FALSE, controls_only = TRUE, reduce_size = FALSE, loglin_age = "loglin", loglin_sexm = 15)


# test = m1$testpreds |>
#   filter(if_any(preds1:preds5, ~ !is.na(.x))) |>
#   mutate(predage = rowMeans(pick(preds1:preds5), na.rm = TRUE))
# test$age = sample_info$age
# test$control = sample_info$control
# test = test |> mutate(age_class = case_when(
#   age <= 2 ~ "Calf",
#   age <= 15 ~ "Subadult",
#   age <= 25 ~ "Adult",
#   age > 25 ~ "Older"
# ))
# test |> group_by(age_class) |> summarize(cor = cor(age, predage), mae = median(abs(predage - age)))
# test |> group_by(control) |> summarize(cor = cor(age, predage), mae = median(abs(predage - age)))
# 
# train = m0$trainpreds |> 
#   filter(if_any(preds1:preds5, ~ !is.na(.x))) |> 
#   mutate(predage = rowMeans(pick(preds1:preds5), na.rm = TRUE)) 
# train$age = unlist(sample_info[sample_info$control == 1, "age"])
# plot(train$age, train$predage)
# cor(train$age, train$predage)
# 
# df = data.frame(trueage = test$age, testpredage = test$predage, trainpredage = train$predage, control = sample_info$control)
# cor(df)
# pairs(df)
# 
# 
# # 
# # 
# # # elastic net
# # 
# # ## 
# # 
# # en.fit <- cv.glmnet(x, y, foldid = my.foldid, alpha = 0.5, weights= sample_info$weights)
# # en.fit.optalpha <- cv.glmnet(x, y, foldid = my.foldid, alpha = 0.01, weights= sample_info$weights) # determined elsewhere
# # 
# # # random forest
# # rf.fit <- randomForest(x, y, mtry = 4)
# # rf.opt <- tuneRF(x,y, mtryStart = 512, ntreeTry = 1000, stepFactor = 2)
# # #determines the optimal number of variables to sample at each split as 6 (approx sqrt(40))
# # rf.fit <- randomForest(x, y, mtry = 2048, importance = TRUE, ntree = 1000)