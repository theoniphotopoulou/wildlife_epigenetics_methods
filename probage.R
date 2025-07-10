# # remove comment if want to run in standalone mode rather than through run_all.R and go to line 104
# library(dplyr)
# source("age-transformations-zoller.R")
# source("probage_fns.R")
# set.seed(47)
# load("output/dolphin_data.Rdata")

# function to run probage
run_probage_nocatch = function(x, sample_info, n_cpg = 1000, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE, 
                       loglin_age = FALSE, loglin_sexm = NULL){
  # downsample matches sample size to what is returned if controls_only == TRUE. Only used if controls_only == FALSE
  
  nfolds = max(sample_info$kf.foldid)
  
  test.pred.acc.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  train.pred.acc.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  test.pred.bias.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  train.pred.bias.df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID)
  
  for(i in 1:nfolds){
    
    if (loglin_age == "loglin") { 
      y = fun_llin3.trans(sample_info$age, loglin_sexm)
      y = scale_to_range(y, max(0.01, min(sample_info$age, na.rm = TRUE)), max(sample_info$age, na.rm = TRUE))
    } else if (loglin_age == "log"){
      y = log(pmax(sample_info$age, 0.01))
      y = scale_to_range(y, max(0.01, min(sample_info$age, na.rm = TRUE)), max(sample_info$age, na.rm = TRUE))
    } else if (loglin_age == "lin") {
      y = sample_info$age
    } else { stop("Unknown age transformation, must be one of log, loglin, lin")}
    
    # inclusion indicator for training model (all except fold i, taking filters in account)
    if (controls_only) {
      inc <- sample_info$kf.foldid != i & sample_info$control == 1
    } else if (reduce_size) {
      inc <- sample_info$kf.foldid != i & sample_info$kf.redfoldid > 0
    } else {
      inc <- sample_info$kf.foldid != i
    }
    
    xt <- x[inc, ]
    yt <- y[inc]
    
    # inclusion indicator for testing (all fold i observations)
    inc_pred <- sample_info$kf.foldid == i

    if(use_weights){ 
      wts <- sample_info$weights[inc]
      fit <- probage_fit_sites(betas = xt, age = yt, n_cpg = n_cpg, weights = wts)
    } else {
      fit <- probage_fit_sites(betas = xt, age = yt, n_cpg = n_cpg)
    } 
    
    # predict
    pred <- probage_fit_animals(betas = x, age = y, site_results = fit)
    
    pred.acc <- pred$acc
    pred.bias <- pred$bias
    
    test.pred.acc <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc_pred], preds = pred.acc[inc_pred])
    train.pred.acc <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc], preds = pred.acc[inc])
    
    test.pred.bias <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc_pred], preds = pred.bias[inc_pred])
    train.pred.bias <- data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID[inc], preds = pred.bias[inc])
    
    test.pred.acc.df <- left_join(test.pred.acc.df, test.pred.acc, by = join_by(DOLPHIN.ID))
    train.pred.acc.df <- left_join(train.pred.acc.df, train.pred.acc, by = join_by(DOLPHIN.ID))
    
    test.pred.bias.df <- left_join(test.pred.bias.df, test.pred.bias, by = join_by(DOLPHIN.ID))
    train.pred.bias.df <- left_join(train.pred.bias.df, train.pred.bias, by = join_by(DOLPHIN.ID))
    
  }
  
  metric1_name = paste0("probage_mle_acc_", n_cpg)
  metric2_name = paste0("probage_mle_bias_", n_cpg)
  test.pred.df <- rbind(test.pred.acc.df |> mutate(metric = metric1_name), test.pred.bias.df |> mutate(metric = metric2_name))
  train.pred.df <- rbind(train.pred.acc.df |> mutate(metric = metric1_name), train.pred.bias.df |> mutate(metric = metric2_name))
  colnames(test.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  colnames(train.pred.df) <- c("DOLPHIN.ID", paste0("preds", 1:nfolds), "metric")
  
  return(list(testpreds = test.pred.df, trainpreds = train.pred.df))
  
}

run_probage = function(x, sample_info, n_cpg = 1000, use_weights = FALSE, controls_only = FALSE, reduce_size = FALSE, 
                       loglin_age = "lin", loglin_sexm = NULL){
  
  tmp <- try(
    run_probage_nocatch(x, sample_info, n_cpg = n_cpg, use_weights = use_weights, controls_only = controls_only, reduce_size = reduce_size,  
              loglin_age = loglin_age, loglin_sexm = loglin_sexm), 
    silent = TRUE)
  if(class(tmp) != "try-error"){
    tmp = tmp
  } else {
    metric1_name = paste0("probage_mle_acc_", n_cpg)
    metric2_name = paste0("probage_mle_bias_", n_cpg)
    empty_df = data.frame(DOLPHIN.ID = sample_info$DOLPHIN.ID, setNames(as.list(rep(NA, nfolds)), paste0("preds", 1:nfolds)))
    tmp = list(testpreds = rbind(empty_df |> mutate(metric = metric1_name), empty_df |> mutate(metric = metric2_name)),
               trainpreds = rbind(empty_df |> mutate(metric = metric1_name), empty_df |> mutate(metric = metric2_name)))
  }
  tmp
}

# m3 <- run_probage(x = x, sample_info = sample_info, use_weights = FALSE, controls_only = FALSE, reduce_size = TRUE)
