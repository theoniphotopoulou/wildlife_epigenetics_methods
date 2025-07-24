get_performance_stats = function(valdata, testdata, opt_threshold = TRUE){
  
  # F1 score per par_set: threshold from validation, metrics from test
  f1_points_list <- list()
  par_sets <- unique(testdata$par_set)
  
  for (ps in par_sets) {
    val_sub <- valdata %>% filter(par_set == ps)
    test_sub <- testdata %>% filter(par_set == ps)
    
    # ROC on validation set
    roc_obj_val <- roc(response = val_sub$obs, predictor = val_sub$predy, quiet = TRUE, direction = "<")
    auc_val = as.numeric(auc(roc_obj_val))
    thresholds <- roc_obj_val$thresholds
    
    if(opt_threshold){
      
      best_f1 <- -Inf
      best_threshold <- NA
      
      # Find threshold with best F1 on validation set
      for (t in thresholds) {
        preds_val <- as.integer(val_sub$predy >= t)
        TP <- sum(preds_val == 1 & val_sub$obs == 1, na.rm = TRUE)
        FP <- sum(preds_val == 1 & val_sub$obs == 0, na.rm = TRUE)
        FN <- sum(preds_val == 0 & val_sub$obs == 1, na.rm = TRUE)
        
        precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
        recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA
        
        if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
          f1 <- 2 * precision * recall / (precision + recall)
          if (f1 > best_f1) {
            best_f1 <- f1
            best_threshold <- t
          }
        }
      }
      
    } else { 
      best_threshold = 0.5 
    }
    
    # Now apply the best threshold to the test set
    preds_test <- as.integer(test_sub$predy >= best_threshold)
    
    TP <- sum(preds_test == 1 & test_sub$obs == 1, na.rm = TRUE)
    FP <- sum(preds_test == 1 & test_sub$obs == 0, na.rm = TRUE)
    FN <- sum(preds_test == 0 & test_sub$obs == 1, na.rm = TRUE)
    TN <- sum(preds_test == 0 & test_sub$obs == 0, na.rm = TRUE)
    
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
    recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA
    f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA
    }
    
    acc <- (TP + TN) / (TP + TN + FP + FN)
    
    sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA
    
    # AUC on test set
    roc_obj_test <- roc(response = test_sub$obs, predictor = test_sub$predy, quiet = TRUE, direction = "<")
    auc_test = as.numeric(auc(roc_obj_test))
    
    f1_points_list[[ps]] <- tibble(
      par_set = ps,
      auc_val = auc_val,
      auc_test = auc_test,
      threshold = best_threshold,
      f1 = f1,
      acc = acc,
      sensitivity = sens,
      specificity = spec,
      precision = precision, 
      recall = recall,
      TP = TP,
      FN = FN,
      TN = TN,
      FP = FP
    )
  }
  
  # Combine and join with parameter table
  f1_points <- bind_rows(f1_points_list) 
  
  return(f1_points)
}