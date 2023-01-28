#' Plot the results of the n-fold cross-validation
#' 
#'
#' @param results A list from \url{cv_cgl}.
#' @import ggplot2
#' @export
plot_cv = function(results){
  
  lambda_k = results$lambda_k
  k_range = results$k_range
  n_k = length(k_range)
  n = length(results$reord)
  mean = NULL
  se = NULL
  npred = NULL
  apmse_min = array(0, dim = c(n_k, 2))
  colnames(apmse_min) = c("log_lambda", "val")
  rownames(apmse_min) = k_range
  min = Inf
  idx_min = Inf
  idx = rep(k_range, each = 100)
  n_fold = max(results$fold)
  apsse = results$apsse
  apmse = array(sapply(1:n_fold, function(i) apsse[,,i]/sum(results$folds == i)),
                dim = c(100, n_k, n_fold))
  for(j in 1:n_k){
    mean_j = rowSums(apsse[,j,])/n
    mean = c(mean, mean_j)
    se = c(se, apply(apmse[,j,], 1, sd)/sqrt(n_fold))
    npred_j = results$res_data[[j]]$npred
    npred = c(npred, npred_j)
    idx_min_j = which(mean_j == min(mean_j))
    min_j = mean_j[idx_min_j]
    apmse_min[j,] = c(log(lambda_k[idx_min_j, j]), min_j)
    if(min_j < min){
      min = min_j
      npred_min = npred_j[idx_min_j]
    }
  }
  idx_min1se = which((mean < min+se[which(mean == min)]) & (npred <= npred_min))[1]
  k_min1se = idx[idx_min1se]
  k_min = k_range[apmse_min[,2] == min(apmse_min[,2])]
  min1se = list(log_lambda = log(as.vector(lambda_k)[idx_min1se]), val = mean[idx_min1se])
  min = list(log_lambda = apmse_min[apmse_min[,2] == min(apmse_min[,2]), 1], 
             val = apmse_min[apmse_min[,2] == min(apmse_min[,2]), 2])
  
  df = data.frame(mean = mean, lower = mean-se,
                  upper = mean+se, 
                  log_lambda = as.vector(log(lambda_k)),
                  label = idx)
  df2 = data.frame(log_lambda = apmse_min[,1], text_pos = max(mean),
                   label = k_range, 
                   text = "min")
  plot = ggplot(df, aes(x = log_lambda, y = mean, group = label)) + 
    xlab("log(lambda)") + ylab("APMSE") +
    facet_wrap(~label) +
    geom_errorbar(aes(ymin = lower, ymax = upper), col = "gray60") + 
    geom_point(col = "red") +
    geom_vline(data = df2, mapping = aes(xintercept = log_lambda), linetype = "dashed") +
    geom_text(data = df2, mapping = aes(x = log_lambda, y = text_pos, label = text),
              angle = 90, vjust = 1) +
    geom_vline(data.frame(x = min1se$log_lambda, label = k_min1se), 
               mapping = aes(xintercept = x), linetype = "dashed") +
    geom_text(data.frame(x = min1se$log_lambda, label = k_min1se, text = "min1se"), 
              mapping = aes(x = min1se$log_lambda, y = max(mean), label = text),
              angle = 90, vjust = 1) 
  
  return(list(apmse_min = apmse_min, k_min1se = k_min1se, k_min = k_min,
              min1se = min1se, min = min, plot = plot))
}

