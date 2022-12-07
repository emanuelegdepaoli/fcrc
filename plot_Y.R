# plot functions

plot_Y = function(results, matrices, names, Y, Z, t, Zc = NULL, sep = F, devnew = F){
  n = dim(Z)[3] # observations to plot can be different from matrices$n
  n_t = length(t)
  p = matrices$p
  k = matrices$k
  if(is.null(Zc)) Zc = array(1, dim = c(n_t, 1, n))
  B_est = matrix(results$beta, ncol = k, byrow = T)
  Bc_est = matrix(results$betac, ncol = k, byrow = T)
  Y_est = apply(Z, 3, function(x) rowSums(x*matrices$Phi%*%t(B_est))) + 
    apply(Zc, 3, function(x) rowSums(x*matrices$Phi%*%t(Bc_est)))
  # estimate null model
  Bc_null = matrix(solve(matrices$M)%*%matrices$P, ncol = k, byrow = T)
  Y_est_null = apply(Zc, 3, function(x) rowSums(x*matrices$Phi%*%t(Bc_null)))
  df_est = data.frame(year = rep(t, n), 
                      names = rep(names, each = n_t),
                      Lifexp = matrix(Y_est, ncol = 1),
                      label = "est")
  df_est_null = data.frame(year = rep(t, n), 
                           names = rep(names, each = n_t),
                           Lifexp = matrix(rep(Y_est_null, n), ncol = 1),
                           label = "est_null")
  df_obs = data.frame(year = rep(t, n), 
                      names = rep(names, each = n_t),
                      Lifexp = matrix(Y, ncol = 1),
                      label = "obs")
  
  obj = ggplot(df_est) + facet_wrap(~ names) + 
    geom_point(mapping = aes(x = year, y = Lifexp, group = label, colour = label), 
               data = df_obs, size = 0.3) + 
    geom_line(mapping = aes(x = year, y = Lifexp, group = label, colour = label),
              data = df_est) +
    geom_line(mapping = aes(x = year, y = Lifexp, group = label, colour = label),
              data = df_est_null)
  
  list = lapply(1:length(names), function(i)
    ggplot(df_est[df_est$names == names[i],]) + 
      geom_point(mapping = aes(x = year, y = Lifexp, group = label, colour = label), 
                 data = df_obs[df_obs$names == names[i],], size = 0.3) + 
      geom_line(mapping = aes(x = year, y = Lifexp, group = label, colour = label),
                data = df_est[df_est$names == names[i],]) +
      geom_line(mapping = aes(x = year, y = Lifexp, group = label, colour = label),
                data = df_est_null[df_est_null$names == names[i],]) +
      ggtitle(names[i])
  )
  if (sep) {
    for(i in 1:length(list)) {
      if (devnew) {
        dev.new()
      } 
      plot(list[[i]])
    }
  } else {
    obj
  }
}