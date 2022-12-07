sim_nrep = 100
# n, p, ncomp
temp1 = matrix(c(50, 40, 1, 
                 50, 40, 4, 
                 50, 100, 4), ncol = 3, byrow = T)
corr_t = c(0.2, 0.6)
corr_X = c(0.2, 0.6)
snr = c(2, 4)
temp2 = expand.grid(corr_t, corr_X, snr)
grid_sim = do.call(rbind, apply(temp1, 1, function(x) cbind(matrix(rep(x, nrow(temp2)), ncol = ncol(temp1), byrow = T),
                                                  temp2)))
colnames(grid_sim) = c("n", "p", "ncomp", "corr_t", "corr_x", "snr")
res = lapply(1:24, function(i) readRDS(paste0("sim_setup",i,".rds")))
attr(res, "info")

# round keep zeros
myround = function(x, digit) sapply(x, function(x) format(round(x, digits = digit), nsmall = digit))
dol = function(x) paste(rep("$", length(x)), x, rep("$", length(x)))

#order table paper
#$\rho_X$ & $\rho_T$ & $n$ & $p$ & $q$  &    & CGL     & BGL    & Naive&  & CGL     & BGL    & Naive 
# SNR = 2, FNR, FPR
str2print = ""
for(i in corr_X){
  for(j in corr_t){
    count = 0
    for(l in 1:nrow(temp1)){
      row = c(temp1[l,], j, i, 2)
      idx = which(apply(grid_sim, 1, function(x) sum(x == row)== length(row)))
      # if(idx > 16){ # for partial results, omit when simulation will be finished
      #   str_fin = paste(paste(rep(" ", 11), rep(" & ", 11), collapse = ""), "\\\\", collapse = "") # \\ to print \
      # } else {
        if(count == 0){
          str1 = paste(c(dol(grid_sim[idx, c("corr_x", "corr_t", "n","p", "ncomp")]), ""), "& ", collapse = "")
        } else {
          str1 = paste(c(" ", " ", grid_sim[idx, c("n","p", "ncomp")], ""), "& ", collapse = "")
        }
        
        #FPR
        mean = dol(myround(apply(res[[idx]], 3, colMeans)["FPR_t",], 2))
        se = myround(apply(res[[idx]], c(2, 3), sd)["FPR_t",]/sqrt(sim_nrep), 2)
        zeros = (apply(res[[idx]], 3, colMeans)["FPR_t",] == 0)
        se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
        str_se =  dol(paste0(rep("(", 3), se, ")"))
        str2 = paste(mean, str_se, "& ", collapse = "")
        #FNR
        mean = dol(myround(apply(res[[idx]], 3, colMeans)["FNR_t",], 2))
        se = myround(apply(res[[idx]], c(2, 3), sd)["FNR_t",]/sqrt(sim_nrep), 2)
        zeros = (apply(res[[idx]], 3, colMeans)["FNR_t",] == 0)
        se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
        str_se =  dol(paste0(rep("(", 3), se, ")"))
        str3 = paste(mean,str_se, c(rep("& ", length(c(mean))-1), "\\\\"), collapse = "")
        str_fin = paste(c(str1, str2, " & ", str3), collapse = "")
      # }
      
      str2print = paste(str2print, str_fin, sep = " \n")
      count = count + 1
    }
  }
}
cat(str2print)

# SNR = 2,est error, pred error
str2print = ""
for(i in corr_X){
  for(j in corr_t){
    count = 0
    for(l in 1:nrow(temp1)){
      row = c(temp1[l,], j, i, 2)
      idx = which(apply(grid_sim, 1, function(x) sum(x == row)== length(row)))
      # if(idx > 16){ # for partial results, omit when simulation will be finished
      #   str_fin = paste(paste(rep(" ", 11), rep(" & ", 11), collapse = ""), "\\\\", collapse = "") # \\ to print \
      # } else {
        if(count == 0){
          str1 = paste(c(dol(grid_sim[idx, c("corr_x", "corr_t", "n","p", "ncomp")]), ""), "& ", collapse = "")
        } else {
          str1 = paste(c(" ", " ", grid_sim[idx, c("n","p", "ncomp")], ""), "& ", collapse = "")
        }
        
        # pred error
        mean = dol(myround(apply(res[[idx]], 3, colMeans)["pred_error_test",], 2))
        se = myround(apply(res[[idx]], c(2, 3), sd)["pred_error_test",]/sqrt(sim_nrep), 2)
        zeros = (apply(res[[idx]], 3, colMeans)["pred_error_test",] == 0)
        se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
        str_se =  dol(paste0(rep("(", 3), se, ")"))
        str2 = paste(mean, str_se, "& ", collapse = "")
        # est error * 100
        mean = dol(myround(apply(res[[idx]]*100, 3, colMeans)["est_error",], 2))
        se = myround(apply(res[[idx]]*100, c(2, 3), sd)["est_error",]/sqrt(sim_nrep), 2)
        zeros = (apply(res[[idx]]*100, 3, colMeans)["est_error",] == 0)
        se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
        str_se =  dol(paste0(rep("(", 3), se, ")"))
        str3 = paste(mean,str_se, c(rep("& ", length(c(mean))-1), "\\\\"), collapse = "")
        str_fin = paste(c(str1, str2, " & ", str3), collapse = "")
#      }
      
      str2print = paste(str2print, str_fin, sep = " \n")
      count = count + 1
    }
  }
 }
cat(str2print)


# SNR = 4, FNR, FPR
str2print = ""
for(i in corr_X){
  for(j in corr_t){
    count = 0
    for(l in 1:nrow(temp1)){
      row = c(temp1[l,], j, i, 4)
      idx = which(apply(grid_sim, 1, function(x) sum(x == row)== length(row)))
      # if(idx > 16){ # for partial results, omit when simulation will be finished
      #   str_fin = paste(paste(rep(" ", 11), rep(" & ", 11), collapse = ""), "\\\\", collapse = "") # \\ to print \
      # } else {
      if(count == 0){
        str1 = paste(c(dol(grid_sim[idx, c("corr_x", "corr_t", "n","p", "ncomp")]), ""), "& ", collapse = "")
      } else {
        str1 = paste(c(" ", " ", grid_sim[idx, c("n","p", "ncomp")], ""), "& ", collapse = "")
      }
      
      #FPR
      mean = dol(myround(apply(res[[idx]], 3, colMeans)["FPR_t",], 2))
      se = myround(apply(res[[idx]], c(2, 3), sd)["FPR_t",]/sqrt(sim_nrep), 2)
      zeros = (apply(res[[idx]], 3, colMeans)["FPR_t",] == 0)
      se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
      str_se =  dol(paste0(rep("(", 3), se, ")"))
      str2 = paste(mean, str_se, "& ", collapse = "")
      #FNR
      mean = dol(myround(apply(res[[idx]], 3, colMeans)["FNR_t",], 2))
      se = myround(apply(res[[idx]], c(2, 3), sd)["FNR_t",]/sqrt(sim_nrep), 2)
      zeros = (apply(res[[idx]], 3, colMeans)["FNR_t",] == 0)
      se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
      str_se =  dol(paste0(rep("(", 3), se, ")"))
      str3 = paste(mean,str_se, c(rep("& ", length(c(mean))-1), "\\\\"), collapse = "")
      str_fin = paste(c(str1, str2, " & ", str3), collapse = "")
      # }
      
      str2print = paste(str2print, str_fin, sep = " \n")
      count = count + 1
    }
  }
}
cat(str2print)

# SNR = 2,est error, pred error
str2print = ""
for(i in corr_X){
  for(j in corr_t){
    count = 0
    for(l in 1:nrow(temp1)){
      row = c(temp1[l,], j, i, 4)
      idx = which(apply(grid_sim, 1, function(x) sum(x == row)== length(row)))
      # if(idx > 16){ # for partial results, omit when simulation will be finished
      #   str_fin = paste(paste(rep(" ", 11), rep(" & ", 11), collapse = ""), "\\\\", collapse = "") # \\ to print \
      # } else {
      if(count == 0){
        str1 = paste(c(dol(grid_sim[idx, c("corr_x", "corr_t", "n","p", "ncomp")]), ""), "& ", collapse = "")
      } else {
        str1 = paste(c(" ", " ", grid_sim[idx, c("n","p", "ncomp")], ""), "& ", collapse = "")
      }
      
      # pred error
      mean = dol(myround(apply(res[[idx]], 3, colMeans)["pred_error_test",], 2))
      se = myround(apply(res[[idx]], c(2, 3), sd)["pred_error_test",]/sqrt(sim_nrep), 2)
      zeros = (apply(res[[idx]], 3, colMeans)["pred_error_test",] == 0)
      se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
      str_se =  dol(paste0(rep("(", 3), se, ")"))
      str2 = paste(mean, str_se, "& ", collapse = "")
      # est error * 100
      mean = dol(myround(apply(res[[idx]]*100, 3, colMeans)["est_error",], 2))
      se = myround(apply(res[[idx]]*100, c(2, 3), sd)["est_error",]/sqrt(sim_nrep), 2)
      zeros = (apply(res[[idx]]*100, 3, colMeans)["est_error",] == 0)
      se[as.numeric(se)== 0 & !zeros] = "<10^{-2}"
      str_se =  dol(paste0(rep("(", 3), se, ")"))
      str3 = paste(mean,str_se, c(rep("& ", length(c(mean))-1), "\\\\"), collapse = "")
      str_fin = paste(c(str1, str2, " & ", str3), collapse = "")
      #      }
      
      str2print = paste(str2print, str_fin, sep = " \n")
      count = count + 1
    }
  }
}
cat(str2print)

