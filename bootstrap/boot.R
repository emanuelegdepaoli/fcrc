load("data_fullCOD.RData")
library(fcrc)
library(ggplot2)
causes
causes = as.vector(sapply(causes, function(x) gsub("_", " ", x)))
countries
p1 = 7
p2 = 9
p3 = 12
p4 = 12
p = p1+p2+p3+p4
p_vec = c(p1, p2, p3, p4)
L = comp_L(p_vec)
n = dim(Z_M)[3]
t = range_years
n_t = length(t)

ncores = 7
nrep = 500

# bootstrap M
set.seed(211)
resM = array(1, dim = c(nrep, p))
for(i in 1:nrep) {
  idx_boot = sample(1:n, n, T)
  Z_M_boot = Z_M[,,idx_boot]
  Y_M_boot = Y_M[,idx_boot]
  res_cv = cv_cgl(Y_M_boot, Z_M_boot, L, t, n, 4:6, ncores = ncores)
  cv_obj = plot_cv(res_cv)
  k = cv_obj$k_min1se
  lambda = exp(cv_obj$min1se$log_lambda)
  mat = mat_comp(Z_M_boot, Y_M_boot, t, k)
  res_best = alm_cgl_path(mat, L, lambda, eps = 1e-10, abs_tol_int = 1e-8, rel_tol_int = 1e-6)
  resM[i,res_best$index_null == 1] = 0
  print(i)
}

# bootstrap F
set.seed(233)
resF = array(1, dim = c(nrep, p))
for(i in 1:nrep) {
  idx_boot = sample(1:n, n, T)
  Z_F_boot = Z_F[,,idx_boot]
  Y_F_boot = Y_F[,idx_boot]
  res_cv = cv_cgl(Y_F_boot, Z_F_boot, L, t, n, 4:6, ncores = ncores)
  cv_obj = plot_cv(res_cv)
  k = cv_obj$k_min1se
  lambda = exp(cv_obj$min1se$log_lambda)
  mat = mat_comp(Z_F_boot, Y_F_boot, t, k)
  res_best = alm_cgl_path(mat, L, lambda, eps = 1e-10, abs_tol_int = 1e-8, rel_tol_int = 1e-6)
  resF[i,res_best$index_null == 1] = 0
  print(i)
}

saveRDS(resM, "bootM_loocv_40causes.rds")
saveRDS(resF, "bootF_loocv_40causes.rds")

# analysis
resM_boot = readRDS("bootM_loocv_40causes_500rep.rds")
resF_boot = readRDS("bootF_loocv_40causes_500rep.rds")
load("analysis_25countries_40causes.RData")
# reproducing idx bootstrap
set.seed(211)
idx_M = array(1, dim = c(nrep, n))
for(i in 1:nrep) {
  idx_boot = sample(1:n, n, T)
  idx_M[i,] = idx_boot
}

table(unlist(apply(idx_M[resM_boot[, causes == "5-39 EXT"] == 0,], 1, function(x) countries[!1:n%in%x])))/length(resM_boot[, causes == "5-39 EXT"] == 0)
apply(idx_M[resM_boot[, causes == "40-64 CIRC"] == 0,], 1, function(x) countries[!1:n%in%x])


# M
matM_k4 = mat_comp(Z_M, Y_M, t, 4)
cvM_obj = plot_cv(resM_cv_loo)
resM = alm_cgl_path(matM_k4, L, exp(cvM_obj$min1se$log_lambda), eps = 1e-10, abs_tol_int = 1e-8, rel_tol_int = 1e-6)
dfM = data.frame(group = rep(factor(!(1:p)%in%which(resM$index_null == 1), labels = c("not selected", "selected")), each = nrep),
                 Category = factor(rep(causes, each = nrep)))
dfM = dfM[as.vector(resM_boot) != 0,]
pdf("plot/bootM_500rep.pdf")
ggplot(data = dfM) + geom_bar(aes(x = ..count../nrep, y = Category, fill = group), show.legend = F) +
  labs(y = "Causes", x = "Proportion of selection") + scale_x_continuous(breaks = c(0.6, 0.7, 0.8, 0.9)) +
  scale_fill_manual("legend", values = c("not selected" = "black", "selected" = "grey")) + theme_bw() + 
  geom_vline(xintercept = 0.7, linetype = "dashed", col = "grey")
dev.off()

# F
matF_k4 = mat_comp(Z_F, Y_F, t, 4)
cvF_obj = plot_cv(resF_cv_loo)
resF = alm_cgl_path(matF_k4, L, exp(cvF_obj$min1se$log_lambda), eps = 1e-10, abs_tol_int = 1e-8, rel_tol_int = 1e-6)
dfF = data.frame(group = rep(factor(!(1:p)%in%which(resF$index_null==1), labels = c("not selected", "selected")), each = nrep),
                 Category = factor(rep(causes, each = nrep)))
dfF = dfF[as.vector(resF_boot) != 0,]
pdf("plot/bootF_500rep.pdf")
ggplot(data = dfF) + geom_bar(aes(x = ..count../nrep, y = Category, fill = group), show.legend = F) +
  labs(y = "Causes", x = "Proportion of selection") + scale_x_continuous(breaks = c(0.6, 0.7, 0.8, 0.9)) +
  scale_fill_manual("legend", values = c("not selected" = "black", "selected" = "grey")) + theme_bw() + 
  geom_vline(xintercept = 0.7, linetype = "dashed", col = "grey")
dev.off()

# M & F
df = rbind(dfM, dfF)
df$Sex = factor(c(rep("M", nrow(dfM)), rep("F", nrow(dfF))))
pdf("plot/boot_500rep.pdf")
ggplot(data = df) + geom_bar(aes(x = ..count../nrep, y = Category, fill = group), show.legend = F) +
  facet_wrap(~Sex, ) + labs(y = "Causes", x = "Proportion of selection") + scale_x_continuous(breaks = c(0.6, 0.7, 0.8, 0.9)) +
  scale_fill_manual("legend", values = c("not selected" = "black", "selected" = "grey")) + theme_bw() + 
  geom_vline(xintercept = 0.7, linetype = "dashed", col = "grey")
dev.off()

table(dfM$group, dfM$Category)[2,][order(table(dfM$group, dfM$Category)[2,])]
table(dfF$group, dfF$Category)[2,][order(table(dfF$group, dfF$Category)[2,])]
