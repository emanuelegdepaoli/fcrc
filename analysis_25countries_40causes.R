load("data_fullCOD.RData")
library(fcrc)
causes
countries
q = 4
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

# M 
resM_cv_loo = cv_cgl(Y_M, Z_M, L, t, n, 4:7)
plotM_obj_loocv = plot_cv(resM_cv_loo)

matM_k4 = mat_comp(Z_M, Y_M, t, 4)

resM_loo = alm_cgl_path(matM_k4, L, exp(plotM_obj_loocv$min1se$log_lambda), 
                   eps = 1e-10, abs_tol_int = 1e-8, rel_tol_int = 1e-6)

# F 
resF_cv_loo = cv_cgl(Y_F, Z_F, L, t, n, 4:7)
plotF_obj_loocv = plot_cv(resF_cv_loo)

matF_k4 = mat_comp(Z_F, Y_F, t, 4)

resF_loo = alm_cgl_path(matF_k4, L, exp(plotF_obj_loocv$min1se$log_lambda), 
                   eps = 1e-10, abs_tol_int = 1e-8, rel_tol_int = 1e-6)

# plot coefficients 
library(ggplot2)
library(ggrepel)
library(dplyr)
pred_names = sapply(causes, function(x) gsub(".*_", "", x))
comp_names =  c("Age class 0-4", "Age class 5-39", "Age class 40-64", "Age class 65+")
idxM = (1:matM_k4$p)[-which(resM_loo$index_null==1)]
idxF = (1:matF_k4$p)[-which(resF_loo$index_null==1)]
idx_comp = unlist(sapply(1:length(p_vec), function(x) rep(x, each = p_vec[x])))
B_estM = matrix(resM_loo$beta, ncol = matM_k4$k, byrow = T)
Bc_estM = matrix(resM_loo$betac, ncol = matM_k4$k, byrow = T)
B_estF = matrix(resF_loo$beta, ncol = matF_k4$k, byrow = T)
Bc_estF = matrix(resF_loo$betac, ncol = matF_k4$k, byrow = T)
curvesM = matM_k4$Phi%*%t(B_estM[idxM,])
curvesF = matF_k4$Phi%*%t(B_estF[idxF,])
dfM_coef = data.frame(Year = rep(t, length(idxM)), 
                      Causes = rep(pred_names[idxM], each = n_t),
                      Coefficient = matrix(curvesM, ncol = 1),
                      Composition = rep(idx_comp[idxM], each = n_t))
dfF_coef = data.frame(Year = rep(t, length(idxF)), 
                      Causes = rep(pred_names[idxF], each = n_t),
                      Coefficient = matrix(curvesF, ncol = 1),
                      Composition = rep(idx_comp[idxF], each = n_t))
df_coef = rbind(dfM_coef, dfF_coef)
df_coef$Sex = as.factor(c(rep("M", nrow(dfM_coef)), rep("F", nrow(dfF_coef))))

# labels all causes 
df_coef = df_coef  %>%
  mutate(lab = if_else(Year == max(Year) & Composition == 3,
                       as.character(Causes), NA_character_))
df_coef = df_coef  %>%
  mutate(lab = if_else(Year == max(Year) & Composition == 4, 
                       as.character(Causes), lab))
df_coef = df_coef  %>%
  mutate(lab = if_else(Year == max(Year) & Composition == 1,
                       as.character(Causes), lab))
df_coef = df_coef  %>%
  mutate(lab = if_else(Year == max(Year) & Composition == 2,
                       as.character(Causes), lab))
df_coef$Composition = factor(df_coef$Composition,levels = 1:4, labels =  comp_names)

ageclass = c("0-4", "5-39","40-64", "65+")
library(ggrepel)
library(RColorBrewer)
brewer.pal.info[brewer.pal.info$category == 'qual',]
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
names(mycolors) = unique(df_coef$Causes)

for(i in 1:4) {
  pdf(paste0("plot/coef", ageclass[i], ".pdf"), height=6, width=10)
  obj = ggplot(df_coef %>% subset(Composition == comp_names[i]), aes(x = Year, y = Coefficient, colour = Causes)) +
    geom_line() + 
    geom_hline(yintercept = 0, lty = 2) + ggtitle(comp_names[i]) + 
    facet_wrap(~ Sex) +
    scale_color_manual(values = mycolors) + 
    geom_label_repel(aes(label = lab),
                     nudge_x = 1,
                     na.rm = TRUE, show.legend = F) +
    theme(legend.position="none") 
  plot(obj)
  dev.off()
}

# single plots by sex and causes
for(i in 1:4) {
  for(sex in c("M", "F")){
    pdf(paste0("plot/coef", ageclass[i],sex, ".pdf"), height=6, width=6)
    obj = ggplot(df_coef %>% subset(Composition ==  comp_names[i] & Sex == sex), 
                 aes(x = Year, y = Coefficient, colour = Causes)) +
      geom_line() + 
      geom_hline(yintercept = 0, lty = 2) + ggtitle(paste(comp_names[i], sex)) + 
      scale_color_manual(values = mycolors) + 
      geom_label_repel(aes(label = lab),
                       nudge_x = 1,
                       na.rm = TRUE, show.legend = F) +
      theme(legend.position="none") 
    plot(obj)
    dev.off()
  }
}

for(sex in c("M", "F")){
  pdf(paste0("plot/coef", sex, ".pdf"), height=6, width=6)
  obj = ggplot(df_coef %>% subset(Sex == sex), 
               aes(x = Year, y = Coefficient, colour = Causes)) +
    geom_line() + facet_wrap(~ Composition, scales = "free") +
    geom_hline(yintercept = 0, lty = 2) + 
    scale_color_manual(values = mycolors) + 
    geom_label_repel(aes(label = lab),
                     nudge_x = 1,
                     na.rm = TRUE, show.legend = F) +
    theme(legend.position="none") 
  plot(obj)
  dev.off()
}

# fitted curves
plot_Y(resM_loo, matM_k4, countries, Y_M, Z_M, t, sep = F, devnew = T)
plot_Y(resF_loo, matM_k4, countries, Y_F, Z_F, t, sep = F, devnew = T)

# in-sample error 
error_comp(resM_loo$beta, resM_loo$beta0, Y_M, Z_M, matM_k4$Phi, t)
error_comp(resF_loo$beta, resF_loo$beta0, Y_F, Z_F, matF_k4$Phi, t)

# importance plots
library(fda)
library(RColorBrewer)
estM_fdobj = fd(coef = t(B_estM), 
                basisobj = create.bspline.basis(rangeval = c(min(t), max(t)), nbasis = 4))
energy_matM = matrix(0, ncol = p, nrow = n_t)
for(i in 1:(n_t-1)) energy_matM[i+1,] = diag(inprod(estM_fdobj, estM_fdobj, rng = c(min(t)+i-1, min(t)+i)))
idx_start = c(1, cumsum(head(p_vec, -1)) + 1)
idx_end = cumsum(p_vec)
contrM = sapply(1:q, function(i) rowSums(energy_matM[,idx_start[i]:idx_end[i]]))
contrM_norm = contrM[-1,]/rowSums(contrM[-1,])
data_toplotM = data.frame(cum_energy = as.vector(contrM_norm), 
                          Class = rep(c("0-4", "5-39", "40-64", "65+"), each = n_t-1),
                          Year = rep(1966:2012, q))
estF_fdobj = fd(coef = t(B_estF), 
                basisobj = create.bspline.basis(rangeval = c(min(t), max(t)), nbasis = 4))
energy_matF = sapply(t, function(x) diag(inprod(estF_fdobj, estF_fdobj, rng = c(min(t), x))))
idx_start = c(1, cumsum(head(p_vec, -1)) + 1)
idx_end = cumsum(p_vec)
contrF = sapply(1:q, function(i) colSums(energy_matF[idx_start[i]:idx_end[i],]))
contrF_norm = contrF[-1,]/rowSums(contrF[-1,])
data_toplotF = data.frame(cum_energy = as.vector(contrF_norm), 
                          Class = rep(c("0-4", "5-39", "40-64", "65+"), each = n_t-1),
                          Year = rep(1966:2012, q))

data_toplot = rbind(data_toplotM, data_toplotF)
data_toplot$Sex = rep(c("M", "F"), each = nrow(data_toplotM))
pdf("plot/rel_energy.pdf", height = 6, width = 10)
ggplot() + geom_bar(aes(y = cum_energy, x = Year, fill = Class, group = Sex), data = data_toplot,
                    stat="identity",  width = 1) + ylab("Relative magnitude") + facet_wrap(~Sex) +
  scale_fill_discrete(name = "Age classes") + scale_y_continuous(n.breaks = 10) +theme_bw()  +
  scale_fill_brewer(palette = "Accent")
dev.off()

# save results
save(list=c("resM_cv_loo", "resF_cv_loo"), 
     file = "analysis_25countries_40causes.RData")
