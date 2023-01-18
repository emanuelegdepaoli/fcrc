# Functional Concurrent Regression with Compositional Covariates 

This is the repository associated with the paper “[Functional Concurrent Regression with Compositional Covariates 
and its application to the time-varying effect of causes of death on human longevity](https://arxiv.org/abs/2301.06333) by E. G. Depaoli, M. Stefanucci and S. Mazzuco.

## Abstract
Multivariate functional data that are cross-sectionally compositional data are attracting increasing interest in the statistical modeling literature, a major example being trajectories over time of compositions derived from cause-specific mortality rates. In this work, we develop a novel functional concurrent regression model in which independent variables are functional compositions. This allows us to investigate the relationship over time between life expectancy at birth and compositions derived from cause-specific mortality rates of four distinct age classes, namely 0--4, 5--39, 40--64 and 65+. A penalized approach is developed to estimate the regression coefficients and select the relevant variables. Then an efficient computational strategy based on an augmented Lagrangian algorithm is derived to solve the resulting optimization problem. The good performances of the model in predicting the response function and estimating the unknown functional coefficients are shown in a simulation study. The results on real data confirm the important role of neoplasms and cardiovascular diseases in determining life expectancy emerged in other studies and reveal several other contributions not yet observed. 

## Contents

- the folder `fcrc` includes the R package implementing the methods.
- the folder `data` includes the data used for the analysis.
- the file `preprocessing.R` contains the code to preprocess the data.
- the file `analysis_25countries_40causes.R` contains the code to reproduce the analysis and the plots. 
- the folders `bootstrap` and `simulations` contain the code to reproduce the results of the simulations.
- the folder `plot` contain the plots showed in the paper


