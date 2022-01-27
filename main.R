#Welcome to grassland phenology modlling project#####################################
# thie repo fits the Choler-Hufkens to observed cover################################
#By Jinyan Yang (jyyoung@live.com)###################################################

# known issues####
# 1. I have not yet made the dependent packages install automatically;
# 2. fitting take a long time so may as well just use the fitted reasults 
# (i.e., skip step II in the following);
# 3. in the code model v1 is refered as v13 and v0 refered as v10
# This is because i can and I want to make whoever reads the code suffer =P
# It's acutally because various versions of model tested are not presented
# 4. there is a day.lay term which was used to determ the legency effect of water and T but
# the term is not currently useful; please ingore while using 
# what it does is simply ingnore the number of days of data

# Step I: load environment####
source('r/read_spc_nm.R')
source('r/load.R')

# Step II: fit model to data####
source('r/fit_all.R')#v13
source('r/fit_q1_qs0.R')#v10
# get best fits
source('r/get_best_fit.R')
# get Ci for fitting
source('r/get_precit_ci.R')

# Step III: make plots#####
#  check if par differ significantly
source('r/compare_chain.R')
# predict according to fitted param
source('makeFigure/plot_fit.R')

# fig 1 theoretical plot
source('makeFigure/hypotheses.R')
# fig 3 ts for v10
source('makeFigure/plot_v10_scatter.R')
# fig 4 ts and scatter v13
source('makeFigure/plot_fit_TS_scatter.R')
# fig 5 beta
source('makeFigure/fw_ambient.R')
# fig 6 par distributions
source('makeFigure/par_vioPlot.R')
# fig 7 brown-down thresholds
source('makeFigure/plot_threshold_ambient.R')
# fig 8 drought
# first make prediction
source('makeFigure/plot_predc_drt.R')
# then make bar plot
source('makeFigure/plot_drt.R')

# fig s1 ts for v10
source('makeFigure/plot_v10_ts.R')

# table 2  stats
source('r/get_bic.R')
# table s2 fitted par value
source('r/get_fittedPar.R')


