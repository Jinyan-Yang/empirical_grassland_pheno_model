#Welcome to grassland phenology modlling project#####################################
# thie repo fits the Choler-Hufkens to observed cover################################
#By Jinyan Yang (jyyoung@live.com)###################################################

# known issues####
# 1. I have not yet made the dependent packages install automatically;
# 2. fitting take a long time so may as well just use the fitted reasults 
# (i.e., skip step II in the following);


# Step I: load environment####
source('r/read_spc_nm.R')
source('r/load.R')

# Step II: fit model to data####
source('r/fit_all.R')
# get Ci for fitting
source('r/get_precit_ci.R')

# Step III: make plots#####
source('makeFigure/plot_fit.R')


