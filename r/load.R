day.lag <- 3
source('r/pace_data_process.R')
source('r/ym_data_process.R')
source('r/v13_common_fun.R')
source('models/hufkens/hufkensV13.R')
source('r/functions_mcmc_V12.R')
source('r/process_paddock_gcc_met.R')
source('r/plot.mcmc.r')
source('r/function_get_ci.R')
source('r/get_bic.R')

library(zoo)
library(foreach)
library(doParallel)
library(DEoptim)
library(mvtnorm)
library(lubridate)
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")

