source('r/plot.mcmc.r')
source('r/read_spc_nm.R')
# 
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
library(doBy)
library(lubridate)

# ts with ci####
pdf('figures/v10_ts.pdf',width = 5*2,height = 5*5*.618)
par(mfrow=c(5,2))
par(mar=c(5,5,1,5))
for (i in seq_along(species.vec)) {
  fn <-  paste0('tmp/pred.smv0.chain.',species.vec[i],'.Control.Ambient.rds')
  plot.ts.ci.func(fn)
 
}

dev.off()
