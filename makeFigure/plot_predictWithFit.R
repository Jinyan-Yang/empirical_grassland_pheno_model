day.lag <- 3
source('r/pace_data_process.R')
source('r/ym_data_process.R')
source('r/v13_common_fun.R')
source('models/hufkens/hufkensV13.R')
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
library(zoo)
source('r/plot.mcmc.r')

ym.df <- get.ym.func('Drought')
# gcc.met.con.df <- get.paddock.func('control')
# gcc.met.iri.df <- get.paddock.func('irrigated')
species.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','Pha','ym')
plot.nm.vec <- c('Bis','Dig','Luc','Fes','Rye','Kan','Rho','Pha','YM')
# species.vec <- 'Pha'
pdf('figures/makepredictions.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec)) {
  
  pred.con <- 'Drought'
  if(species.vec[i]=='ym'){
    df = ym.df
    # 
    c.wd <- getwd()
    setwd('c:/repo/dn_gcc/')
    ym.met.df <- readRDS('cache/ym.met.rds')
    setwd(c.wd)
    # 
    swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    swc.cap = round(swc.ym.con[[2]]*10)/10
    swc.wilt = round(swc.ym.con[[1]]*100)/100
    bucket.size=1000
  }else{
    df = gcc.met.pace.df
    swc.cap = 0.13
    swc.wilt = 0.05
    bucket.size=300
  }


  plot.mcmc.func.2q(df,species.vec[i],pred.con,'Ambient',
                    my.fun = phenoGrass.func.v13,
                    nm.note='smv13.2q.',use.smooth = TRUE,day.lag = 3,
                    swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,
                    bucket.size = bucket.size,do.predict = 'Control')
  plot.title.func(plot.nm.vec[i]) 
  
}
dev.off()

# # ###########
# plot.mcmc.func.2q(gcc.met.iri.df,'flux','Irrigated','Ambient',
#                   my.fun = phenoGrass.func.v13,
#                   nm.note='v13.2q',use.smooth = TRUE,day.lag = 3,
#                   swc.in.cap =  swc.q.iri[[2]],swc.in.wilt =  swc.q.iri[[1]],bucket.size = bucket.size,do.predict = 'Control')
# plot.title.func('Flux_irrigated')
# 
# # #
# # fn <- sprintf('cache/smv13.2qchain.%s.Irrigated.Ambient.rds','flux')
# # chain.3.ls = readRDS(fn)
# # lapply(chain.3.ls, plot.check.mcmc.func,species.in='Flux-irrigated')
# # 
# # par(mfrow=c(3,2),mar=c(5,5,1,1))
# # for(par.num in 1:6){
# # 
# #   start.row <- nrow(chain.3.ls[[1]]) / 4*3
# # 
# #   plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(start.row:nrow(chain.3.ls[[1]])))
# # 
# # }
# 
# dev.off()