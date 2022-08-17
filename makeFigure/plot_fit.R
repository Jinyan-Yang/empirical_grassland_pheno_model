####################################
# plot fitted vs obs
####################################
# 
ym.18.df <- get.ym.func(18)
ym.drt.df <- get.ym.func('Drought')
sd.df <- summaryBy(GCC ~ Date,
                   data = ym.drt.df,FUN=c(sd),na.rm=TRUE,keep.names = F)
ym.18.df$GCC.sd <- sd.df$GCC.sd

gcc.met.con.df <- get.paddock.func('control')

pdf('figures/obs_pred.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec)) {
  
  # use different soil water cap and wilt for different site
  if(species.vec[i]=='ym'){
    df = ym.18.df
    # 
    # c.wd <- getwd()
    # setwd('c:/repo/dn_gcc/')
    ym.met.df <- readRDS('cache/ym/ym.met.rds')
    # 
    swc.ym.con <- quantile(ym.met.df$swc,na.rm=T,probs = c(0.01,0.99))
    swc.cap = round(swc.ym.con[[2]]*10)/10
    swc.wilt = round(swc.ym.con[[1]]*100)/100
    bucket.size=1000
  }else if(species.vec[i]=='flux'){
    df = gcc.met.con.df
    swc.q.con <- quantile(gcc.met.con.df$vwc,na.rm=T,probs = c(0.01,0.99))
    swc.cap =  0.3#round(swc.q.con[[2]]*10)/10
    swc.wilt = 0.01#round(swc.q.con[[1]]*100)/100
    bucket.size=1000
  }else{
    df = gcc.met.pace.df
    swc.cap = 0.13
    swc.wilt = 0.05
    bucket.size=300
  }
  # 
  par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
    f.t.opt = c(10,25,40,NA,NA,NA),
    f.extract = c(0.2,1.5,8,NA,NA,NA),
    f.sec = c(0.01,0.15,0.5,NA,NA,NA),
    f.growth = c(0.01,0.15,0.5,NA,NA,NA),
    q = c(0.1,3,15,NA,NA,NA),
    q.s = c(0.1,1,15,NA,NA,NA))
  row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
  
  # do ploting and prediction for v11
  plot.mcmc.func.2q(df,species.vec[i],
                    prep.in='Control',temp.in='Ambient',
                    my.fun = phenoGrass.func.v13,
                    nm.note='v1.2q.',use.smooth = TRUE,day.lag = 3,
                    swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,
                    bucket.size = bucket.size)
  plot.title.func(species.vec[i])
  
  
  #
  par.df <- data.frame(#f.h = c(200,220,240,NA,NA),
    f.t.opt = c(10,25,40,NA,NA,NA),
    f.extract = c(0.2,1.5,8,NA,NA,NA),
    f.sec = c(0.01,0.15,0.5,NA,NA,NA),
    f.growth = c(0.01,0.15,0.5,NA,NA,NA))
  row.names(par.df) <- c('min','initial','max','fit','stdv','prop')
# do the predict for v10
plot.mcmc.func.2q(df = df,
                  species.in=species.vec[i],
                  prep.in='Control',temp.in='Ambient',
                  my.fun = phenoGrass.func.v13,
                  nm.note='v0.',use.smooth = TRUE,
                  day.lag = 3,
                  swc.in.cap = swc.cap,swc.in.wilt = swc.wilt,
                  bucket.size = bucket.size,
                  q.s.in=0,q.in=1)

plot.title.func(species.vec[i])
}
dev.off()

# # This bit is not necessary but useful for debuging####
# #daisgnostic plot
# pdf('figures/diag.pdf',width = 8,height = 8*.618)
# for (i in seq_along(species.vec)) {
# 
#   fn <- sprintf('cache/smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
#   chain.3.ls = readRDS(fn)
#   lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i])
# 
#   par(mfrow=c(3,2),mar=c(5,5,1,1))
#   for(par.num in 1:6){
# 
#     start.row <- nrow(chain.3.ls[[1]]) / 4*3
# 
#     plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(start.row:nrow(chain.3.ls[[1]])))
# 
#   }
# }
# dev.off()

# ts with ci####
pdf('figures/v11_ts.pdf',width = 5*2,height = 5*5*.618)
par(mfrow=c(5,2))
par(mar=c(5,5,1,5))
for (i in seq_along(species.vec)) {
  fn <-  paste0('tmp/pred.smv1.2q.chain.',species.vec[i],'.Control.Ambient.rds')
  plot.ts.ci.func(fn)
  par(new=T)
  plot(0,pch=NA,ann=F,axes=F)
  legend('topleft',
         legend = sprintf('(%s) %s',letters[i],species.vec.nm[i]),
         bty='n')
}

dev.off()
