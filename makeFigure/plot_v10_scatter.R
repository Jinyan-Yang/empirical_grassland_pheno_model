source('r/read_spc_nm.R')
# 
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")

palette(c(col.df$iris,col.df$daisy))
# plot
png('figures/obs_fit_v10_scatter.png',height = 600,width = 600)
par(mar=c(5,5,5,5))
for (i in seq_along(species.vec)){
  fn <- sprintf('tmp/pred.smv13.q1.qs0.chain.%s.Control.Ambient.rds',species.vec[i])
  hufken.pace.pred <- readRDS(fn)
  # 
  ci.fm <- sprintf('tmp/ci.smv13.q1.qs0.chain.%s.Control.Ambient.rds',species.vec[i])
  ci.m <- readRDS(ci.fm)
  # hufken.pace.pred$cover.05 <- ci.m[1,]
  # hufken.pace.pred$cover.95 <- ci.m[2,]
  hufken.pace.pred$cover.50 <- ci.m[3,]
  # 
  if(i == 1){
    plot(GCC.norm~cover.50,data = hufken.pace.pred,
         xlim=c(0,1),ylim=c(0,1),
         xlab='Modelled cover',ylab = 'Observed cover',pch=16,col=i)
    
  }else{
    points(GCC.norm~cover.50,data = hufken.pace.pred,
           xlim=c(0,1),ylim=c(0,1),
           pch=16,col=i)
  }
  legend('bottomright',legend = species.vec.nm,col=palette(),
         pch=16,bty='n')
  legend('topleft',legend = '(b)',bty='n',col = 2)
  abline(a=0,b=1,lty='dashed',col='grey',lwd=2)
  
}
dev.off()