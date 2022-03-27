####################################
# plot fitted vs obs
####################################
source('r/read_spc_nm.R')
source('r/load.R')
# 
y.nm.vec <- c(expression(T[opt]),expression(r[extract]),
              expression(r[senescence]),expression(r[growth]),
              expression(q[growth]),expression(q[senescence]))

# This bit is not necessary but useful for debuging####
#daisgnostic plot
palette(col.df$cityNight[c(3:5)])
# palette(col.df$iris[c(1,1,2,2,2,3,3,3,4,4)])

pdf('figures/chain.V1.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec)) {
  fn <- sprintf('cache/smv1.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  chain.3.ls = readRDS(fn)

  # plot.check.mcmc.func(chain.3.ls[[1]],species.in=species.vec[i],
  #                      nm.vec = y.nm.vec,plot.par.no=1:6,col.in = t_col(palette()[1],percent = 50))
  plot.hist.func(chain.3.ls,nm.vec = y.nm.vec,plot.par.no=1:6)
  # lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i],
  #        nm.vec =y.nm.vec,plot.par.no=1:6,col.in = i)
  plot.title.func(species.vec.nm[i],where.to=0.5)
  # par(mfrow=c(3,2),mar=c(5,5,1,1))
  # for(par.num in 1:6){
  # 
  #   start.row <- nrow(chain.3.ls[[1]]) / 4*3
  # 
  #   plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(start.row:nrow(chain.3.ls[[1]])))
  # 
  # }
}
dev.off()
# 
palette(col.df$cityNight[c(3:5)])
pdf('figures/chain.V0.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec)) {
  fn <- sprintf('cache/smv0.chain.%s.Control.Ambient.rds',species.vec[i])
  chain.3.ls = readRDS(fn)
  plot.hist.func(chain.3.ls,nm.vec = y.nm.vec,plot.par.no=1:4)
  # lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i],
  #        nm.vec =y.nm.vec,plot.par.no=1:6,col.in = i)
  plot.title.func(species.vec.nm[i],where.to=0.5)
  # lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i],
  #        nm.vec =y.nm.vec,plot.par.no=1:4)
  # plot.title.func(species.vec.nm[i],where.to=0.5)
  # plot.title.func(paste0('(',letters[i],') ',species.vec.nm[i]),where.to=0)
  # par(mfrow=c(3,2),mar=c(5,5,3.5,1))
  # for(par.num in 1:6){
  #   
  #   start.row <- nrow(chain.3.ls[[1]]) / 4*3
  #   
  #   plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(start.row:nrow(chain.3.ls[[1]])),
  #                       nm.vec = y.nm.vec)
  # 
  # }
  # par(new=T)
  # par(mfrow=c(1,1),mar=c(5,5,3,1))
  # plot(1,pch=NA,ann=F,axes=F)
  # title(species.vec.nm[i],line=3.5)
}
dev.off()
