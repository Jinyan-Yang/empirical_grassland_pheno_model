source('r/read_spc_nm.R')
spc.vec <- species.vec[3:9]#[3:9]
# prepare significance data####
# read in significant data

# make plot####
library(vioplot)
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
palette(c(col.df$iris))
col.nm.vec <- c(1,1,2,2,2,3,3,3,4,4)

plot.violin.func <- function(fn.crtl,fn.drt,burin.frac=0.75,log.y=F,y.range=F,spc.nm){
  
  # tmp.ls <- list()
  # for (i in seq_along(spc.vec)) {
   
    # fn <- sprintf('cache/smv13.qs1.chain.%s.Control.Ambient.rds',spc.vec[i])

    in.chain =  readRDS(fn.crtl)
    burnIn = 1
    chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round(nrow(m.in)* (burin.frac)):nrow(m.in),])

    tmp.ls <- do.call(rbind,chain.3.ls.new)
    # 
    in.chain.drt =  readRDS(fn.drt)
    burnIn = 1
    chain.3.ls.drt = lapply(in.chain.drt,function(m.in)m.in[round(nrow(m.in)* (burin.frac)):nrow(m.in),])
    
    tmp.ls.drt <- do.call(rbind,chain.3.ls.drt)
    
    # fn.1 <- sprintf('cache/v13.2q.chain.%s.bestfit.rds',spc.vec[i])
    # chain.fes <- readRDS(fn.1)
    # tmp.ls[[i]] <- data.frame(spc = species.vec.nm[i],
    #                           par.val = chain.fes[,col2plot])
    
  # }
  
  plot.df <- tmp.ls
  plot.df$trt <- 'ctrl'
  
  plot.df.drt <- tmp.ls.drt
  plot.df.drt$trt <- 'drt'
  
  plot.all.df <- rbind(plot.df,plot.df.drt)
  # plot.df$spc[plot.df$spc =='flux'] <- 'Flux Tower'
  # plot.all.df$trt <- as.factor(plot.df$trt)
  # plot.all.df[plot.all.df$trt=='ctrl',]

  # the if statement is not used for now
  par.vec <- names(plot.all.df)
  par(mfrow=c(3,2))
  for (i in 1:6) {
    tmp.df <- plot.all.df[,c(par.vec[i],'trt')]
    names(tmp.df) <- c('par.val','trt')
    if(i %in% 5:6){
      vioplot(par.val~trt,tmp.df,col=col.nm.vec,
              xlab='',ylab=par.vec[i],las = 2,ylog = T)
    }else{
      vioplot(par.val~trt,tmp.df,col=col.nm.vec,
              xlab='',ylab=par.vec[i],las = 2,ylog = F)
    }
    
    
  }
  
}

pdf('figures/par.compare.pdf',width = 4*2,height = 4*.618*3)


# labels
y.nm.vec <- c(expression(T[opt]),expression(r[extract]),
              expression(r[senescence]),expression(r[growth]),
              expression(q[growth]),expression(q[senescence]))

var.vec <- c(1,2,4,3,5,6)
for (spc.i in seq_along(spc.vec)) {
  # par(mfrow=c(2,1))
  par(mfrow=c(3,2))
  par(mar=c(5,5,3,1))
  
  fn <- sprintf('cache/%s.chain.%s.Control.Ambient.rds','smv1.2q',spc.vec[spc.i])
  fn.drt <- sprintf('cache/%s.chain.%s.Drought.Ambient.rds','smv1.drt',spc.vec[spc.i])
  plot.violin.func(fn.crtl = fn,
                   fn.drt = fn.drt,spc.nm = spc.vec[spc.i])
  par(new=T,mfrow=c(1,1))
  plot(0,ann=F,axes=F,xlab='',ylab='',pch=NA
      )
  
  title(spc.vec[spc.i])


}


dev.off()
