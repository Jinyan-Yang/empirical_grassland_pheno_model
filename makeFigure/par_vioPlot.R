# function to get mean
get.fit.value.func <- function(fn,burin.frac=0.75){
  in.chain =  readRDS(fn)
  burnIn = 1
  chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round(nrow(m.in)* (1-burin.frac)):nrow(m.in),])
  
  chain.fes <- do.call(rbind,chain.3.ls.new)
  
  return(colMeans(chain.fes[burnIn:nrow(chain.fes),]))
}
# function to get CI
get.fit.ci.func <- function(fn,burin.frac=0.75){
  in.chain =  readRDS(fn)
  burnIn = 1
  chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round(nrow(m.in)* (1-burin.frac)):nrow(m.in),])
  
  chain.fes <- do.call(rbind,chain.3.ls.new)
  
  out.df <- data.frame(f.t.opt = quantile(chain.fes[,1],probs = c(0.05,0.95)),
                       f.extract = quantile(chain.fes[,2],probs = c(0.05,0.95)),
                       f.sec = quantile(chain.fes[,3],probs = c(0.05,0.95)),
                       f.growth= quantile(chain.fes[,4],probs = c(0.05,0.95)),
                       q = quantile(chain.fes[,5],probs = c(0.05,0.95)),
                       q.s = quantile(chain.fes[,6],probs = c(0.05,0.95))
  )
  
  return(out.df)
}

# 
# loop through all params####
tmp.ls <- list()

spc.vec <-c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM','Flux')
# spc.vec <-c('Kan','YM','Flux')
for (spc.i in seq_along(spc.vec)) {
  fn <- sprintf('cache/smsmv13.2q.chain.%s.Control.Ambient.rds',spc.vec[spc.i])
  
  v13.chain <- get.fit.value.func(fn)
  
  v13.chain.ci <- get.fit.ci.func(fn)
  
  tmp.ls[[spc.i]] <- data.frame(model='v1.1',
                                site = spc.vec[spc.i],
                                
                                f.t.opt = v13.chain[1],
                                f.t.opt.05 = v13.chain.ci[1,1],
                                f.t.opt.95 = v13.chain.ci[2,1],
                                
                                f.extract = v13.chain[2] ,
                                f.extract.05 = v13.chain.ci[1,2],
                                f.extract.95 = v13.chain.ci[2,2],
                                
                                
                                f.sec = v13.chain[3] ,
                                f.sec.05 = v13.chain.ci[1,3],
                                f.sec.95 = v13.chain.ci[2,3],
                                
                                f.growth= v13.chain[4],
                                f.growth.05 = v13.chain.ci[1,4],
                                f.growth.95 = v13.chain.ci[2,4],
                                
                                q = v13.chain[5],
                                q.05 = v13.chain.ci[1,5],
                                q.95 = v13.chain.ci[2,5],
                                
                                q.s =v13.chain[6],
                                q.s.05 = v13.chain.ci[1,6],
                                q.s.95 = v13.chain.ci[2,6]
                                
                                
  )
}

out.df = do.call(rbind,tmp.ls)

# prepare significance data####
# read in significant data
all.var.ls <- readRDS('cache/compare.var.rds')

# loop through all pars
out.ls <- list()
for (par.i in seq_along(all.var.ls)) {
  # create empty df to store info
  spc.condition.df <- data.frame(spc = spc.vec,
                                 condition.1 = '',
                                 condition.2 = '',
                                 condition.3 = '',
                                 condition.4 = '',
                                 condition.5 = '',
                                 condition.6 = '',
                                 condition.7 = '',
                                 condition.8 = '',
                                 condition.9 = '')
  
  
  # subset to parameter
  tmp.df <- all.var.ls[[par.i]]
  
  tmp.m <- matrix(NA,ncol=11,nrow=10)
  rownames(tmp.m)<- spc.vec
  colnames(tmp.m) <- c(spc.vec,'sig')
  for(spc.nm in seq_along(spc.vec)[-length(spc.vec)]){
    
    # subset the 
    index.spc <- grep(paste0(spc.vec[spc.nm],'-'),tmp.df$pars)
    
    tmp.m[(spc.nm+1):10,spc.nm] <- tmp.df$sig[index.spc]
  }
  
  tmp.m[,11] <- ''
  
  # here's the logic is to find the common ones and assign the same letter
  i= 1
  letter.nm <- 1
  while (i<=10){
    
    index.ns <- which(tmp.m[,i]==0)
    
    # consider different situation
    if(length(index.ns) == 0){
      final.index <- NULL
    }
    
    if(length(index.ns) == 1){
      final.index <- c(i,index.ns)
    }
    
    if(length(index.ns) >1){
      tmp.ls <- list()
      for (tmp.i in seq_along(index.ns)){
        tmp.ls[[tmp.i]] <- intersect(which(tmp.m[,index.ns[tmp.i]]==0),index.ns)
        
        if(length(tmp.ls[[tmp.i]])>0){
          tmp.ls[[tmp.i]][length(tmp.ls[[tmp.i]])+1] <- index.ns[tmp.i]
        }
      }
      final.index <- unique(unlist(tmp.ls))
      final.index[length(final.index)+1] <- i
    }
    
    final.index <- unique(final.index)
    # assign significant letter
    if(length(final.index)>0){
      tmp.m[final.index,11] <- paste0(tmp.m[final.index,11],letters[letter.nm])
      letter.nm <- letter.nm+1
    }

    if(length(final.index)==0 & tmp.m[i,11]==''){
      tmp.m[i,11] <- letters[letter.nm]
      letter.nm <- letter.nm+1
    }
    
    i <- i+1
 
  }
  
  # 
  out.ls[[par.i]] <- as.data.frame(tmp.m)
}

# make plot####
library(vioplot)
devtools::source_url("https://github.com/Jinyan-Yang/colors/blob/master/R/col.R?raw=TRUE")
palette(c(col.df$iris))
col.nm.vec <- c(1,1,2,2,2,3,3,3,4,4)

plot.box.func <- function(spc.vec,col2plot,burin.frac=0.75,y.nm){
  
  tmp.ls <- list()
  for (i in seq_along(spc.vec)) {
    fn <- sprintf('cache/smsmv13.2q.chain.%s.Control.Ambient.rds',spc.vec[i])
    # fn <- sprintf('cache/smv13.qs1.chain.%s.Control.Ambient.rds',spc.vec[i])
    
    in.chain =  readRDS(fn)
    burnIn = 1
    chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round(nrow(m.in)* (1-burin.frac)):nrow(m.in),])
    
    chain.fes <- do.call(rbind,chain.3.ls.new)
    
    tmp.ls[[i]] <- data.frame(spc = spc.vec[i],
                              par.val = chain.fes[,col2plot])
    
  }
 
  plot.df <- do.call(rbind,tmp.ls)
  plot.df$spc[plot.df$spc =='flux'] <- 'Flux Tower'
  plot.df$spc <- factor(plot.df$spc,levels=spc.vec)
  vioplot(par.val~spc,plot.df,col=col.nm.vec,
          xlab='',ylab=y.nm)
  

  # sig.df <- out.ls[[col2plot]]
  
# add significance letter
    # for (par.i in 1:10) {
    #   # find the position 
    #   adj.val <- 0.095 * (par.i-2)+0.172
    #   # correct for uneven start and end 
    #   if(par.i==1){
    #     adj.val <- 0.07
    #   }
    #   if(par.i==2){
    #     adj.val <- 0.172
    #   }
    #   if(par.i==9){
    #     adj.val <- 1-0.172
    #   }
    #   if(par.i==10){
    #     adj.val <- 1-0.07
    #   }
    #   mtext(sig.df[par.i,11],side = 3,line = 1+par.i%%2,adj = adj.val,cex=0.8)
    # }

  
}

pdf('figures/par.vioplot.pdf',width = 4*2,height = 4*.618*3)
par(mfrow=c(3,2))
par(mar=c(5,5,1,1))

# var.nm.vec <- c('f.t.opt', 'f.extract',
#                 'f.sec','f.growth', 
#                 'q','q.s')

y.nm.vec <- c(expression('Opt. T '*(degree*C)),
              expression('Max. tran. rate '*(mm~d^-1)),
              expression('Senescence rate '*(cover~d^-1)),
              expression('Growth rate '*(cover~d^-1)),
              expression('q (growth)'),
              expression(q[s]~(Senescence)))

var.vec <- c(1,2,3,4,6,5)
for (plot.var.nm in seq_along(var.vec)) {
  
  
  plot.box.func(spc.vec,
                col2plot = var.vec[var.vec],
                y.nm = y.nm.vec[var.vec[var.vec]])
  
  legend('topleft',legend = sprintf('(%s)',letters[plot.var.nm]),
         bty='n')
  
}


dev.off()
# # plot significance.
pdf('figures/significance.pdf',width = 4*2,height = 4*3)
y.nm.vec <- c(('Opt. T '),
              ('Max. tran. rate '),
              ('Senescence rate '),
              ('Growth rate '),
              ('q[Growth]'),
              ('q[Senescence]'))

library(raster)
par(mfrow=c(3,2))
par(mar=c(5,6,1,1))
index.nm <- 1
for(plot.var.nm in c(1,2,3,4,6,5)){
  # subset for only the par needed
  plot.df <- out.ls[[plot.var.nm]]
  plot.m <- as.matrix(plot.df[2:10,1:9])
  rownames(plot.m)[rownames(plot.m) == 'Flux'] <- 'Flux Tower'
  x.nm <- colnames(plot.m)
  y.nm <- rownames(plot.m)
  # conver to a matrix
  plot.m <- matrix(as.numeric(plot.m),ncol=9)
  # plot as a raster
  plot(raster(plot.m),breaks=c(-1,0.5,1.1),col=c('grey','black'),legend=F,
       ann=F,axes=F,main = y.nm.vec[plot.var.nm])
  legend('topleft',legend =paste0('(',letters[index.nm],") " ,
                                  y.nm.vec[plot.var.nm]),
         bty='n')
  # add grid lines
  abline(h=(0:10)/9,lty='dotted')
  abline(v=(0:10)/9,lty='dotted')
  # 
  axis(1,at = (1:9)/9-0.05,labels = x.nm)
  axis(2,at = (1:9)/9-0.05,labels = rev(y.nm),las = 1)
  # 
  index.nm=index.nm+1
  
  
}

dev.off()

