source('r/read_spc_nm.R')
# function to get the 5% and 95% quantiles of the difference between the par vals of two fit
get.quantile.func <- function(sm.nm.1,sm.nm.2,var.index){
  # read chain
  fn.1=paste0('cache/smsmv13.2q.chain.',sm.nm.1,'.Control.Ambient.rds')
  fn.2=paste0('cache/smsmv13.2q.chain.',sm.nm.2,'.Control.Ambient.rds')
  # set burn in
  burin.frac=0.75
  # combine chains
  in.chain.1 =  readRDS(fn.1)
  chain.1= lapply(in.chain.1,function(m.in)m.in[round(nrow(m.in)* (1-burin.frac)):nrow(m.in),])
  chain.1.full <- do.call(rbind,chain.1)
  
  in.chain.2 =  readRDS(fn.2)
  chain.2= lapply(in.chain.2,function(m.in)m.in[round(nrow(m.in)* (1-burin.frac)):nrow(m.in),])
  chain.2.full <- do.call(rbind,chain.2)
  
  # reorder chains
  set.seed(1965)
  chain.1.random <- sample(chain.1.full[,var.index])
  set.seed(1965)
  chain.2.random <- sample(chain.2.full[,var.index])
  
  # get difference
  diff.1.2 <- chain.1.random - chain.2.random
  # get quantile
  bond <- quantile(diff.1.2,probs = c(0.05,0.95))
  # output
  out.df <- data.frame(pars = paste0(sm.nm.1,'-',sm.nm.2),
                       quantile.05 = bond[[1]],
                       quantile.95 = bond[[2]])
  
  return(out.df)
}

# function to do all species and params
get.sig.all.var.func <- function(x){
  spc.vec <- species.vec#c('Bis','Luc','Dig','Kan','Rho','Fes','Pha','Rye','YM','Flux')
  nm.var.tot<-  length(spc.vec)
  
  out.ls <- list()
  for (spc.1.nm in seq_along(spc.vec)){
    
    tmp.spc.ls <- list()
    for (spc.2.nm in (spc.1.nm):nm.var.tot) {
      if(spc.2.nm != spc.1.nm){
        sm.nm.1 <- spc.vec[spc.1.nm]
        sm.nm.2 <- spc.vec[spc.2.nm]
        tmp.df <- get.quantile.func(sm.nm.1,sm.nm.2,x)
        tmp.spc.ls[[spc.2.nm]] <- tmp.df
      }
      
    }
    tmp.spc.df <- do.call(rbind,tmp.spc.ls)
    out.ls[[spc.1.nm]] <- tmp.spc.df
  }
  
  all.spc.df <- do.call(rbind,out.ls)
  
  
  all.spc.df$sig <- 0
  all.spc.df$product <- all.spc.df$quantile.05 * all.spc.df$quantile.95
  
  all.spc.df$sig[all.spc.df$product > 0 ] <- 1
  return(all.spc.df)
}

all.var.ls <- lapply(1:6,get.sig.all.var.func)
saveRDS(all.var.ls,'cache/compare.var.rds')
# index.spc <- grep(paste0(spc.vec[1],'-'),all.spc.df$pars)
# df.sub <- all.spc.df[index.spc,]

