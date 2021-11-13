
find.par.func <- function(par.num,set.value){
  
  tmp.ls <- list()
  for (i in seq_along(species.vec)) {
    fn <- sprintf('cache/smsmv13.2q.chain.%s.Control.Ambient.rds',species.vec[i])
    # fn <- sprintf('cache/smv13.qs1.chain.%s.Control.Ambient.rds',spc.vec[i])
    
    in.chain =  readRDS(fn)
    burnIn = 1
    chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round(nrow(m.in)* (1-0.75)):nrow(m.in),])
    
    chain.fes <- do.call(rbind,chain.3.ls.new)
    
    par.val = chain.fes[,par.num]
    
    new.chain <- (par.val-set.value)
    
    x <- quantile(new.chain,probs=c(0.05,0.95))
    
    if(x[[1]] * x[[2]]>0){
      sig = 1
    }else{
      sig = 0
    }
    
    tmp.ls[[i]] <- data.frame(spc = species.vec.nm[i],
                              par.num = par.num,
                              sig = sig)
    
  }
  
  q.df <- do.call(rbind,tmp.ls)
  # plot.df$spc[plot.df$spc =='flux'] <- 'Flux Tower'
  q.df$spc <- factor(q.df$spc,levels=unique(q.df$spc))
  return(q.df)
}

q.df <- find.par.func(5,set.value=1)
qs.df <- find.par.func(6,set.value=0)

