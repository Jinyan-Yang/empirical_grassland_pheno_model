# ###############################################################
# test whether the fitted par are different each other
# ###############################################################

# function to get the 5% and 95% quantiles of the difference between the par vals of two fit
get.quantile.func <- function(sm.nm.1,sm.nm.2,var.index){
  
  
  # read only the best fits that were selected
  fn.1 <- sprintf('cache/v13.2q.chain.%s.bestfit.rds',sm.nm.1)
  fn.2 <- sprintf('cache/v13.2q.chain.%s.bestfit.rds',sm.nm.2)
  
  # # # read chain
  in.chain.1 =  readRDS(fn.1)
  in.chain.2 =  readRDS(fn.2)
  
  # get difference
  diff.1.2 <- in.chain.1[,var.index] - in.chain.2[,var.index]
  # get quantile
  bond <- quantile(diff.1.2,probs = c(0.05,0.95))
  # output
  out.df <- data.frame(pars = paste0(sm.nm.1,'-',sm.nm.2),
                       quantile.05 = bond[[1]],
                       quantile.95 = bond[[2]])
  
  return(out.df)
}

# wrap function to do all species and params
get.sig.all.var.func <- function(x){
  spc.vec <- species.vec
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
  
  all.spc.df$sig[all.spc.df$product > 0 & all.spc.df$quantile.05 >0] <- 1
  all.spc.df$sig[all.spc.df$product > 0 & all.spc.df$quantile.05 <0] <- -1
  
  return(all.spc.df)
}

# do for all par and spc
all.var.ls <- lapply(1:6,get.sig.all.var.func)
saveRDS(all.var.ls,'cache/compare.var.rds')


