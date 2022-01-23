# #######################################
# read fitted par valyes and make table
# #######################################

# function to get quantiles
get.fit.ci.func <- function(fn,burin.frac=0.75){
  in.chain =  readRDS(fn)
  # burnIn = 1
  chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round(nrow(m.in)* (burin.frac)):nrow(m.in),])
  
  chain.fes <- do.call(rbind,chain.3.ls.new)
  hist(chain.fes[,3])
  out.df <- data.frame(f.t.opt = quantile(chain.fes[,1],probs = c(0.05,0.95,.5)),
                       f.extract = quantile(chain.fes[,2],probs = c(0.05,0.95,.5)),
                       f.sec = quantile(chain.fes[,3],probs = c(0.05,0.95,.5)),
                       f.growth= quantile(chain.fes[,4],probs = c(0.05,0.95,.5)),
                       q = quantile(chain.fes[,5],probs = c(0.05,0.95,.5)),
                       q.s = quantile(chain.fes[,6],probs = c(0.05,0.95,.5))
  )
  
  return(out.df)
}


# loop for all spc and pars
tmp.ls <- list()

spc.vec <- species.vec

for (spc.i in seq_along(spc.vec)) {
  fn <- sprintf('cache/smv1.2q.chain.%s.Control.Ambient.rds',spc.vec[spc.i])
  
  # v13.chain <- get.fit.value.func(fn)
  
  v13.chain.ci <- get.fit.ci.func(fn)
  
  tmp.ls[[spc.i]] <- data.frame(model='v1.1',
                                site = spc.vec[spc.i],
                                
                                f.t.opt = v13.chain.ci[3,1],
                                f.t.opt.05 = v13.chain.ci[1,1],
                                f.t.opt.95 = v13.chain.ci[2,1],
                                
                                f.extract = v13.chain.ci[3,2] ,
                                f.extract.05 = v13.chain.ci[1,2],
                                f.extract.95 = v13.chain.ci[2,2],
                                
                                
                                f.sec = v13.chain.ci[3,3] ,
                                f.sec.05 = v13.chain.ci[1,3],
                                f.sec.95 = v13.chain.ci[2,3],
                                
                                f.growth= v13.chain.ci[3,4],
                                f.growth.05 = v13.chain.ci[1,4],
                                f.growth.95 = v13.chain.ci[2,4],
                                
                                q = v13.chain.ci[3,5],
                                q.05 = v13.chain.ci[1,5],
                                q.95 = v13.chain.ci[2,5],
                                
                                q.s =v13.chain.ci[3,6],
                                q.s.05 = v13.chain.ci[1,6],
                                q.s.95 = v13.chain.ci[2,6]
                                
                                
  )
}

# save results
out.df = do.call(rbind,tmp.ls)

write.csv(out.df,'cache/fittedParValue.csv',row.names = F)


