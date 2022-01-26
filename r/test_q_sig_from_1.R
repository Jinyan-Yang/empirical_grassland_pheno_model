source('r/read_spc_nm.R')

# process ym data and pace data
pace.q.ls <- list()

for (i in seq_along(species.vec)) {
  fn.con.11 <- sprintf('cache/smv1.2q.chain.%s.Control.Ambient.rds',
                       species.vec[i])
  
  par.con.11 <- readRDS(fn.con.11)
  
  chain.3.ls.new = lapply(par.con.11,function(m.in)m.in[round((0.75)*nrow(m.in)):nrow(m.in),])
  
  chain.fes <- do.call(rbind,chain.3.ls.new)
  
  q.d.vec <- chain.fes[,5] - 1
  
  quants <- quantile(q.d.vec,probs = c(0.05,0.95))

  sig <- quants[[1]]*quants[[2]]
  
  if(sig<0){
    sig.val=0
  }else{
    sig.val=1
  }
  
  pace.q.ls[[i]] <-  data.frame(spc = species.vec[i],
                                q.sig = sig.val)
  
}

pace.q.df <- do.call(rbind,pace.q.ls)

