source('r/read_spc_nm.R')
source('r/load.R')

# # 
#daisgnostic plot
pdf('figures/diagV1.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec)) {

  fn <- sprintf('cache/smv1.2q.chain.%s.Control.Ambient.rds',species.vec[i])
  chain.3.ls = readRDS(fn)
  lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i])

  par(mfrow=c(3,2),mar=c(5,5,1,1))
  for(par.num in seq_along(colnames(chain.3.ls[[1]]))){

    start.row <- nrow(chain.3.ls[[1]]) / 4*3

    plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(start.row:nrow(chain.3.ls[[1]])),
                        nm.vec = colnames(chain.3.ls[[1]])[par.num])

  }
}
dev.off()

pdf('figures/diagV0.pdf',width = 8,height = 8*.618)
for (i in seq_along(species.vec)) {
  
  fn <- sprintf('cache/smv0.chain.%s.Control.Ambient.rds',species.vec[i])
  chain.3.ls = readRDS(fn)
  lapply(chain.3.ls, plot.check.mcmc.func,species.in=species.vec[i])
  
  par(mfrow=c(3,2),mar=c(5,5,1,1))
  for(par.num in 1:5){
    
    start.row <- nrow(chain.3.ls[[1]]) / 4*3
    
    plot.line.mcmc.func(chain.3.ls,par.num,range.iter =  round(start.row:nrow(chain.3.ls[[1]])))
    
  }
}

dev.off()