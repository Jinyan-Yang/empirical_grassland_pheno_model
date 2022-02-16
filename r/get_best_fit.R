# read fitted chains for each species and getthe best fit par############
#########################################################################

# function to extract par
get.par.func <- function(fn,iter){
  # read file
  in.chain =  readRDS(fn)
  # 
  chain.fes <- do.call(rbind,in.chain)
  chain.fes <- chain.fes[!duplicated(chain.fes),]
  best.1000 <- chain.fes[order(chain.fes$ll,decreasing=TRUE),]
  # select
  chain.sample <- best.1000[1:iter,]
  return(chain.sample)
}

# set the num of iteration to extract
sample.size <- 1000
# loop through all spc
for(spc.i in seq_along(species.vec)){
  species.in <- species.vec[spc.i]
  # get param for v`13`
  fn.v13=sprintf('cache/smv1.2q.chain.%s.Control.Ambient.rds',
            species.in)
  # 
  par.v13 <- get.par.func(fn=fn.v13,iter = sample.size)
  # 
  fn.v13.out <- sprintf('cache/v13.2q.chain.%s.bestfit.rds',
                        species.in)
  # save output
  saveRDS(par.v13,fn.v13.out)
  
  # get param for v`10`
  fn.v10=sprintf('cache/smv0.chain.%s.Control.Ambient.rds',
                 species.in)
  # 
  par.v10 <- get.par.func(fn=fn.v10,iter = sample.size)
  # 
  fn.v10.out <- sprintf('cache/v13.q1.qs0.chain.%s.bestfit.rds',
                        species.in)
  # save output
  saveRDS(par.v10,fn.v10.out)
}

