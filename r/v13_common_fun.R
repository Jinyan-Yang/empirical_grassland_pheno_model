###################################################################
# functions to do the mcmc fitting fot V13 version of the model####
###################################################################

# source('r/functions_mcmc_v12.r')

# wrpped function to fit mcmc to a df
fit.mcmc.2q.func <- function(df = gcc.met.pace.df,
                             species.in = 'Luc',prep.in = 'Control', 
                             temp.in ='Ambient',subplot =NA,
                             my.fun = phenoGrass.func.v13,
                             out.nm.note = '',use.smooth = FALSE,
                             day.lag = 3,
                             bucket.size = 300,
                             swc.capacity = 0.13,
                             swc.wilt = 0.05,
                             n.iter = 10000,
                             norm.min.max=NULL,
                             cal.initial=F,
                             par.df,q.given =NULL,q.s.given=NULL,
                             use.mcmc=TRUE,
                             de.note=''){
  
  # check for error in input
  if(ncol(par.df)==6 & !is.null(q.s.given)){
    stop('qs not set correctly. ')
  }

  if(ncol(par.df)==5 & !is.null(q.given)){
    stop('q not set correctly. ')
  }
  
  # 
  s.time <- Sys.time()
  gcc.met.pace.df.16 <<- get.pace.func(df,
                                      species.in =species.in,
                                      prep.in = prep.in,
                                      temp.in =temp.in,
                                      subplot = subplot,
                                      norm.min.max = norm.min.max)
  
  my.fun <<-  my.fun
  use.smooth <<-  use.smooth
  day.lag <<- day.lag
  bucket.size <<-  bucket.size
  swc.capacity <<- swc.capacity
  swc.wilt <<-  swc.wilt
  n.iter <<-  n.iter
  
  # get deoptim to fit initial values
  if(cal.initial){
  
    initial.vec <- get.ini.func(par.df = par.df,q.given =q.given,q.s.given=q.s.given)
     # out.nm <- paste0('tmp/deopt_',q.given,q.s.given,de.note,species.in,prep.in,temp.in,'.rds')
     # saveRDS(initial.vec,out.nm)
    # initial.vec <-c(19.4,1,0.198,0.1)
    }else{
    initial.vec<-NULL
  }
# decide if mcmc should be use
  if(use.mcmc){
    # prepare param
    if(is.null(initial.vec)){
      par.df['initial',] <- c(20,1.5,0.005,0.15,3,0.5)[seq_along(par.df['initial',])]
    }else{
      par.df['initial',] <- initial.vec
    }
    par.df['stdv',] <- (par.df['max',] - par.df['min',])/100
    
    par.df <<- par.df
    
    #setup parallel backend to use many processors
    cores=3#detectCores(logical = FALSE)
    cl <- makeCluster(cores[1]) #not to overload your computer
    registerDoParallel(cl)
    
    on.exit(stopCluster(cl))
    
    chain.fes <- foreach(i=1:3,
                         .packages = 'mvtnorm',
                         .export=ls(envir=globalenv())) %dopar% {
                           
                           chain.tmp = mh.MCMC.func.2q(iterations=n.iter,
                                                       par.df = par.df,
                                                       gcc.met.pace.df.16 = gcc.met.pace.df.16,
                                                       bucket.size = bucket.size,
                                                       day.lay = day.lag,
                                                       swc.capacity = swc.capacity,
                                                       swc.wilt = swc.wilt,
                                                       my.fun = phenoGrass.func.v13,
                                                       use.smooth = T,q.given =q.given,q.s.given=q.s.given)
                         }
    
    
    
    # save file
    if(use.smooth==TRUE){
      smooth.nm='sm'
    }else{
      smooth.nm=''
    }
    
    if(is.na(subplot)){
      out.name <- sprintf('cache/%s%schain.%s.%s.%s.rds',smooth.nm,out.nm.note,species.in,prep.in,temp.in)
    }else{
      out.name <- sprintf('cache/%schain.%s.rds',out.nm.note,subplot)
    }
    
    saveRDS(chain.fes,out.name)
    
    time.used <- (Sys.time() - s.time)
    
    df <- data.frame(site = species.in,
                     time = time.used,
                     when = Sys.time())
    
    write.table(df, file = "time.used.txt", sep = "\t",append = T,
                row.names = FALSE)
  }
  
}

# generic function to do mcmc in paralle
mh.MCMC.func.2q <- function(iterations,par.df,
                            gcc.met.pace.df.16,
                            bucket.size,swc.wilt ,
                            swc.capacity,day.lay,
                            my.fun,
                            use.smooth,q.given,q.s.given){


  
  # get prior 
  prior.prob <- prior.func(par.df)
  
  # start MC ####
  # intial
  chain = array(dim = c(iterations+1,ncol(par.df)))
  # chain[1,] = as.numeric(par.df['initial',])
  chain <- as.data.frame(chain)
  names(chain) <- names(par.df)
  chain[1,] = (par.df['initial',])
  
  # chain move on
  ll.vec <- c()
 
  for (i in 1:iterations){
    # prpose a set of par values based on previous chain value
    proposal = proposal.func(chain[i,],par.df)
  
    # deal with missing q 
    if(is.null(q.given)){
      q.val = proposal[5]
      q.past <- chain$q[i]
    }else{
      q.past = q.val = q.given
    }
    
    if(is.null(q.s.given)){
      q.s.val = proposal[6]
      q.s.past <- chain$q.s[i]
    }else{
      q.s.past = q.s.val = q.s.given
    }

    # prior.prob,data,data.sd,bucket.size = 300,...
    
    proposal.ll <- posterior.func(prior.prob,FUN = my.fun,
                                  gcc.df = gcc.met.pace.df.16,
                                  f.h = 222,
                                  f.t.opt = proposal[1],
                                  f.extract = proposal[2],
                                  f.sec = proposal[3],
                                  f.growth = proposal[4] ,
                                  q = q.val ,
                                  q.s = q.s.val ,
                                  t.max = 45,
                                  day.lay = day.lay,
                                  bucket.size = bucket.size,
                                  swc.wilt = swc.wilt ,
                                  swc.capacity = swc.capacity,
                                  use.smooth = use.smooth)
    former.ll <- posterior.func(prior.prob,FUN = my.fun,
                                gcc.df = gcc.met.pace.df.16,
                                f.h = 222,
                                f.t.opt = chain[i,1],
                                f.extract = chain[i,2],
                                f.sec = chain[i,3],
                                f.growth = chain[i,4] ,
                                q =  q.past,
                                q.s = q.s.past ,
                                t.max = 45,
                                day.lay = day.lay,
                                swc.wilt = swc.wilt ,
                                swc.capacity = swc.capacity,
                                bucket.size = bucket.size,
                                use.smooth = use.smooth)
    # save loglikelhood
    if(i==1){
      ll.vec[1] <- former.ll
    }
    
    probab = exp(proposal.ll- former.ll)
    if (runif(1) < probab){
      chain[i+1,] = proposal
      ll.vec[i+1] <- proposal.ll
    }else{
      ll.vec[i+1] <- former.ll
      chain[i+1,] = chain[i,]
    }
  }
  
  chain$ll <- ll.vec
  
  return(chain)
}
