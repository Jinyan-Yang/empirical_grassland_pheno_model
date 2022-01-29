###################################################################
# calculate CI of the fitted model by sampling from the chains####
###################################################################
# 

get.mod.ci.func <-  function(df = gcc.met.pace.df,
                             species.in,prep.in,temp.in,subplot=NULL,
                             nm.note='',use.smooth=FALSE,
                             my.fun = phenoGrass.func.v11,
                             swc.in.cap = 0.13,swc.in.wilt = 0.05,bucket.size =300,
                             norm.min.max=NULL,
                             day.lag=3,
                             do.predict =NULL,
                             burn.proportion = 0.75,q.s.in=NULL,q.in=NULL,sample.size =100){
  # data processing
  if(is.null(subplot)){
    gcc.met.pace.df.16 <- get.pace.func(df,
                                        species.in = species.in,
                                        prep.in = prep.in,
                                        temp.in =temp.in,
                                        norm.min.max=norm.min.max)
    
    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    if(is.null(do.predict)){
      fn=paste0('cache/',sm.nm,nm.note,'chain.',
                species.in,'.',prep.in,'.',temp.in,'.rds')
      rds.nm = paste0('tmp/ci.',sm.nm,nm.note,'chain.',species.in,'.',prep.in,'.',temp.in,'.rds')
    }else{
      fn=paste0('cache/',sm.nm,nm.note,'chain.',
                species.in,'.',do.predict,'.',temp.in,'.rds')
      rds.nm = paste0('tmp/ci.',sm.nm,nm.note,'chain.',species.in,'.',do.predict,'.predict.',temp.in,'.rds')
    }
    
    
    # fn=paste0('cache/chain.',species.in,'.',prep.in,'.',temp.in,'.rds')
    
    
  }else{
    species.in = subplot
    prep.in = ''
    temp.in =''
    gcc.met.pace.df.16 <- get.pace.func(gcc.met.pace.df,subplot = subplot)
    
    if(use.smooth){
      sm.nm='sm'
    }else{
      sm.nm=''
    }
    
    fn=paste0('cache/',sm.nm,nm.note,'chain.',subplot,'.rds')
    
    rds.nm = paste0('tmp/ci.',sm.nm,nm.note,'chain.',subplot,'.rds')
  }
  
  
  # read chains 
  # read only the best fits that were selected
  tmp.str <- gsub('sm',replacement = '',nm.note)
  fn <- sprintf('cache/%schain.%s.bestfit.rds',tmp.str,species.in)
  print(paste0('par file used: ',fn))
  chain.sample <- readRDS(fn)
  chain.sample <-  subset(chain.sample,select=-c(ll))
  # in.chain =  readRDS(fn)
  # 
  # # if(is.list(in.chain)){
  # #   # assuming 1/3 burn in
  # #   burnIn = 1
  # #   chain.3.ls.new = lapply(in.chain,function(m.in)m.in[round((burn.proportion)*nrow(m.in)):nrow(m.in),])
  # #   
  # #   chain.fes <- do.call(rbind,chain.3.ls.new)
  # # }else{
  # #   burnIn = nrow(in.chain)/3
  # #   chain.fes <- in.chain
  # # }
  # # 
  # # # chain.sub <- chain.fes[burnIn:nrow(chain.fes),]
  # # # 
  # # # set.seed(1935)
  # 
  # chain.fes <- do.call(rbind,in.chain)
  # best.1000 <- chain.fes[order(chain.fes$ll,decreasing=TRUE),]
  # # # sample chain
  # # len.chain <- nrow(chain.sub)
  # # sample.index <- sample(1:len.chain,sample.size)
  # 
  # # chain.sample <- as.data.frame(chain.sub[sample.index,])
  #   chain.sample <- best.1000[1:sample.size,]
  # deal with no q
  if(!('q' %in% names(chain.sample))){
    chain.sample$q <- q.in
    print(paste0('sensitivities of growth and senesence set to ',
                 q.in))
  }
  if(!('q.s' %in% names(chain.sample))){
    chain.sample$q.s <- q.s.in
    print(paste0('sensitivitiy of senesence set to ',q.s.in))
  }
  # make prediction with ci#####
  gcc.met.pace.df.16 <- gcc.met.pace.df.16[order(gcc.met.pace.df.16$Date),]
  pred.vec <- list()

  # bic.df <- data.frame(f.t.opt = chain.sample[,1],
  #                      f.extract = chain.sample[,2],
  #                      f.sec= chain.sample[,3],
  #                      f.growth = chain.sample[,4],
  #                      q =  chain.sample[,5],
  #                      q.s =  chain.sample[,6],
  #                      bic=NA)

  # bic.ls <- list()

  for(i in 1:sample.size){
    fit.par.vec <- chain.sample[i,]#unlist(unname(chain.sample[i,]))

    hufken.pace.pred <- my.fun(gcc.met.pace.df.16,
                               f.h = 222,
                               
                               f.t.opt = fit.par.vec$f.t.opt,
                               f.extract = fit.par.vec$f.extract,
                               f.sec= fit.par.vec$f.sec,
                               f.growth = fit.par.vec$f.growth,
                               q =  fit.par.vec$q,
                               q.s =  fit.par.vec$q.s,
                               
                               bucket.size = bucket.size,
                               swc.wilt = swc.in.wilt ,
                               swc.capacity = swc.in.cap ,
                               t.max = 45,
                               day.lay = day.lag,use.smooth = use.smooth)

    pred.vec[[i]] <- hufken.pace.pred$cover.hufken

    # bic.ls[[i]] <- data.frame(f.t.opt = fit.par.vec[1],
    #                           f.extract = fit.par.vec[2],
    #                           f.sec= fit.par.vec[3],
    #                           f.growth = fit.par.vec[4],
    #                           q =  fit.par.vec[5],
    #                           q.s =  fit.par.vec[6],
    #                           bic = get.bic.func(model.vec = hufken.pace.pred$cover.hufken,
    #                                               data.vec = hufken.pace.pred$cover,n.fd = 6))


  }

  # save results
  pred.m <- do.call(rbind,pred.vec)
  # bic.df <- do.call(rbind,bic.ls)
  # quantile(pred.m[1,],probs = c(0.05,0.95),na.rm=T)
  pred.ci.m <- apply(pred.m,2,quantile,probs = c(0.05,0.95,0.5),na.rm=T)

  # save prediction for future use
  saveRDS(pred.ci.m,rds.nm)
  # rds.nm.bic <- paste0('tmp/bic.',sm.nm,nm.note,'chain.',species.in,'.',prep.in,'.',temp.in,'.rds')
  # saveRDS(bic.df,rds.nm.bic)
}