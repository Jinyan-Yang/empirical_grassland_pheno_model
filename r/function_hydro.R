##### 
# calculate soil psi
psi.s.func <- function(swc,
                       psi.e = -0.74e-3,#KPa
                       b = 6.4, 
                       swc.sat = 0.24){
  psi.e * (swc / swc.sat)^-b
}


# calculate soil conductivity
k.soil.func <- function(swc,
                        swc.sat = 0.24,
                        psi.e = -0.74e-3,#KPa
                        b = 6.4, 
                        k.sat = 11.4){
  
  psi.soil <- psi.s.func(swc = swc,
                         psi.e = psi.e,#KPa
                         b = b, 
                         swc.sat = swc.sat)
  
  k.s <- k.sat * (psi.e / psi.soil)^(2+3/b)
  
  # set max to k.sat
  k.s[which(k.s> k.sat)] <- k.sat
  return(k.s)
}

# calculate soil conductance
k.s.p.func <- function(swc,
                       r.l=0.48,
                       lai=8,
                       l.v=1490,
                       r.root=0.0005,
                       swc.sat = 0.13,
                       psi.e = -0.36e-3,#KPa
                       b = 4.26, 
                       k.sat = 79.8){
  ks <- r.l*2*pi*
    k.soil.func(swc,swc.sat,psi.e,b,k.sat)/
    lai/log10(1/sqrt(pi * l.v)/r.root)
  
  return(ks)
}

# calculate Emax 
e.frac.func <- function(swc,
                        k.plant = 0.02,#plant hydro conductance; mmol m–2 s–1 MPa–1
                        swc.sat = 0.13,
                        psi.e = -0.74e-3,#KPa
                        b = 4.26, 
                        k.sat = 79.8,#soil sat hydo conductivity; molm–1 s–1MPa–1
                        psi.min = -2,#min plant psi; Mpa
                        r.l = 0.48,#Root length index; m m–2
                        lai = 8,
                        l.v = 1490,#Root length density; m m–3
                        r.root = 0.0005#Mean radius of water absorbing roots; m 
){
  # get psi.soil
  psi.soil <- psi.s.func(swc = swc,
                         psi.e = psi.e,#KPa
                         b = b, 
                         swc.sat = swc.sat)
  
  # psi.soil <- max(psi.soil, psi.min)
  
  # get resistance
  R.p = 1 / k.plant
  R.s = 1 / k.s.p.func(swc,swc.sat=swc.sat,psi.e=psi.e,b=b,k.sat=k.sat,
                       r.l=r.l,
                       lai=lai,
                       l.v=l.v,
                       r.root=r.root)
  # get conductance
  k.l <- 1 / (R.p + R.s)
  # get max E under swc
  e.max <- k.l * (psi.soil - psi.min)
  e.max <- max(e.max,0)
  # get resistance under sat swc
  R.s.sat <- 1 / k.s.p.func(swc.sat,swc.sat=swc.sat,psi.e=psi.e,b=b,k.sat=k.sat,
                            r.l=r.l,
                            lai=lai,
                            l.v=l.v,
                            r.root=r.root)
  
  k.l.sat <- 1/(R.p + R.s.sat)
  # get max e under swc.sat
  e.max.sat <- k.l.sat * (psi.e - psi.min)
  
  return(e.max/e.max.sat)
}

