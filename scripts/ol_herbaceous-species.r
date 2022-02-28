library(jagsUI)

data <- read.csv("ocup_ol.csv", header = T)
cova <- read.csv("covas_herb.csv", header = T)

str(win.data <- list(y = y,rr = rr,ac = ac,na = na,sh = sh,ms = ms,
                     nsite=  dim(y)[1],nrep = dim(y)[2],nspec= dim(y)[3],nyear= dim(y)[4]
                     
))
#model BUGS
sink("ol_herb.txt")
cat("model {
    # model for missing covariates
    for (t in 1:nyear){ 
      for(j in 1:nsite){
    
     as[j,t] ~ dnorm(mu.as, tau.as)
     hl[j,t] ~ dnorm(mu.hl, tau.hl)
     pp[j,t] ~ dnorm(mu.pp, tau.pp)
     pl[j,t] ~ dnorm(mu.pl, tau.pl)
     ra[j,t] ~ dnorm(mu.ra, tau.ra)
    
    }
    }
    tau.as <-pow(sd.as,-2)
    sd.as ~ dunif(0,1)
    mu.as ~ dnorm(0,1)
    
    tau.hl <-pow(sd.hl,-2)
    sd.hl ~ dunif(0,1)
    mu.hl ~ dnorm(0,1)
    
    tau.pp <-pow(sd.pp,-2)
    sd.pp ~ dunif(0,1)
    mu.pp ~ dnorm(0,1)
    
    tau.pl <-pow(sd.pl,-2)
    sd.pl ~ dunif(0,1)
    mu.pl ~ dnorm(0,1)
    
    tau.ra <-pow(sd.ra,-2)
    sd.ra ~ dunif(0,1)
    mu.ra ~ dnorm(0,1)
    

    mu.lpsi    ~ dnorm(0, 1/1.25^2)    
    mu.lp      ~ dnorm(0, 1/1.25^2)
    betalpsi1  ~ dnorm(0, 1/1.25^2)
    betalpsi2  ~ dnorm(0, 1/1.25^2)
    betalpsi3  ~ dnorm(0, 1/1.25^2)
    betalpsi4  ~ dnorm(0, 1/1.25^2)
    betalpsi5  ~ dnorm(0, 1/1.25^2)
    
    
    tau.lpsi <- pow(sd.lpsi, -2)   
    tau.lp   <- pow(sd.lp, -2) 
    sd.lpsi  ~ dunif(0,1)  
    sd.lp    ~ dunif(0,1) 
    
 # priors for species specific effects in occupancy and detection 
  
    for (t in 1:nyear){
      lpsi[t] ~ dnorm(mu.lpsi, tau.lpsi) 
      lp  [t] ~ dnorm(mu.lp,   tau.lp) 
    
    }
    
    # priors for species-specific and random year effects
    for (i in 1:nspec){
     for (t in 1:nyear){
     
      lpsi  [i,t] ~ dnorm (mu.sp.lpsi [t], tau.sp.lpsi [t])
      
    } 
    } 
        
    # Ecological model, process model (true occurrence at site j) 
    ## year as random effect
    for (t in 1:nyear){                                     
     for (j in 1:nsite){                                         
    
        logit(psi[j,t]) <- lpsi[t] + betalpsi1 * as[j,t]
                                   + betalpsi2 * hl[j,t]
                                   + betalpsi3 * pp[j,t]
                                   + betalpsi4 * pl[j,t]
                                   + betalpsi5 * ra[j,t]
    
    z[j,i,e] ~ dbern(psi[j,i,e])
    
    }
    }
    }
    
    # Observation model for site j, replicate nrep=k, nspec=i
    for (t in 1:nyear) {
    logit(p[t]) <- lp[t]      
    
    for (t in 1:nyear) {
     for (j in 1:nsite)  {           
      for (k in 1:nrep)   { 
      
       mu.p[j,k,t] <- z[j,e]*p[t]  
       y[j,k,t] ~ dbern(mu.p[j,k,t])
    
    }
    }
    } 
    }
    
    # Derived quantities
    
    N <- sum(z)
    } #model
    
    ",fill=TRUE)
sink()
zst <- apply(y,c(1,3,4),max)
inits <- function(){list(z=zst)}

params1 <- c("mu.sp.lpsi","lpsi","lp","betalpsi1","betalpsi2","betalpsi3","betalpsi4","betalpsi5")

ni <- 500000    
nt <- 10  
nb <- 250000
nc <- 3

out1 = jags(win.data, inits, params1, "ol_herb.txt", n.chains=nc, 
             n.iter=ni, n.burnin=nb, n.thin=nt)

