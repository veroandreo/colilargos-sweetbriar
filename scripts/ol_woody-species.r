library(jagsUI)

data <- read.csv("ocup_ol.csv", header = T)
cova <- read.csv("covas_woody.csv", header = T)

str(win.data <- list(y = y,rr = rr,ac = ac,na = na,sh = sh,ms = ms,
                     nsite=  dim(y)[1],nrep = dim(y)[2],nspec= dim(y)[3],nyear= dim(y)[4]
                     
))
#model BUGS
sink("ol_woody.txt")
cat("model {
    # model for missing covariates
    for (t in 1:nyear){ 
      for(j in 1:nsite){
    
     rr[j,t] ~ dnorm(mu.rr, tau.rr)
     ac[j,t] ~ dnorm(mu.ac, tau.ac)
     na[j,t] ~ dnorm(mu.na, tau.na)
     ot[j,t] ~ dnorm(mu.ot, tau.ot)
     ms[j,t] ~ dnorm(mu.ms, tau.ms)
    
    }
    }
    tau.rr <-pow(sd.r,-2)
    sd.rr ~ dunif(0,1)
    mu.rr ~ dnorm(0,1)
    
    tau.ac <-pow(sd.ac,-2)
    sd.ac ~ dunif(0,1)
    mu.ac ~ dnorm(0,1)
    
    tau.na <-pow(sd.na,-2)
    sd.na ~ dunif(0,1)
    mu.na ~ dnorm(0,1)
    
    tau.ot <-pow(sd.ot,-2)
    sd.ot ~ dunif(0,1)
    mu.ot ~ dnorm(0,1)
    
    tau.ms <-pow(sd.ms,-2)
    sd.ms ~ dunif(0,1)
    mu.ms ~ dnorm(0,1)
    

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
    
        logit(psi[j,t]) <- lpsi[t] + betalpsi1 * rr[j,t]
                                   + betalpsi2 * ac[j,t]
                                   + betalpsi3 * na[j,t]
                                   + betalpsi4 * ot[j,t]
                                   + betalpsi5 * ms[j,t]
    
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

out1 = jags(win.data, inits, params1, "ol_woody.txt", n.chains=nc, 
             n.iter=ni, n.burnin=nb, n.thin=nt)

