
# 1. Read the data, derive the log odds ratio 
# comparing GC to G and corresponding standard error.

d <- read.csv('material\\CochranePancreaticData12months.csv')
d$logOR <- with(d, log(d.GC*(N.G-d.G)/(d.G*(N.GC-d.GC))))
d$varLogOR <- with(d, 1/d.GC + 1/(N.GC-d.GC) + 1/d.G + 1/(N.G-d.G))


# 2. Make a forest plot of the log odds ratio values. 

d <- d[order(d$logOR),]
with(d,
     forest(x=logOR, vi=varLogOR, 
            slab=study, 
            atransf=exp,
            at=log(c(0.1,0.25,0.5,1,4,16)),
            alim=log(c(0.03,24)),
            refline=0, 
            xlab="Odds Ratio",
            main="Odds Ratio for 12 month mortality\nGemcitabine+Chemo to Chemo")
)
text(log(0.0005), nrow(d)+1.5, "Study", cex=1, adj=0, font=2)
text(log(20), nrow(d)+1.5, "Odds Ratio [95%CI]",cex=1, adj=0, font=2)




# 3. Fit a random effect model

cat("model{
  for (i in 1:nStudies) {
    thetahat[i] ~ dnorm(theta[i], precThetaHat[i])
    theta[i] ~ dnorm(meanDiff, prec)
    precThetaHat[i] <- pow(s[i],-2)
  }
  meanDiff ~ dnorm(0,1E-04)
  tau ~ dunif(0,100)
  prec <- pow(tau,-2)
}", file='pancreatic.txt')

bugsdata <- list(thetahat=d$logOR, s=sqrt(d$varLogOR), 
                 nStudies=length(d$study))




parameters = c("meanDiff","tau","theta")

inits = list(list(meanDiff=rnorm(1, -5, 5), tau=runif(1,0,5)), 
             list(meanDiff=rnorm(1, -5, 5), tau=runif(1,0,5)),
             list(meanDiff=rnorm(1, -5, 5), tau=runif(1,0,5)))

parameters =  c("meanDiff","tau","theta")     # The parameter(s) to be monitored.
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 3                  # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.


jagsModel = jags.model( "pancreatic.txt" , data=bugsdata ,  inits=inits ,
                        n.chains=nChains , n.adapt=adaptSteps )

cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
samples = coda.samples( jagsModel , variable.names=parameters ,
                        n.iter=nIter , thin=thinSteps )

# look at convergence diagnostics
plot(samples)
summary(samples, start = burnInSteps)

plot(jagsModel)
gelman.diag(samples)
gelman.plot(samples)





# Add Bayes estimate to plot
with(d,
     forest(x=logOR, vi=varLogOR, 
            slab=study, 
            atransf=exp,
            at=log(c(0.1,0.25,0.5,1,4,16)),
            alim=log(c(0.03,24)),
            refline=0, 
            xlab="Odds Ratio",
            main="Odds Ratio for 12 month mortality\nGemcitabine+Chemo to Chemo")
)
text(log(0.0005), nrow(d)+1.5, "Study", cex=1, adj=0, font=2)
text(log(20), nrow(d)+1.5, "Odds Ratio [95%CI]",cex=1, adj=0, font=2)

# Adding additional estimates

x1 <- samples[[1]]
x2 <- samples[[2]]
x3 <- samples[[3]]

addpoly(x=mean(c(x1[,1], x2[,1],x3[,1])), 
        vi=var(c(x1[,1], x2[,1],x3[,1])),
        row=0, col='red', mlab='Bayesian RE model')
addpoly(x=mean(reBugs1$sims.list$meanDiff), 
        vi=var(reBugs1$sims.list$meanDiff),
        row=-2, col='red', mlab='Bayesian RE model',
        atransf=exp, cex=1)
