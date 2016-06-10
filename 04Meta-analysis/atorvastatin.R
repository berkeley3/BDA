#--- Load libraries ---#
library(metafor)	# This is the meta-analysis library in R
library(rjags)



# 1. Read the data, derive the mean difference from placebo, and 
# corresponding standard error. 

d <- read.csv('material\\AtorvastatinStudies.csv')
ator <- subset(d, dose1>0)
pbo <- subset(d, dose1==0, select=c("trial","n","ldlPcfb","ldl0",
                                    "seLdl0","seLdlPcfb"))
names(pbo) <- paste(names(pbo),".pbo",sep='')

ator <- merge(ator, pbo, by.x="trial",by.y="trial.pbo", all=T)

ator$diff <- ator$ldlPcfb - ator$ldlPcfb.pbo


#--- Calculate the difference and associated standard error ---#
ator$seDiff <- sqrt(ator$seLdlPcfb^2 + ator$seLdlPcfb.pbo^2)

ator10 <- ator[ator$dose1==10,]

# 2. Make a forest plot of the data.

# Table of data 
ator10 = ator10[(order(ator10$diff)),]

with(ator10,
     forest(x=diff, sei=seDiff, 
            slab=trial, refline=NA, 
            at=seq(-80,-10,by=10),
            ylim=c(0,11),
            #            xlim=c(-120,40),
            xlab="Mean Difference",
            main="Difference in %CFB LDL for 10 mg Atorvastatin")
)
text(-90.5, 9.5, "Study", cex=0.95, adj=0, font=2)
text(-15, 9.5, "Mean Diff. + 95%CI",cex=0.95, adj=0, font=2)

# 3. Fit a random effects model 

cat("
    model{
    for (i in 1:nStudies) {
    thetahat[i] ~ dnorm(theta[i], precThetaHat[i])
    theta[i] ~ dnorm(meanDiff, prec)
    precThetaHat[i] <- pow(s[i],-2)
    }
    meanDiff ~ dnorm(0,1E-04)
    tau ~ dunif(0,100)
    prec <- pow(tau,-2)
    }", file="atorvastatin.txt")

bugsdata = with(ator10, 
                list(thetahat=diff, s=seDiff, 
                     nStudies=length(trial))
)

parameters = c("meanDiff","tau","theta")

inits = list(list(meanDiff=rnorm(1, -20, 20), tau=runif(1,0,5)), 
             list(meanDiff=rnorm(1, -20, 20), tau=runif(1,0,5)),
             list(meanDiff=rnorm(1, -20, 20), tau=runif(1,0,5)))

parameters =  c("meanDiff","tau","theta")     # The parameter(s) to be monitored.
adaptSteps = 500              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 3                  # Number of chains to run.
numSavedSteps=50000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.


jagsModel = jags.model( "atorvastatin.txt" , data=bugsdata ,  inits=inits ,
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

with(ator10,
     forest(x=diff, sei=seDiff, 
            slab=trial, refline=NA, 
            at=seq(-80,-10,by=10),
            ylim=c(0,11),
            #            xlim=c(-120,40),
            xlab="Mean Difference",
            main="Difference in %CFB LDL for 10 mg Atorvastatin")
)
text(-90.5, 9.5, "Study", cex=0.95, adj=0, font=2)
text(-15, 9.5, "Mean Diff. + 95%CI",cex=0.95, adj=0, font=2)

x1 <- samples[[1]]
x2 <- samples[[2]]
x3 <- samples[[3]]


addpoly(x=mean(c(x1[,1], x2[,1],x3[,1])), 
        vi=var(c(x1[,1], x2[,1],x3[,1])),
        row=0, col='red', mlab='Bayesian RE model')
