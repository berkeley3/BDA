library(maptools)
library(sp)
library(spdep)
library(R2OpenBUGS)
library(RColorBrewer)
library(classInt)

# imort the map------------
load('data/ITA_adm3.RData') #url: http://gadm.org/country

vcmap=gadm[gadm$NAME_2=='Vercelli',]
#plot(vcmap)

#import db with observed and expected cases----
db=read.csv('data\\db.csv',sep=';')

## workout to merge municipalities  on your db with municipalities on the map object---
db$comune=tolower(db$COMUNE)

comuni=vcmap[,'NAME_3']@data
comuni[,1]=tolower(as.character(comuni[,1]))
a=setdiff(sort(tolower(as.character(comuni[,1]))),db$comune)
comuni[match(a,comuni[,1]),1]="borgo d'ale"
db.comuni=merge(comuni, db, by.y='comune', by.x='NAME_3', sort=F)

### create spatial dataframe of observed and expected----
vcmap=SpatialPolygonsDataFrame(vcmap,db.comuni,match.ID=F)

### - Create adjacency matrix-----
vcnb=poly2nb(vcmap)

### Create weights------
vcWBweights=nb2WB(vcnb)

## define and run model------
d=list(O=vcmap$OSS,E=vcmap$ATT, N=86)
inits1=list(u=rep(0,86),v=rep(0,86),precu=1,precv=1, alpha=0.5)
#inits2=list(u=rep(1,86),v=rep(1,86),precu=0.1,precv=0.1)

model=function() {
  for(i in 1:N){
    O[i]~dpois(mu[i])
    mu[i]<-theta[i]*E[i]
    log(theta[i])<- alpha+u[i]+v[i]
    u[i]~dnorm(0,precu)
    SMR[i]<-O[i]/E[i]
    prob[i]<-step(theta[i]-1)
  }
  v[1:N]~car.normal(adj[],weights[],num[],precv)
  precu~dgamma(0.001,0.001)
  precv~dgamma(0.001,0.001)
  alpha ~ dflat()
  sigmau<-1/precu
  sigmav<-1/precv
}
mfile=paste(getwd(), "/model.txt", sep = "", collapse = "")
tdir=paste(getwd(), "/Endoutput", sep = "", collapse = "")
dir.create(tdir)
write.model(model,mfile)
res= bugs(data = c(d, vcWBweights),
          inits = list(inits1),
          parameters.to.save = c("u", "v", "theta","prob", "sigmau", "sigmav"),
          model.file = mfile, debug=T,
          n.thin = 3,
          n.chains = 1,
          n.iter = 60000, 
          n.burnin = 20000) 

### add results to map object------
vcmap$prob=res$mean$prob
vcmap$theta=res$mean$theta
vcmap$u=res$mean$u
vcmap$v=res$mean$v
logfile=paste(getwd(), "/Endoutput/log.txt", sep = "", collapse = "")
reslog= bugs.log(file = logfile)

### plot the map-----
fbreaks=c(min(vcmap$theta)-0.1, quantile(vcmap$theta, probs=c(0.2,0.4,0.6,0.8)), max(vcmap$theta)+0.1)
thetaint=classIntervals(vcmap$theta, n=5, style="fixed",
                        fixedBreaks=fbreaks) #fbreaks: parametro modificabile
cols=brewer.pal(5, 'Blues') # parametro modificabile
cktheta=list(labels=as.character(formatC(thetaint$brks, digits=3)),
             at=(0:5)/5, height=.5)
print(spplot(vcmap, "theta", col.regions= cols,
             at=thetaint$brks, axes=TRUE, colorkey=cktheta,
             main=title)) ## Title: parametro modificabile

