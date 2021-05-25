install.package(tseries)
library(forecast)
library(stats)
library(tseries)
library(lmtest)
library(nlme)

#auto.arima()
sncf<-read.table("sncf",header=TRUE,row.names=1)

#
X<-ts(sncf,start = c(1955, 1),freq=12)
plot(X,ylab="Nombre de voyageurs")
abline(v=1961,col="red")

n=length(X)

#Application moyenne mobile
y=NULL
y=c(1/2,rep(1,11),1/2)
T=filter(X,y)
plot(T)
abline(v=72,col="red")


#On s'intÃ©resse a la fenetre 95-215, car on souhaite prÃ©dire 216.
X<-ts(sncf[73:(n-12),],start = c(1961, 1),freq=12)

plot(X)
lines(T,col="red")
layout(1)
#############   ESTIMATION AVEC DECOMPOSE

Decomp=decompose(X)
Res=na.omit(decompose(X)$random)
plot(Res)
acf(Res,50)
pacf(Res,50)
###### Retour d'un (2,0,3)(0,0,1)(12)
auto.arima(Res)

###### Choix aléatoire bon AIC

arima(Res,order=c(0,0,3),seasonal=c(0,1,1))
Residus=residuals(arima(Res,order=c(0,0,3),seasonal=c(0,1,1)))

acf(Residus)
pacf(Residus)


#### Bruit blanc
p.val=NULL
for(i in 1:24) p.val[i]=Box.test(Residus, lag = i, type = c("Box-Pierce"))$p.value
barplot(p.val)
abline(h=0.05,col="red",lty=3)


#### Prédiction SARIMA

PredDec=predict(arima(Res,order=c(0,0,3),seasonal=c(0,1,1)), n.ahead = 12)
plot(na.omit(Decomp$trend),type="l")

#### Approximation linéaire de la tendance
Reg1<-lm(Decomp$trend~time(Decomp$trend))
summary(Reg1)

#### Modèle non réel, puisque les résidus sont auto corrélés
plot(Reg1$residuals)
qqnorm(Reg1$residuals)

acf(Reg1$residuals,100)
pacf(Reg1$residuals,100)
kpss.test(Reg1$residuals)

auto.arima(Reg1$residuals)
#ARIMA(0,1,2)

ResTrend=residuals(arima(Reg1$residuals,order=c(0,1,2)))
acf(ResTrend)

###################################################
#				Résidus régression plus partie stationnaire
###################################################
Res=na.omit(decompose(X)$random)
Reg1<-lm(na.omit(Decomp$trend)~time(na.omit(Decomp$trend)))
Res1=residuals(Reg1)
#
Y=X-coef(Reg1)[1]-coef(Reg1)[2]*time(X)
#
Y=Y-Decomp$seasonal
#
Z=Res+Decomp$random
auto.arima(Z)
#######
plot(diff(Z,1))
adf.test(Z)

#######
arima(Z,order=c(0,0,2))

######


plot(diff(Y,1))
acf(diff(Y,1),50)
pacf(diff(Y,1),50)
adf.test(diff(Y,1))
###
Test<-arima(diff(Y,1),order=c(1,0,0),seasonal=c(1,0,0))
#
plot(residuals(Test))

x=rnorm(1000)
adf.test(x[2:1000]+x[1:999])
###
T=diff(Y,12)
par(mfrow=c(1,2))
plot(diff(Y,1))
plot(diff(T,1))

adf.test(diff(T,1))
kpss.test(diff(T,1))
  ###)
par(mfrow=c(1,2))
plot(Res)
plot(Res1,type="l")
#Résidus du modèle
layout(1)
Z=Res1+Res
mean(Y)



par(mfrow=c(1,2))
plot(Y,ylab="z(t)")
plot(diff(Y),ylab="z(t)-z(t-1)")

#Tests stationnarité

kpss.test(Y)
adf.test(Y)

###
gls(na.omit(Decomp$trend)~time(na.omit(Decomp$trend)),corARMA(c(1.2172,0.3372),0,2))
###

par(mfrow=c(1,2))
plot(diff(Z,1))
acf(diff(Z,1),50)
pacf(diff(Z,1),50)


kpss.test(diff(Y,1))
adf.test(diff(Y,1))
#ACF PACF
par(mfrow=c(1,2))
acf(Y,50)
pacf(Y,50)

auto.arima(Y)

mean(diff(Y))
var(diff(Y))

#Test de différentes valeurs de SARIMA
plot(diff(Z,1))
adf.test(diff(Z,1))
arima(Z,order=c(0,0,2),seasonal=list(order=c(0,0,1)))

#Récupération résidus
ResBB=residuals(arima(Z,order=c(0,0,2),seasonal=list(order=c(0,0,1),period=12),include.mean=FALSE))
plot(ResBB)
acf(ResBB,100)
pacf(ResBB,100)


#### Bruitblanc
layout(1)
p.val=NULL
for(i in 1:24) p.val[i]=Box.test(ResBB, lag = i, type = c("Box-Pierce"))$p.value
barplot(p.val)
abline(h=0.05,col="red",lty=3)


#### Bruitblanc
p.val=NULL
for(i in 1:24) p.val[i]=Box.test(ResBB, lag = i, type = c("Ljung"))$p.value
barplot(p.val)
abline(h=0.05,col="red",lty=3)


####PREDICTION
ArModRes=arima(Z,order=c(1,1,3),include.mean=FALSE)
Prediction=forecast(ArModRes,18,level=95)
summary(Prediction)

Saison=Decomp$seasonal[7:18]
SansTendance = rep(Saison,len=18) + Prediction$mean
plot(SansTendance)
Tend=coef(Reg1)[1]+coef(Reg1)[2]*(1971+(7:24)/12)

Prev=SansTendance+Tend


plot(Prev)
#################### IC 95% ##################
Splus = rep(Saison,len=18) + Prediction$upper
Smoins=rep(Saison,len=18) + Prediction$lower

ICup=Splus + Tend
IClow=Smoins + Tend

#################### GRAPHES #################*
X11()
plot(ts(sncf[145:(n-11),],start = c(1967, 1),freq=12),xlim=c(1970,1973),ylim=c(1500,4500),col="blue",ylab="Nombre d'usagers")
points(1971+(6:23)/12,Prev,col="red",pch=4)
lines(1971+(12:23)/12,ts(ICup[-(1:6),1],freq=12),col="red",lty=2)
lines(1971+(12:23)/12,ts(IClow[-(1:6),1],freq=12),col="red",lty=2)
lines(1971+(12:23)/12,ts(sncf[(n-11):n,],freq=12),col="blue",lty=3)
abline(v=1972)
title("Prévisions à partir d'un ARIMA(1,1,3)")
help(legend)

##########################ETUDE DE X_t - X_(t-12)
X<-ts(sncf[73:(n-12),],start = c(1961, 1),freq=12)
auto.arima(X)
### DIFF



Res=residuals(arima(X,order=c(1,1,1), seasonal=c(0,1,1)))
plot(Res)
acf(Res)
pacf(Res)
layout(1)
kpss.test(Res)
adf.test(Res)

#####

arima(X,order=c(1,1,1),seasonal=c(0,1,1))

#####  
layout(1)
plot(decompose(X)$seasonal)

HW=HoltWinters(X)
HW_forecast=forecast(HW,h=12,level=95)
X11()
plot(ts(sncf[169:(n-11),],start = c(1969, 1),freq=12),xlim=c(1970,1973),ylim=c(1500,4500),col="blue",ylab="Nombre d'usagers")
points(1971+(12:23)/12,HW_forecast$mean,col="red",pch=4)
lines(1971+(12:23)/12,ts(HW_forecast$upper,freq=12),col="red",lty=2)
lines(1971+(12:23)/12,ts(HW_forecast$lower,freq=12),col="red",lty=2)
lines(1971+(12:23)/12,ts(sncf[(n-11):n,],freq=12),col="blue",lty=3)
abline(v=1972)
title("Estimation par lissage de HoltWinters")

#r=12
#cste=1.8
##### Prévision non paramétrique
valid=validation(as.ts(X),hprev=12,where.r=11:13,where.cte=seq(1,2,0.1),methprev=3)


PrevNP=np(X,kprev=12,r=12,cte=1.8,methprev=3)

NP=PrevNP$prevision
ICplus=PrevNP$CS
ICmoins=PrevNP$CI


X11()
plot(ts(sncf[145:(n-11),],start = c(1967, 1),freq=12),xlim=c(1970,1973),ylim=c(1500,4500),col="blue",ylab="Nombre d'usagers")
points(1971+(12:23)/12,NP,col="red",pch=4)
lines(1971+(12:23)/12,ts(ICplus,freq=12),col="red",lty=2)
lines(1971+(12:23)/12,ts(ICmoins,freq=12),col="red",lty=2)
lines(1971+(12:23)/12,ts(sncf[(n-11):n,],freq=12),col="blue",lty=3)
abline(v=1972)
title("Estimation non paramétrique")


auto.arima(X)
S=arima(X,order=c(1,1,1),seasonal=c(0,1,1))
S
plot(residuals(S))
mean(residuals(S))

PredSarima=forecast(S,h=12,level=95)
PredSarima$mean
#RMSEP ARIMA
sqrt(sum((Prev[-(1:6)]-ts(sncf[(n-11):n,]))^2)/12)
#RMSEP HW
sqrt(sum((HW_forecast$mean-sncf[(n-11):n,])^2)/12)
#RMSEP SARIMA 
sqrt(sum((PredSarima$mean-sncf[(n-11):n,])^2)/12)

sqrt(sum((NP-sncf[(n-11):n,])^2/12))
