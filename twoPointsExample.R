#Script to explore the dependence between two points 
#Assume x1,x2\sim N(0,sigma) with sigma=[1 gamma;gamma 1]


library(MASS)

gamma=0.99
covMatrix=matrix(c(1,gamma,gamma,1),ncol=2)


nSamples=100000
samples=mvrnorm(nSamples,mu=c(0,0),Sigma=covMatrix)


#Errors
linearApprox=gamma*samples[,1]
e=samples[,2]-linearApprox

#\|x2-gamma*x1\|^2=\int(x2-x1)^2dP
f=sum(e^2)/nSamples
print(f)

#Estimating the covariance
covEstimate=samples[,2]/samples[,1]
covEstimate=matrix(rep(0,4),ncol=2)
for (j in 1:(nSamples/10)){
	covEstimate2=t(t(samples[j,]))%*%samples[j,]+covEstimate2
}
covEstimate2=10/nSamples*covEstimate2
