#Script to experiment with covariances known.



#Script to apply the three point gradient descent
library(MASS)
#set.seed(7)
#Defining the values of gamma
gamma1=0.1;gamma2=0.1;gamma3=-0.9

covMatrix=matrix(c(1,gamma1,gamma2,gamma1,1,gamma3,gamma2,gamma3,1),ncol=3)

samples=mvrnorm(1,rep(0,3),covMatrix)
#plot(samples)

##By doing the analysis we have 3 system of 2 by 2
#Equations to solve for the values of a,b,c
#See file GaussianCovFreeIdeas.xoj page 2


#Creating the systems of equations to solve for a b c

M1=matrix(c(1,gamma3,gamma3,1),ncol=2,byrow=T) #Matrix for a,b
b1=matrix(c(gamma1,gamma2),ncol=1)

M2=matrix(c(1,gamma2,gamma2,1),ncol=2,byrow=T) #Matrix for c,d
b2=matrix(c(gamma1,gamma3),ncol=1)

M3=matrix(c(1,gamma1,gamma1,1),ncol=2,byrow=T) #Matrix for e,f
b3=matrix(c(gamma2,gamma3),ncol=1)

#Finding the values for a,b,c,...
ab=solve(M1,b1)
cd=solve(M2,b2)
ef=solve(M3,b3)

#linear estimates
x1=ab[1,1]*samples[2]+ab[2,1]*samples[3]
x2=cd[1,1]*samples[1]+cd[2,1]*samples[3]
x3=ef[1,1]*samples[1]+ef[2,1]*samples[2]


#################Up to this point we haven't used any information of GP #######################





#Plotting the real samples with the linear estimates
plot(samples,ylim=c(-3,3))
points(c(x1,x2,x3),pch=10,col='blue')
