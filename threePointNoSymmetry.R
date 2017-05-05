#Script to apply the three point gradient descent
library(MASS)
set.seed(6)
#Defining the values of gamma
gamma1=0.9;gamma2=0.9;gamma3=0.9

covMatrix=matrix(c(1,gamma1,gamma2,gamma1,1,gamma3,gamma2,gamma3,1),ncol=3)

samples=mvrnorm(1,rep(0,3),covMatrix)
plot(samples)

#Gradient descent
#We want to minimize \|(A(gamma)-I)*samples\|

s1=samples[1];s2=samples[2];s3=samples[3]
xx=matrix(samples,ncol=1)



#Building h'
dh=matrix(c(s2,s3,0,0,0,0,0,0,s1,s3,0,0,0,0,0,0,s1,s2),ncol=6,byrow=T)

#Initializing the algorithm
#g1=g2=g3=0.5
#G0=c(g1,g2,g3)
set.seed(Sys.time())
G0=runif(6,-2,2)
#print('The initial guess is')
#print(G0)


A=function(X){
	x1=X[1];x2=X[2];x3=X[3];x4=X[4];x5=X[5];x6=X[6]
	return(matrix(c(0,x1,x2,x3,0,x4,x5,x6,0),ncol=3,byrow=T))
}

norm=function(x){
	return(sqrt(sum(x^2)))
}

niter=0
maxiter=100000
h=0.1
tol=1e-5

while(T){
	
	gradf=t((A(G0)-diag(3))%*%xx)%*%dh

	Gnew=G0-h*gradf

	if(sqrt(sum(gradf^2))<tol||niter==maxiter){
		break
	}
	
	niter=niter+1
	G0=Gnew
}


####Solving for the two linear systems
##View T:R^3\rightarrow R^6

x1=G0[1];x2=G0[2];x3=G0[3];x4=G0[4];x5=G0[5];x6=G0[6]

M=matrix(c(1,0,-x2,0,1,-x1,1,-x4,0,0,-x3,1,-x6,1,0,-x5,0,1),ncol=3,byrow=T)
a=matrix(c(x1,x2,x3,x4,x5,x6),ncol=1)
b=t(M)%*%a
gammaEstimate=solve((t(M)%*%M),b)

#The estimate for gamma is
print('The estimate for gamma is')
print(gammaEstimate)

##Second view T:R^3\rightarrow R^2


estGamma1=x1-gamma3*x2
estGamma2=+x2-gamma3*x1


















