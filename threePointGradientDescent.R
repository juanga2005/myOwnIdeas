#Script to apply the three point gradient descent
library(MASS)
set.seed(3)
#Defining the values of gamma
gamma1=0.9;gamma2=0.8;gamma3=0.95

covMatrix=matrix(c(1,gamma1,gamma2,gamma1,1,gamma3,gamma2,gamma3,1),ncol=3)

samples=mvrnorm(1,rep(0,3),covMatrix)


#Gradient descent
#We want to minimize \|(A(gamma)-I)*samples\|

x=samples[1];y=samples[2];z=samples[3]
xx=matrix(samples,ncol=1)



#Building h transpose
ht=matrix(c(y,z,0,x,0,z,0,x,y),ncol=3,byrow=T)

#Initializing the algorithm
#g1=g2=g3=0.5
#G0=c(g1,g2,g3)
set.seed(Sys.time())
G0=runif(3,-2,2)
print('The initial guess is')
print(G0)


A=function(X){
	x=X[1];y=X[2];z=X[3]
	return(matrix(c(0,x,y,x,0,z,y,z,0),ncol=3,byrow=T))
}

norm=function(x){
	return(sqrt(sum(x^2)))
}

niter=0
maxiter=1000
h=0.1
tol=1e-5

while(T){
	
	gradf=t((A(G0)-diag(3))%*%xx)%*%t(ht)

	Gnew=G0-h*gradf

	if(sqrt(sum(gradf^2))<tol||niter==maxiter){
		break
	}
	
	niter=niter+1
	G0=Gnew
}


####Solving for the two linear systems
##View T:R^3\rightarrow R^6

a12=G0[1];a13=G0[2];a23=G0[3]

M=matrix(c(1,0,-a13,0,1,-a12,1,-a23,0,0,-a12,1,-a23,1,0,-a13,0,1),ncol=3,byrow=T)
a=matrix(c(a12,a13,a12,a23,a13,a23),ncol=1)
b=t(M)%*%a
gammaEstimate=solve((t(M)%*%M),b)

##View T:R^3\rightarrow R^2





















