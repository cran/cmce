 #    gpfuns.r: Supporting functions for EOF Bayesian Calibration model
 #    Copyright (C) 2016  Matthew T. Pratola <mpratola@stat.osu.edu>

 #    This program is free software: you can redistribute it and/or modify
 #    it under the terms of the GNU Affero General Public License as published
 #    by the Free Software Foundation, either version 3 of the License, or
 #    (at your option) any later version.

 #    This program is distributed in the hope that it will be useful,
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #    GNU Affero General Public License for more details.

 #    You should have received a copy of the GNU Affero General Public License
 #    along with this program.  If not, see <http://www.gnu.org/licenses/>.





######################################################################
# Supporting Functions

# rhomat:
# Calculate the correlation matrix for the power exponential/Gaussian model.
#
# The correlation parameters are rho_1,...,rho_p for a p-dim space, each in [0,1] and
# we will have p distance (geoD) matrices.  We construct these in a list of lists, see
# for example the following:
# l1=list(m1=design.distmat.dim1)
# l2=list(m2=design.distmat.dim2)
# l=list(l1=l1,l2=l2)
#
rhomat<-function(geoD,rho,alpha=2)
{
	rho=as.vector(rho)
	if(!is.vector(rho)) stop("non-vector rho!")
	if(any(rho<0)) stop("rho<0!")
	if(any(rho>1)) stop("rho>1!")
	if(any(alpha<1) || any(alpha>2)) stop("alpha out of bounds!")
	if(!is.list(geoD)) stop("wrong format for distance matrices list!")
	if(length(geoD)!=length(rho)) stop("rho vector doesn't match distance list")

	R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=nrow(geoD$l1$m1))
	for(i in 1:length(rho))
		R=R*rho[i]^(abs(geoD[[i]][[1]])^alpha)

	return(list(R=R))
}

# mroot:
# A wrapper for calculating the square root of a matrix using the Cholesky decomposition
mroot<-function(R,eps=1e-10)
{
	nr=dim(R)[1]
	const=1
	if(sum(diag(R))!=nr) #not a correlation matrix
		const=max(abs(R))
	R=R/const
	cR=chol(R+diag(nr)*eps)
	cR=sqrt(const)*cR
	return(t(cR))  # chol gives upper triangular but the definition is lower triangular, hence the transpose
}

# minv:
# A wrapper for calculating the inverse of a matrix using the Choleksy decomposition
minv<-function(R,eps=1e-10)
{
	nr=dim(R)[1]
	const=1
	if(sum(diag(R))!=nr) #not a correlation matrix
		const=max(abs(R))
	R=R/const
	Ri=chol2inv(chol(R+diag(nr)*eps))
	Ri=1/const*Ri
	return(Ri)
}

# mtrace:
# A wrapper for performing the trace of a square matrix
mtrace<-function(R)
{
	return(sum(diag(R)))
}



######################################################################
# Log Likelihoods for GP models

# gp.lp:
# Calculate log likelihood.  Defaults of lambdaf=Inf corresponds to model with no nugget.
gp.lp<-function(data,l.v,mu,lambda,rhos,lambdaf=Inf,eps=1e-10,from="")
{
	n=length(data)
	if(dim(l.v[[1]][[1]])[1]!=n) stop(paste(from,"gp.lp: Dimension of vector ``data'' doesn't match dimension of matrix ``R''\n",sep=""))

	R=rhomat(l.v,rhos)$R+diag(n)*(lambda/lambdaf)
	cR=mroot(R,eps)
	Ri=minv(R,eps)
	Ei=lambda*Ri
	logdetE=-n*log(lambda)+sum(diag(cR))
	lp=-0.5*logdetE-0.5*t(data-mu)%*%Ei%*%(data-mu)

	return(lp)
}


# gp.realize.R:
# Draw a realization from a GP given the correlation matrix R with mean mu, process precision lambda.
gp.realize.R<-function(R,mu,lambda,eps=1e-10,from="")
{
	n=dim(R)[1]
	L=1/sqrt(lambda)*mroot(R,eps)
	u=rnorm(n)
	realization=L%*%u+mu

	return(realization)
}

# gp.realize:
# Draw a realization from a GP at location distances specified by l.v, with mean mu, process precision
# lambda, correlation rhos and error precision lambdaf.
gp.realize<-function(l.v,mu,lambda,rhos,lambdaf=Inf,eps=1e-10,from="")
{
	n=dim(l.v[[1]][[1]])[1]
	if(length(mu)==1) mu=rep(mu,n)
	if(length(mu)!=n) stop(paste(from,"gp.realize: Dimension of mean vector doesn't match dimension of matrix ``R''\n",sep=""))

	R=rhomat(l.v,rhos)$R+diag(n)*(lambda/lambdaf)
	realization=gp.realize.R(R,mu,lambda,eps,from="gp.realize->")

	return(realization)
}

# gp.crealize:
# Draw a conditional realization from a GP at location distances specified by the first m data location
# distances specified by l.v with mean mu, process precision lambda, correlation rhos and error 
# precision lambdaf given data.
# l.v is (m+n)x(m+n), data is nx1, mu is (m+n)x1.
gp.crealize<-function(l.v,data,mu,lambda,rhos,lambdaf=Inf,eps=1e-10,from="")
{
	n=length(data)
	mpn=dim(l.v[[1]][[1]])[1]
	m=mpn-n
	mp1=m+1
	if(length(mu)==1) mu=rep(mu,mpn)
	if(length(mu)!=mpn) stop(paste(from,"gp.crealize: Dimension of mean vector doesn't match dimension of matrix ``R''\n",sep=""))

	R=rhomat(l.v,rhos)$R+diag(mpn)*(lambda/lambdaf)
	mu.1=mu[1:m]
	mu.2=mu[mp1:mpn]
	R11=R[1:m,1:m]
	R22=R[mp1:mpn,mp1:mpn]
	R12=R[1:m,mp1:mpn]
	R22i=minv(R22,eps)
	mu1.2=mu.1+R12%*%R22i%*%(data-mu.2)
	realization=mu1.2+rnorm(m,sd=1/sqrt(lambdaf))

	return(realization)
}



# gp.c2realize2:
# Like gp.crealize except we draw the last m components rather than the first m components.
# Draw a conditional realization from a GP at location distances specified by the last m data location
# distances specified by l.v with mean mu, process precision lambda, correlation rhos and error 
# precision lambdaf given data.
# l.v is (n+m)x(n+m), data is nx1, mu is (n+m)x1.
gp.crealize2<-function(l.v,data,mu,lambda,rhos,lambdaf=Inf,eps=1e-10,from="")
{
	n=length(data)
	npm=dim(l.v[[1]][[1]])[1]
	m=npm-n
	np1=n+1
	if(length(mu)==1) mu=rep(mu,npm)
	if(length(mu)!=npm) stop(paste(from,"gp.crealize: Dimension of mean vector doesn't match dimension of matrix ``R''\n",sep=""))

	R=rhomat(l.v,rhos)$R+diag(npm)*(lambda/lambdaf)
	mu.1=mu[1:n]
	mu.2=mu[np1:npm]
	R11=R[1:n,1:n]
	R22=R[np1:npm,np1:npm]
	R21=R[np1:npm,1:n]
	R11i=minv(R11,eps)
	mu2.1=mu.2+R21%*%R11i%*%(data-mu.1)
	realization=mu2.1+rnorm(m,sd=1/sqrt(lambdaf))

	return(realization)
}




# gp.cisrealize:
# Draw a conditional realization from a GP at location distances specified by the in-sample points
# specified by l.v with mean mu, process precision lambda, correlation rhos and error precision
# lambdaf given data.  In otherwords, predict the underlying process given noisy data.
# l.v is nxn, data is nx1, mu is nx1.
gp.cisrealize<-function(l.v,data,mu,lambda,rhos,lambdaf,eps=1e-10,from="")
{
	n=dim(l.v[[1]][[1]])[1]
	if(length(mu)==1) mu=rep(mu,n)
	if(length(mu)!=n) stop(paste(from,"gp.crealize: Dimension of mean vector doesn't match dimension of matrix ``R''\n",sep=""))
	if(length(data)!=n) stop(paste(from,"gp.crealize: Dimension of data vector doesn't match dimension of matrix ``R''\n",sep=""))

	R=rhomat(l.v,rhos)$R+diag(n)*(lambda/lambdaf)
	Rp=rhomat(l.v,rhos)$R
	R1.2=Rp-Rp%*%minv(R,eps)%*%Rp
	Ri=minv(R,eps)
	mu1.2=mu+Rp%*%Ri%*%(data-mu)
	realization=gp.realize.R(R1.2,mu1.2,lambda,eps,from="gp.cisrealize->")

	return(realization)
}



# makedistlist
# Make list of distance matrices
makedistlist<-function(design)
{
	design=as.matrix(design)
	p=ncol(design)
	l.d=vector("list",p)
	for(i in 1:p) {
		l.d[[i]]=list(outer(design[,i],design[,i],"-"))
		names(l.d[[i]])=paste("m",i,sep="")
		names(l.d)[[i]]=paste("l",i,sep="")
	}
		
	return(l.d)
}

# getranges
# Return ranges of design inputs
getranges<-function(design)
{
	r=t(apply(design,2,range))
	r[,1]=floor(r[,1])
	r[,2]=ceiling(r[,2])

	return(r)
}

# scaledesign
# Rescale the design to the [0,1] hypercube.
scaledesign<-function(design,r)
{
	p=ncol(design)
	for(i in 1:p)
		design[,i]=(design[,i]-r[i,1])/(r[i,2]-r[i,1])

	return(design)
}

# unscale
unscalemat<-function(mat,r)
{
	p=ncol(mat)
	for(i in 1:p)
		mat[,i]=(mat[,i]*(r[i,2]-r[i,1]))+r[i,1]

	return(mat)
}
