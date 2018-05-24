 #    calfuns.r: Supporting functions for EOF Bayesian Calibration model
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
# Support for correlation parameter Gibbs and MH steps

# dlp.rho:
# This function is for calculating the difference in log posterior for the draw of
# the rho's for the l'th component of the V's.
# dlp.rho = lp.rho(new)-lp.rho(old)
dlp.rho<-function(l,vl,lambdav,rhol,rhol.new,pinfo)
{
	l.v=pinfo$l.v
	for(i in 1:pinfo$q)
		l.v[[i]][[1]]=l.v[[i]][[1]][1:pinfo$m,1:pinfo$m]
	return( lp.rho(l,vl,lambdav,rhol.new,l.v,pinfo)-lp.rho(l,vl,lambdav,rhol,l.v,pinfo) )
}

# lp.rho:
# Calculate the part of the log posterior needed for the draw of the rho's
lp.rho<-function(l,vl,lambdav,rhol,l.v,pinfo)
{
	R=rhomat(l.v,rhol)$R
	cR=mroot(R,pinfo$eps)
	Ri=minv(R,pinfo$eps)
	logdetR=2*sum(log(diag(cR)))

	logpost = (length(vl))/2.0*log(lambdav) - 0.5*lambdav*t(vl)%*%Ri%*%(vl) - 0.5*logdetR 
			+ (pinfo$rho$a[l]-1)*log(rhol[l]) + (pinfo$rho$b[l]-1)*log(1-rhol[l])

	return(logpost)
}




######################################################################
# Support for precision parameters (Gibbs steps)

# d.lambdav:
# Draw the marginal precision of the right eigenfunctions (Gibbs step)
d.lambdav<-function(V,theta,rhos,pinfo)
{
	lambdav=rep(0,pinfo$nc)
	l.v=pinfo$l.v
	for(i in 1:pinfo$q)
		l.v[[i]][[1]]=l.v[[i]][[1]][1:pinfo$m,1:pinfo$m]

	for(l in 1:pinfo$nc)
	{
		R=rhomat(l.v,rhos[l,])$R
		Ri=minv(R,pinfo$eps)
		lambdav[l]=rgamma(1,pinfo$lambdav$a[l]+(pinfo$m)/2,
							pinfo$lambdav$b[l]+0.5*t(V[1:pinfo$m,l])%*%Ri%*%V[1:pinfo$m,l])
	}

	return(lambdav)
}




######################################################################
# Drawing the theta's

# dlp.theta:
# This function is for calculating the difference in log posterior for the draw of theta
# theta is q x 1
# theta.new is q x 1
# lambdav is nc x 1
# rhos is nc x q
# lambda is scalar
# u is n x nc
# v is m+1 x nc
# v.new is m+1 x nc
# lambdaf is scalar
# phi is the n x m+1 matrix of simulations
# yf is the n x 1 observation vector
# delta is the n x 1 discrepancy
dlp.theta<-function(theta,theta.new,lambdav,rhos,lambda,u,v,v.new,phi,phi.new,pinfo)
{
	return( lp.theta(theta.new,lambdav,rhos,lambda,u,v.new,phi.new,pinfo) - 
			lp.theta(theta,lambdav,rhos,lambda,u,v,phi,pinfo) )
}
# Calculate the log posterior needed to draw theta
lp.theta<-function(theta,lambdav,rhos,lambda,u,v,phi,pinfo)
{
	lp=0
	# first the V_.,k parts
	for(l in 1:pinfo$nc)
		lp=lp+lp.v(l,lambdav[l],rhos[l,],lambda,u,v,phi,theta,pinfo)

	# and then the theta prior part
	lp=lp+log(theta>pinfo$theta$l && theta<pinfo$theta$u)

	return(lp)
}

# pred.v.cond:
# This function is for drawing a realization for V_0,k|V_1:m,k where k is the k \in 1:nc is the k'th component
# Really only used as part of the MH proposal for theta...
# lambdav is scalar
# rhos is q x 1 corrlation parameters for the k'th component
# lambda is scalar
# U is n x nc
# V is m+1 x nc
# phi is the n x m+1 matrix of simulations
# theta is q x 1
pred.v.cond<-function(k,lambdav,rhos,lambda,U,V,phi,theta,pinfo)
{
	m=pinfo$m
	Vbar=apply(V,1,mean)

	d=phi-t(V[,-k]%*%t(U[,-k]))
	dbar=apply(U[,k]*d,2,mean)
	Sbar=mean(U[,k]^2)
	dtilde=dbar/Sbar

	R=rhomat(set.lv(theta,pinfo$l.v),rhos)$R
	Ri=minv(R,pinfo$eps)
	Ei=lambdav*Ri+pinfo$n*lambda*Sbar*diag(m+1)
	E=minv(Ei,pinfo$eps)
	mu=E%*%(dtilde*pinfo$n*lambda*Sbar+lambdav*Ri%*%Vbar)

 	# conditional mean and variance
 	Emi=minv(E[1:m,1:m],pinfo$eps)
	mu.v=mu[m+1]+E[m+1,1:m]%*%Emi%*%(V[1:m,k]-mu[1:m])
	E.v=E[m+1,m+1]-E[m+1,1:m]%*%Emi%*%E[1:m,m+1]

	# draw realization
	E.v=max(0,E.v) # we might sometimes get underflow, so bound it below by 0.
	v.cond=mu.v+rnorm(1,sd=sqrt(E.v))

	return(v.cond)
}

# pred.vm.cond:
# This function is for drawing a realization for V_1:m,k|V_0,k where k is the k \in 1:nc is the k'th component
# lambdav is scalar
# rhos is q x 1 corrlation parameters for the k'th component
# lambda is scalar
# U is n x nc
# V is m+1 x nc
# phi is the n x m+1 matrix of simulations
# theta is q x 1
pred.vm.cond<-function(k,lambdav,rhos,lambda,U,V,phi,theta,pinfo)
{
	m=pinfo$m
	Vbar=apply(V,1,mean)

	d=phi-t(V[,-k]%*%t(U[,-k]))
	dbar=apply(U[,k]*d,2,mean)
	Sbar=mean(U[,k]^2)
	dtilde=dbar/Sbar

	R=rhomat(set.lv(theta,pinfo$l.v),rhos)$R
	Ri=minv(R,pinfo$eps)
	Ei=lambdav*Ri+pinfo$n*lambda*Sbar*diag(m+1)
	E=minv(Ei,pinfo$eps)
	mu=E%*%(dtilde*pinfo$n*lambda*Sbar+lambdav*Ri%*%Vbar)

 	# conditional mean and variance
 	Emi=matrix(1/E[m+1,m+1],1,1)
	mu.v=mu[1:m]+E[1:m,m+1]%*%Emi%*%(V[m+1,k]-mu[m+1])
	E.v=E[1:m,1:m]-E[1:m,m+1]%*%Emi%*%E[m+1,1:m]

	# draw realization
	#vm.cond=mu.v
	vm.cond=gp.realize.R(E.v,mu.v,1,eps=1e-10,from="pred.vm.cond->")

	return(vm.cond)
}

# lp.v:
# Calculate the log probability for V[,k] given theta and everything else.
# Used in the draw of theta, for instance.
lp.v<-function(k,lambdav,rhos,lambda,U,V,phi,theta,pinfo)
{
	m=pinfo$m
	Vbar=apply(V,1,mean)

	d=phi-t(V[,-k]%*%t(U[,-k]))
	dbar=apply(U[,k]*d,2,mean)
	Sbar=mean(U[,k]^2)
	dtilde=dbar/Sbar

	R=rhomat(set.lv(theta,pinfo$l.v),rhos)$R
	Ri=minv(R,pinfo$eps)
	Ei=lambdav*Ri+pinfo$n*lambda*Sbar*diag(m+1)
	E=minv(Ei,pinfo$eps)
	mu=E%*%(dtilde*pinfo$n*lambda*Sbar+lambdav*Ri%*%Vbar)

	cE=mroot(E,pinfo$eps)
	logdetE=2*sum(log(diag(cE)))
	lp=-0.5*logdetE-0.5*t(V[,k]-mu)%*%Ei%*%(V[,k]-mu)

	return(lp)
}

# difference in log probability of yf after marginalizing out delta
dlp.yfmd<-function(yf,phi0,phi0.new,mu.delta,lambdaf,lambdad,psis,pinfo)
{
	return( lp.yfmd(yf,phi0.new,mu.delta,lambdaf,lambdad,psis,pinfo) - 
			lp.yfmd(yf,phi0,mu.delta,lambdaf,lambdad,psis,pinfo) )	
}

# log probability for yf given rest but after marginalizing out delta with respect to the
# GP prior on delta.
lp.yfmd<-function(yf,phi0,mu.delta,lambdaf,lambdad,psis,pinfo)
{
	n=pinfo$n
	I=diag(n)
	lamf=rep(0,n)
	for(j in 1:pinfo$ns)
		lamf[pinfo$six[[j]]]=lambdaf[j]

	R=I
	for(i in 1:pinfo$ns)
		R[pinfo$six[[i]],pinfo$six[[i]]]=(1/lambdad[i])*R[pinfo$six[[i]],pinfo$six[[i]]]%*%rhomat(pinfo$l.d,psis[i,])$R[,,drop=FALSE]
			+ (1/pinfo$lambdakap[i])*phi0[pinfo$six[[i]]]%*%t(phi0[pinfo$six[[i]]])
#		R[pinfo$six[[i]],pinfo$six[[i]]]=(1/lambdad[i])*R[pinfo$six[[i]],pinfo$six[[i]]]%*%rhomat(pinfo$l.d,psis[i,])$R[pinfo$six[[i]],pinfo$six[[i]],drop=FALSE]
# bugfix


#	E=(1/lambdaf)*I + R
	E=diag(1/lamf)+R
	Ei=minv(E,pinfo$eps)
#	mu=mu.delta+phi0
	mu=rep(0,n)
	for(i in 1:pinfo$ns)
		mu[pinfo$six[[i]]]=mu.delta[pinfo$six[[i]]]+pinfo$mukap[i]*phi0[pinfo$six[[i]]]
	lp=-0.5*t(yf-mu)%*%Ei%*%(yf-mu)

	return(lp)
}

# difference in log probability of the prior on one calibration parameter, theta.
# k denotes which calibration parameter to make the calculation for.
dlp.th<-function(k,thetanew,theta,pinfo)
{
	return(pinfo$thetaprior[[k]](thetanew,pinfo$ranges) -
		   pinfo$thetaprior[[k]](theta,pinfo$ranges))
}


######################################################################
# Drawing the U functions and precisions.

# pred.u:
# This function is drawing a realization of U_1:n,k|rest where k is the k \in 1:nc is the k'th component
# lambdau is scalar
# lambda is scalar
# u is n x nc
# v is m+1 x nc
# phi is the n x m+1 matrix of simulations
pred.u<-function(k,lambdau,lambda,U,V,phi,pinfo)
{
	m=pinfo$m
	n=pinfo$n
	Ubar=apply(U,1,mean)
	d=phi-U[,-k]%*%t(V[,-k])
	dbar=apply(d*(t(V[,k])%x%rep(1,n)),1,mean)
	Sbar=mean(V[,k]^2)

	I=diag(1,n)
	E=1/(lambda*(m+1)*Sbar+lambdau)*I
	L=1/sqrt(lambda*(m+1)*Sbar+lambdau)*I

	# draw realization
	mu=E%*%(lambda*(m+1)*dbar+lambdau*Ubar)
	z=rnorm(n)
	u.pred=mu+L%*%z

	return(u.pred)
}




######################################################################
# Drawing the U functions and precisions (Gibbs step).

# d.lambda:
# Draw the lambda parameter.  This is really supposed to be numerical error
# so lambda should be very large (small error variance).
# Here we again only take the simulator outputs to estimate this to ensure stability.
d.lambda<-function(phi,U,V,pinfo)
{
	n=pinfo$n
	m=pinfo$m
	lambda=rgamma(1,pinfo$lambda$a+n*m/2,pinfo$lambda$b
			+0.5*sum((phi[,1:m]-U%*%t(V[1:m,]))^2) )

	return(lambda)
}




######################################################################

# pred.phi0:
# This function is drawing a realization of Phi0 =equiv= Phi(theta0)
# lambda is scalar
# lambdaf is scalar
# v is n x nc
# u is m+1 x nc
# yf is n x 1 field observations
# delta is n x 1 discrepancy vector
pred.phi0<-function(lambda,lambdaf,U,V,yf,delta,pinfo)
{
	n=pinfo$n
	m=pinfo$m
	df=yf-delta
	I=diag(n)
	lamf=rep(0,n)
	for(j in 1:pinfo$ns)
		lamf[pinfo$six[[j]]]=lambdaf[j]
	If=diag(lamf)

#	E=1/(lambda+lambdaf)*I
	E=diag(1/(lambda+lamf))

	# draw realization
#	mu=E%*%(lambda*t(V[m+1,]%*%t(U)) + lambdaf*df)
	mu=E%*%(lambda*t(V[m+1,]%*%t(U)) + If%*%df)
	L=sqrt(E)
	z=rnorm(n)
	phi0.pred=mu+L%*%z

	return(phi0.pred)
}

dlp.phi0<-function(phi0.new,phi0,lambda,lambdaf,U,V.new,V,yf,delta,pinfo)
{
	return( lp.phi0(phi0.new,lambda,lambdaf,U,V.new,yf,delta,pinfo)
			-lp.phi0(phi0,lambda,lambdaf,U,V,yf,delta,pinfo) )
}

lp.phi0<-function(phi0,lambda,lambdaf,U,V,yf,delta,pinfo)
{
	n=pinfo$n
	m=pinfo$m
	df=yf-delta
	I=diag(n)
	E=1/(lambda+lambdaf)*I

	mu=E%*%(lambda*t(V[m+1,]%*%t(U)) + lambdaf*df)
	lp=-0.5*(lambda+lambdaf)*t(phi0-mu)%*%(phi0-mu)

	return(lp)
}

######################################################################
# Handling the discrepancy

# pred.delta:
# Draw realization of discrepancy
pred.delta<-function(phi0,yf,lambdad,psis,lambdaf,kappa,pinfo)
{
	n=pinfo$n
	mu.delta=pinfo$inidelta
#	d=yf-phi0
	d=rep(0,n)
	for(i in 1:pinfo$ns)
		d[pinfo$six[[i]]]=yf[pinfo$six[[i]]]-kappa[i]*phi0[pinfo$six[[i]]]

	I=diag(n)
	lamf=rep(0,n)
	for(j in 1:pinfo$ns)
		lamf[pinfo$six[[j]]]=lambdaf[j]
	If=diag(lamf)

#	R=rhomat(pinfo$l.d,psis)$R
	R=I
	for(i in 1:pinfo$ns)
		R[pinfo$six[[i]],pinfo$six[[i]]]=R[pinfo$six[[i]],pinfo$six[[i]]]%*%rhomat(pinfo$l.d,psis[i,],alpha=pinfo$delta.corrmodel)$R[,,drop=FALSE]
#		R[pinfo$six[[i]],pinfo$six[[i]]]=R[pinfo$six[[i]],pinfo$six[[i]]]%*%rhomat(pinfo$l.d,psis[i,],alpha=pinfo$delta.corrmodel)$R[pinfo$six[[i]],pinfo$six[[i]],drop=FALSE]
# bugfix

	Ri=minv(R,pinfo$eps)
	lam=NULL
	for(i in 1:pinfo$ns)
		lam=c(lam,rep(lambdad[i],length(pinfo$six[[i]])))
	lam=diag(lam)
#	Ei=lambdaf*I+lam%*%Ri
	Ei=If+lam%*%Ri
	E=minv(Ei,pinfo$eps)

	# draw realization
	L=mroot(E,pinfo$eps)
#	mu=E%*%(lambdaf*I%*%d+lam%*%Ri%*%mu.delta)
	mu=E%*%(If%*%d+lam%*%Ri%*%mu.delta)
	z=rnorm(n)
	delta.pred=mu+L%*%z

	return(delta.pred)
}

dlp.delta<-function(delta.new,delta,phi0.new,phi0,yf,lambdad,psis,lambdaf,pinfo)
{
	return( lp.delta(delta.new,phi0.new,yf,lambdad,psis,lambdaf,pinfo)
			-lp.delta(delta,phi0,yf,lambdad,psis,lambdaf,pinfo) )
}

lp.delta<-function(delta,phi0,yf,lambdad,psis,lambdaf,pinfo)
{
	n=pinfo$n
	mu.delta=pinfo$inidelta
	d=yf-phi0
	I=diag(n)

# bugfix
#	R=rhomat(pinfo$l.d,psis,alpha=pinfo$delta.corrmodel)$R
	R=I
	for(i in 1:pinfo$ns)
		R[pinfo$six[[i]],pinfo$six[[i]]]=R[pinfo$six[[i]],pinfo$six[[i]]]%*%rhomat(pinfo$l.d,psis[i,],alpha=pinfo$delta.corrmodel)$R[,,drop=FALSE]

	Ri=minv(R,pinfo$eps)
	Ei=lambdaf*I+lambdad*Ri
	E=minv(Ei,pinfo$eps)
	mu=E%*%(lambdaf*I%*%d+lambdad*Ri%*%mu.delta)
	lp=-0.5*t(delta-mu)%*%Ei%*%(delta-mu)

	return(lp)
}

# d.lambdad:
# Draw realization of the discrepancies precision
d.lambdad<-function(delta,psis,pinfo)
{
	lambdad=rep(0,pinfo$ns)
	for(i in 1:pinfo$ns) {
		n=length(pinfo$six[[i]])
		mu.delta=pinfo$inidelta[pinfo$six[[i]]]
#		R=rhomat(pinfo$l.d,psis[i,],alpha=pinfo$delta.corrmodel)$R[pinfo$six[[i]],pinfo$six[[i]],drop=FALSE]
# bugfix
		R=rhomat(pinfo$l.d,psis[i,],alpha=pinfo$delta.corrmodel)$R[,,drop=FALSE]
		Ri=minv(R,pinfo$eps)
		lambdad[i]=rgamma(1,pinfo$lambdad$a[i]+n/2,
			pinfo$lambdad$b[i]+0.5*t(delta[pinfo$six[[i]]]-mu.delta)%*%Ri%*%(delta[pinfo$six[[i]]]-mu.delta))
	}

	return(lambdad)
}

# dlp.psi:
# This function is for calculating the difference in log posterior for the draw of
# the psi's.
# dlp.psi = lp.psi(new)-lp.psi(old)
dlp.psi<-function(s,p,delta,lambdad,psi,psi.new,pinfo)
{
	return( lp.psi(s,p,delta,lambdad,psi.new,pinfo)-lp.psi(s,p,delta,lambdad,psi,pinfo) )
}

# lp.psi:
# Calculate the part of the log posterior needed for the draw of the psi's
lp.psi<-function(s,p,delta,lambdad,psi,pinfo,eps=pinfo$eps)
{
	n=length(pinfo$six[[s]])
	data=delta[pinfo$six[[s]]]
	mu=pinfo$inidelta[pinfo$six[[s]]]
#	R=rhomat(pinfo$l.d,psi,alpha=pinfo$delta.corrmodel)$R[pinfo$six[[s]],pinfo$six[[s]],drop=FALSE]
# bugfix
	R=rhomat(pinfo$l.d,psi,alpha=pinfo$delta.corrmodel)$R[,,drop=FALSE]
	cR=mroot(R,eps)
	Ri=minv(R,eps)
	Ei=lambdad[s]*Ri
	logdetE=-n*log(lambdad[s])+sum(diag(cR))
	logpost = (pinfo$psis$a[p]-1)*log(psi[p]) + (pinfo$psis$b[p]-1)*log(1-psi[p])
			  -0.5*logdetE -0.5*t(data-mu)%*%Ei%*%(data-mu)

#	logpost = (pinfo$psis$a[p]-1)*log(psi[p]) + (pinfo$psis$b[p]-1)*log(1-psi[p])
#			  + gp.lp(delta,pinfo$l.d,pinfo$inidelta,lambdad,psi,eps=pinfo$eps,from="lp.psi->")

	return(logpost)
}


######################################################################
# Draw the multiplicative discrepancy.  
d.kappa<-function(l,phi0,yf,delta,lambdaf,pinfo)
{
	d=yf[pinfo$six[[l]]]-delta[pinfo$six[[l]]]
	phi0=phi0[pinfo$six[[l]]]

	E=1/(lambdaf[l]*sum(phi0^2)+pinfo$lambdakap[l])
	mu=E*lambdaf[l]*sum(d*phi0)+pinfo$mukap[l]/(lambdaf[l]/pinfo$lambdakap[l]*sum(phi0^2)+1)

	return(rnorm(1,mean=mu,sd=sqrt(E)))
}


######################################################################
# Finally the precision for the field measurement error.

# d.lambdaf:
# Draw Gibbs sample of the field precision.
d.lambdaf<-function(phi0,yf,delta,pinfo)
{
	n=pinfo$n

	lambdaf=rgamma(1,pinfo$lambdaf$a+n/2,pinfo$lambdaf$b
			+ 0.5*t(yf-phi0-delta)%*%(yf-phi0-delta) )

	return(lambdaf)
}

# d.lambdaf:
# Draw Gibbs sample of the field precision.
# This version is for when there are multiple states.
d.lambdaf<-function(j,phi0,yf,delta,kappa,pinfo)
{
	isj=pinfo$six[[j]] #indexing for state j
	n=length(isj)

	lambdaf=rgamma(1,pinfo$lambdaf$a[j]+n/2,pinfo$lambdaf$b[j]
			+ 0.5*t(yf[isj]-kappa[j]*phi0[isj]-delta[isj])%*%(yf[isj]-kappa[j]*phi0[isj]-delta[isj]) )

	return(lambdaf)
}



######################################################################
# Update distance structures given a value for the calibration parameter theta.

# set.lv:
# Set the distance matrices according to value of theta.
set.lv<-function(theta,l.v)
{
	m=nrow(l.v[[1]][[1]])-1
	for(j in 1:length(theta))
	{
		l.v[[j]][[1]][m+1,]=l.v[[j]][[1]][m+1,]+theta[j]
		l.v[[j]][[1]][,m+1]=l.v[[j]][[1]][,m+1]-theta[j]
	}
	return(l.v)
}

# unset.lv:
# Unset the distance matrices according to the previously set value of theta.
unset.lv<-function(theta,l.v)
{
	m=nrow(l.v[[1]][[1]])-1
	for(j in 1:length(theta))
	{
		l.v[[j]][[1]][m+1,]=l.v[[j]][[1]][m+1,]-theta[j]
		l.v[[j]][[1]][,m+1]=l.v[[j]][[1]][,m+1]+theta[j]
	}
	return(l.v)
}


