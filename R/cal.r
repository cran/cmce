#    cal.r: EOF Bayesian Calibration for Stochastic Simulators    
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
# Calibration and Prediction functions

# calibrate:
#' Runs the MCMC algorithm
#'
#' \code{calibrate} runs the MCMC algorithm to calibrate simulator outputs
#' to field observations while estimating unknown settings of calibration
#' parameters and quantifying model discrepancy.  Currently, two forms
#' of discrepancy are supported: additve stationary Gaussian Process discrepancy
#' and scalar multiplicative discrepancy (with a Normal prior).
#'
#' The method can calibrate both stochastic simulator outputs and deterministic
#' simulator outputs, and samples from the posterior distribution of unknowns
#' conditional on the observed field data and simulator outputs.  The method
#' makes use of empirical orthogonal functions to reduce the dimension of
#' the data to make computations feasible.  In addition, there is support for 
#' multi-state simulators, and calibration can be performed when all states
#' or only a subset of states are observed in the field.
#'
#' For more details on the model, see Pratola and Chkrebtii (2016).  Detailed
#' examples demonstrating the method are available at
#' \url{http://www.matthewpratola.com/software}.
#' @param yf A vector of field observations.
#' @param phi A matrix or list of matrices representing simulator outputs
#' @param N The number of MCMC iterations to perform.
#' @param pinfo A list of prior parameter settings.
#' @param mh A list of settings for Metropolis-Hastings proposals.
#' @param last The number of MCMC steps to report (default is 1000).  The first N-last steps are discarded as burn-in.
#' @return A list with "last" samples drawn from the posterior.
# @examples /inst/examples/doc.cal.r
# @examples calibrate(NULL,NULL,1,NULL,NULL)
#' @export
calibrate<-function(yf,phi,N,pinfo,mh,last=1000)
{

	if(N<2000) stop("Minimum number of MCMC iterations required is N>=2000\n")

	if(!is.list(phi)) #coerce it to a list
	{
		phiold=phi
		phi=vector("list",N)
		for(i in 1:N) phi[[i]]=phiold
		rm(phiold)
	}
	else if(length(phi)<N) # we have to resample the simulator outputs
	{
		phiold=phi
		phi=vector("list",N)
		nouts=length(phiold)
		for(i in 1:N)
		{
			phi[[i]]=phiold[[1]]
			samp=sample(1:nouts,m,replace=TRUE)
			for(j in 1:m)
				phi[[i]][,j]=phiold[[samp[j]]][,j]
		}
		rm(phiold)
	}
	if(length(phi)!=N) stop("You don't have N Draws of Phi!\n")

	pinfo$delta.corrmodel=2 #default to gaussian correlation for additive discrepancy
	if(pinfo$delta.corrmodel=="exponential") pinfo$delta.corrmodel=1

	nc=pinfo$nc
	ns=pinfo$ns
	last=last-1
	n=length(yf)
	if(n != length(phi[[1]][,1])) stop("Field grid and model grid don't match!\n")
	m=pinfo$m  # number of simulation model outputs, 1:m.  The m+1'th entry is the unobserved calibrated model.
	q=pinfo$q  # number (dimension) of calibration parameters is theta_1...theta_q
	p=pinfo$p  # dimension of x covariates is x_1...x_p
	l.v=pinfo$l.v

	if(is.null(pinfo$thetaprior)) {
		pinfo$thetaprior=vector("list",q)
		for(i in 1:q) pinfo$thetaprior[[i]]<-function(theta,r){ return(0) }
	}

	## SETUP ######################################################################################################################
	###############################################################################################################################
	# lambda parameters
	draw.lambdaf=matrix(NA,nrow=N,ncol=ns)   # lambdaf parameter for field data
	draw.lambdad=matrix(NA,nrow=N,ncol=ns)   # lambdad parameter for discrepancy
	draw.kappa=matrix(NA,nrow=N,ncol=ns)     # multiplicative discrepancy for the states
	draw.lambdav=matrix(NA,nrow=N,ncol=nc)   # lambdav parameter for each retained right eigenfunction (there are nc of them)
	# rho parameters for each retained right eigenfunction (there are nc of them)
	draw.rhos=vector("list",N)
	for(i in 1:N) draw.rhos[[i]]=matrix(NA,nrow=nc,ncol=q)
	# psi parameters for the discrepancy GP model (there are p of them)
	draw.psis=vector("list",N)
	for(i in 1:N) draw.psis[[i]]=matrix(NA,nrow=ns,ncol=p)
	# theta parameters
	draw.theta=matrix(NA,nrow=N,ncol=q)
	# Posterior draws of discrepancy
	draw.delta=matrix(NA,nrow=N,ncol=n)
	# The U,V draws
	draw.u=vector("list",N)
	for(i in 1:N) draw.u[[i]]=matrix(NA,nrow=n,ncol=nc)
	draw.v=vector("list",N)
	for(i in 1:N) draw.v[[i]]=matrix(NA,nrow=m+1,ncol=nc)


	## INITIALIZE #################################################################################################################
	###############################################################################################################################
	# Initialize U,V
	cat("Calculating Empirical Orthogonal Functions...\n")
	temp=matrix(0,ncol=m,nrow=n*N)
	for(i in 1:N)
		temp[((i-1)*n+1):(i*n),]=phi[[i]][,1:m]
	s=svd(temp)
	for(i in 1:N)
		draw.u[[i]]=s$u[((i-1)*n+1):(i*n),1:nc]%*%diag(sqrt(s$d[1:nc]))
	draw.v[[1]][1:m,]=s$v[,1:nc]%*%diag(sqrt(s$d[1:nc]))
	rm(s)
	
	draw.v[[1]][m+1,]=apply(draw.v[[1]][1:m,],2,mean)

	# Initialize lambda's
	for(i in 1:ns) draw.lambdaf[1,i]=rgamma(1,pinfo$lambdaf$a[i],pinfo$lambdaf$b[i])
	for(i in 1:nc) draw.lambdav[1,i]=rgamma(1,pinfo$lambdav$a[i],pinfo$lambdav$b[i])
	for(i in 1:ns) draw.lambdad[1,i]=rgamma(1,pinfo$lambdad$a[i],pinfo$lambdad$b[i])
	for(i in 1:ns) draw.kappa[1,i]=rnorm(1,pinfo$mukap[i],1/sqrt(pinfo$lambdakap[i]))

	# Initialize GP correlation parameters
	for(i in 1:nc) for(j in 1:q) draw.rhos[[1]][i,j]=0.1
	for(i in 1:ns) for(j in 1:p) draw.psis[[1]][i,j]=0.1

	# Initialize discrepancy
	draw.delta[1,]=pinfo$inidelta  # take the mean from the prior

	# Initialize theta, v(theta)
	draw.theta[1,]=rep(0.5,q)  # we assume the theta's are scaled to 0,1
	for(l in 1:nc)
		draw.v[[1]][m+1,l]=pred.v.cond(l,draw.lambdav[1,l],draw.rhos[[1]][l,],1e15,#draw.lambda[1],
									   draw.u[[1]],draw.v[[1]],phi[[1]],draw.theta[1,],pinfo)

	# Initialize phi(theta), the unknown simulator output
	phi[[1]][,m+1]=pred.phi0(1e15,#draw.lambda[1],
							 draw.lambdaf[1,],draw.u[[1]],draw.v[[1]],yf,
							 draw.delta[1,],pinfo)

	# Stepwidths for correlation parameters and theta
	rp=vector("list",nc)  # steps for the rho parameters for each retained right eigenfunction (there are nc of them)
	for(i in 1:nc) rp[[i]]=rep(mh$rp,q)
	#rr=rep(mh$rr,p)       # steps for the psi parameters for the discrepancy GP (there are p of them)
	rr=vector("list",ns)
	for(i in 1:ns) rr[[i]]=rep(mh$rr,p)
	rth=rep(mh$rth,q)     # steps for the calibration pameters (there are q of them)


	# Accept/reject ratio trackers
	accept.rhos=vector("list",nc)
	for(i in 1:nc) accept.rhos[[i]]=rep(1,q)
	#accept.psis=rep(1,p)
	accept.psis=vector("list",ns)
	for(i in 1:ns) accept.psis[[i]]=rep(1,p)
	accept.theta=rep(1,q)
	lastadapt=0
  
	cat("\n Bayesian Empirical Function Calibration Model")
	cat("\n The last ",last," samples from the posterior will be reported")
	cat("\n Stepwidth for right eigenfunction correlations, rr=",mh$rp)
	cat("\n Stepwidth for discrepancy correlations, rr=",mh$rr)
	cat("\n Stepwidth for calibration parameters, rth=",mh$rth)
	cat("\n ----------------------------------------------------------------\n\n\n")


	## MCMC #######################################################################################################################
	###############################################################################################################################
	# Main MCMC loop
	for(i in 2:N)
	{
		# Draw the correlation parameters for the V's (Metropolis-Hastings steps)
		for(l in 1:nc)
		{
			draw.rhos[[i]][l,]=draw.rhos[[i-1]][l,]
			for(j in 1:q)
			{
				rhos.new=draw.rhos[[i]][l,]
				rhos.new[j]=runif(1,draw.rhos[[i]][l,j]-rp[[l]][j],draw.rhos[[i]][l,j]+rp[[l]][j])
				a=-Inf
				if(rhos.new[j]>0 && rhos.new[j]<1)
					a=dlp.rho(l,draw.v[[i-1]][1:m,l],draw.lambdav[i-1,l],draw.rhos[[i]][l,],rhos.new,pinfo)

				a=min(0,a)
				if(log(runif(1))<a)
				{
					draw.rhos[[i]][l,j]=rhos.new[j]
					accept.rhos[[l]][j]=accept.rhos[[l]][j]+1
				}
		    }
		}


		# Draw the marginal precision of the right eigenfunctions (Gibbs step)
		draw.lambdav[i,]=d.lambdav(draw.v[[i-1]],draw.theta[i-1,],draw.rhos[[i]],pinfo)


		# Draw the calibration parameters
		draw.theta[i,]=draw.theta[i-1,]
		draw.v[[i]]=draw.v[[i-1]]
		phi[[i]][,m+1]=phi[[i-1]][,m+1]
		for(j in 1:q)
		{
			theta.new=draw.theta[i,]
			theta.new[j]=runif(1,draw.theta[i,j]-rth[j],draw.theta[i,j]+rth[j])
			v.new=draw.v[[i]]
			a=-Inf
			if(theta.new[j]>0 && theta.new[j]<1) {
				a=0
				for(l in 1:nc) {
					v.new[m+1,l]=gp.crealize2(set.lv(theta.new,pinfo$l.v),draw.v[[i]][1:m,l],0,draw.lambdav[i,l],draw.rhos[[i]][l,])
				}
				phi0.new=draw.u[[i]]%*%v.new[m+1,]
				a=a+dlp.th(j,theta.new[j],draw.theta[i,j],pinfo)+dlp.yfmd(yf,phi[[i]][,m+1],phi0.new,pinfo$inidelta,draw.lambdaf[i-1,],draw.lambdad[i-1,],draw.psis[[i-1]],pinfo)
			}
			a=min(0,a)
			if(log(runif(1))<a)
			{
				draw.theta[i,j]=theta.new[j]
				draw.v[[i]]=v.new
				phi[[i]][,m+1]=phi0.new
				accept.theta[j]=accept.theta[j]+1
			}
		}


		# Draw delta|rest
		draw.delta[i,]=pred.delta(phi[[i]][,m+1],yf,draw.lambdad[i-1,],draw.psis[[i-1]],draw.lambdaf[i-1,],draw.kappa[i-1,],pinfo)


		# Draw lambda_delta|rest
		draw.lambdad[i,]=d.lambdad(draw.delta[i,],draw.psis[[i-1]],pinfo)


		# Draw psis|rest for k=1..p, these are the discrepancy correlation parameters
		for(l in 1:ns)
		{
			draw.psis[[i]][l,]=draw.psis[[i-1]][l,]
			for(j in 1:p)
			{
				psis.new=draw.psis[[i]][l,]
				psis.new[j]=runif(1,draw.psis[[i]][l,j]-rr[[l]][j],draw.psis[[i]][l,j]+rr[[l]][j])
				a=-Inf
				if(psis.new[j]>0 && psis.new[j]<1)
					a=dlp.psi(l,j,draw.delta[i,],draw.lambdad[i,],draw.psis[[i]][l,],psis.new,pinfo)
				a=min(0,a)
				if(log(runif(1))<a)
				{
					draw.psis[[i]][l,j]=psis.new[j]
					accept.psis[[l]][j]=accept.psis[[l]][j]+1
				}
			}
		}

		# Draw kappa|rest
		for(l in 1:ns)
			draw.kappa[i,l]=d.kappa(l,phi[[i]][,m+1],yf,draw.delta[i,],draw.lambdaf[i-1,],pinfo)

		# Draw lambda_f|rest
		for(j in 1:pinfo$ns)
			draw.lambdaf[i,j]=d.lambdaf(j,phi[[i]][,m+1],yf,draw.delta[i,],draw.kappa[i,],pinfo)


		# Adapt and summary
		cat(i/N*100," percent complete     \r")
		if(i%%(250)==0 && i<(N*.5+1))
		{
			rate.rhos=accept.rhos
			for(l in 1:nc) {
				rate.rhos[[l]]=accept.rhos[[l]]/(i-lastadapt)
			}
		   	cat("\n\nAdapting rates of rhos from ")
			for(l in 1:nc) {
				cat("\n",rate.rhos[[l]])
		    	for(j in 1:q)
					if(rate.rhos[[l]][j]>.49 || rate.rhos[[l]][j]<.39) rp[[l]][j]=rp[[l]][j]*rate.rhos[[l]][j]/.44
			}

		  	rate.theta=accept.theta/(i-lastadapt)
		  	cat("\nAdapting rates of theta from \n",rate.theta,"\n")
			for(j in 1:q)
				if(rate.theta[j]>.49 || rate.theta[j]<.39) rth[j]=rth[j]*rate.theta[j]/.44
			rth[rth<0.001]=0.001

			rate.psis=accept.psis
			for(l in 1:ns) {
				rate.psis[[l]]=accept.psis[[l]]/(i-lastadapt)
			}
		  	cat("Adapting rates of psis from ")
		  	for(l in 1:ns) {
		  		cat("\n",rate.psis[[l]])
			  	for(j in 1:p)
			  		if(rate.psis[[l]][j]>.49 || rate.psis[[l]][j]<.39) rr[[l]][j]=rr[[l]][j]*rate.psis[[l]][j]/.44
			}
			cat("\n")
			lastadapt=i
			for(l in 1:nc) accept.rhos[[l]]=rep(1,q)
			for(l in 1:ns) accept.psis[[l]]=rep(1,p)
			accept.theta=rep(1,q)
		}

	} # end of MCMC loop

	cat("proposal widths:\n")
	cat(" rth:",rth,"\n")
	cat("  rp:",rp[[1]][1], " ; ",rp[[2]][1],"\n")
	cat("  rr:",rr[[1]][1],"\n")

	rx=(N-last):N  # which ones we'll return as after burn-in

	return(list(theta=as.matrix(draw.theta[rx,]),	rhos=draw.rhos[rx],		lambdav=draw.lambdav[rx,],
				v=draw.v[rx],		u=draw.u[rx],	psis=draw.psis[rx],
				lambdad=draw.lambdad[rx,],	phi=phi[rx], delta=draw.delta[rx,],	kappa=draw.kappa[rx,],
				lambdaf=draw.lambdaf[rx,],   yf=yf,	N=N,	pinfo=pinfo, 	mh=mh,	last=last
				))
}
