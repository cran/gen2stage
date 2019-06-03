#library(clinfun)

# single-stage design
gen2single <- function(pu,pa,ep1,ep2,nsoln=5) {
	if(pu<pa){
		soln <- ph2single(pu=pu,pa=pa,ep1=ep1,ep2=ep2,nsoln=nsoln)
	}else if(pu>pa){
		tmp.soln <- ph2single(pu=(1-pu),pa=(1-pa),ep1=ep1,ep2=ep2,nsoln=nsoln)
		soln <- tmp.soln
		soln$r <- tmp.soln$n-tmp.soln$r
	}else{
		stop("!!!Should pu and pa be either pu>pa or pu<pa!!!\n")
	}
	soln

	ph1 <- list()
  	ph1$pu <- pu
  	ph1$pa <- pa
  	ph1$alpha <- ep1
  	ph1$beta <- ep2
  	ph1$out <- soln
  	ph1$nsoln <- nsoln
  	class(ph1) <- "gen2single"
	ph1
}

print.gen2single <- function(x, ...) {
  	xout <- x$out
  	nmax <- x$nsoln
  	xopt <- xout[1,]
  	dimnames(xopt)[[1]] <- c("Optimal")
  	cat("\n Generalized 1-stage Phase II design \n\n")
  	cat("Unacceptable response/toxicity rate: ",x$pu,"\n")
  	cat("Desirable response/toxicity rate: ",x$pa,"\n")
  	cat("Error rates: alpha = ",x$alpha,"; beta = ",x$beta,"\n\n")
 	print(xopt, digits = 4, ...)
  	cat("\n")
}

# two-stage design
gen2simon <- function(pu, pa, ep1, ep2, nmax = 100) {
  	if(nmax > 1000) stop("nmax cannot exceed 1000")

	alpha.beta <- function(data){
		oc.gentwostage.bdry(data[1],data[2],data[3],data[4],data[5],data[6])
	}

	if(pu<pa){
		gph2 <- ph2simon(pu=pu,pa=pa,ep1=ep1,ep2=ep2,nmax=nmax)
		tout <- gph2$out
		t2out <- cbind(rep(pu,dim(tout)[1]),rep(pa,dim(tout)[1]),tout[,c(1:4)])
		t3out <- apply(t2out,1,alpha.beta)
		tmp.alpha <- as.numeric(t3out[1,])
		tmp.beta <- 1-as.numeric(t3out[2,])
		t4out <- tout
		t5out <- cbind(t4out,alpha=tmp.alpha,beta=tmp.beta)
		gph2$out <- t5out
		#print(class(t5out))
		#print(names(t5out))
		#print(head(tmp.alpha))
		#print(head(t5out))
  	}else if(pu>pa){
		tmp.gph2 <- ph2simon(pu=(1-pu),pa=(1-pa),ep1=ep1,ep2=ep2,nmax=nmax)

		tout <- tmp.gph2$out
		t2out <- cbind(rep((1-pu),dim(tout)[1]),rep((1-pa),dim(tout)[1]),tout[,c(1:4)])
		t3out <- apply(t2out,1,alpha.beta)
		tmp.alpha <- as.numeric(t3out[1,])
		tmp.beta <- 1-as.numeric(t3out[2,])
		t4out <- tout
		t5out <- cbind(t4out,alpha=tmp.alpha,beta=tmp.beta)
		tmp.gph2$out <- t5out
		#print(class(t5out))
		#print(names(t5out))
		#print(head(tmp.alpha))
		#print(head(t5out))

		gph2 <- tmp.gph2
		gph2$pu <- pu
		gph2$pa <- pa
		gph2$out[,1] <- tmp.gph2$out[,2]-tmp.gph2$out[,1]
		gph2$out[,3] <- tmp.gph2$out[,4]-tmp.gph2$out[,3]
	}else{
		stop("!!!Should pu and pa be either pu>pa or pu<pa!!!\n")
	}
	ph2 <- gph2
	class(ph2) <- "gen2simon"
  	ph2
}

print.gen2simon <- function(x, ...) {
  	xout <- x$out
  	nmax <- x$nmax
  	n <- nrow(xout)
  	nopt <- ((1:n)[xout[,5]==min(xout[,5])])[1]
  	xopt <- xout[c(nopt,1),]
  	dimnames(xopt)[[1]] <- c("Optimal","Minimax")
  	cat("\n Generalized 2-stage Phase II design \n\n")
  	cat("Unacceptable response/toxicity rate: ",x$pu,"\n")
  	cat("Desirable response/toxicity rate: ",x$pa,"\n")
  	cat("Error rates: alpha = ",x$alpha,"; beta = ",x$beta,"\n\n")
  	print(xopt, digits = 4, ...)
  	cat("\n")
  	if(xopt[1,4]>nmax-10) warning(paste("  Optimal sample size too close to nmax. \n  Try increasing nmax (current value = ",nmax,")\n",sep=""))
}

plot.gen2simon <- function(x, ...) {
  	xout <- x$out
  	n <- nrow(xout)
  	nopt <- ((1:n)[xout[,5]==min(xout[,5])])[1]
  	nopt1 <- min(nopt+5,n)
  	plot(xout[1:nopt1,4],xout[1:nopt1,5],type="l",xlab="Maximum number of patients",ylab="Expected trial size", ...)
  	points(xout[1,4],xout[1,5],pch="M")
  	points(xout[nopt,4],xout[nopt,5],pch="O")
}

# returns the type I and II error rates as well as the probability of early temination
# and expected sample size under pu for a specific boundary
oc.gentwostage.bdry <- function(pu, pa, r1, n1, r, n){
	if(pu<pa){
		out <- oc.twostage.bdry(pu=pu,pa=pa,r1=r1,n1=n1,r=r,n=n)
	}else if(pu>pa){
		tmp.out <- oc.twostage.bdry(pu=(1-pu),pa=(1-pa),r1=(n1-r1),n1=n1,r=(n-r),n=n)
		out <- tmp.out
	}else{
		stop("!!!Should pu and pa be either pu>pa or pu<pa!!!\n")
	}
	out
}

