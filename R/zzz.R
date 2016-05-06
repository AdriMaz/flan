.onLoad <- function(libname, pkgname) {
	loadRcppModules()
# 	local({
# 		.flantol=.Machine$double.eps^0.5
# 		.integrands=list(
# 				CLONE_P0_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)},
# 				CLONE_PK_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)},
# # 				(1-x)^(k-1)*x^rho/(1-delta*x/k)^(k+1)
# 				CLONE_dP0_dr_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)*log(x)},
# 				CLONE_dPK_dr_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)*log(x/k)}
# 				)
# 				},globalenv())
}

.onAttach <- function(libname, pkgname) {
	local({
		.flantol=.Machine$double.eps^0.5
		.flansubd=1e3L
		.tunings=list(z1=0.1,z2=0.9,z3=0.8,z4=0.5,q=0.1)
		.integrands=list(
				CLONE_PGF=function(x,rho,delta) {x^rho/(1+x*delta)},
				CLONE_PGF_dr=function(x,rho,delta) {x^rho/(1+x*delta)*log(x)},
				CLONE_P0_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)},
# 				CLONE_PK_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)},
				CLONE_PK_WD=function(x,rho,delta,k) {x^rho*(1-x)^(k-1)/(1-x*delta)^(k+1)},
				CLONE_dP0_dr_WD=function(x,rho,delta) {(1-x)*x^(rho-1)/(1-delta*x)*log(x)},
# 				CLONE_dPK_dr_WD=function(x,rho,delta,k) {x^rho*(1-x/k)^(k-1)/(1-x/k*delta)^(k+1)*log(x/k)}
				CLONE_dPK_dr_WD=function(x,rho,delta,k) {x^rho*(1-x)^(k-1)/(1-x*delta)^(k+1)*log(x)}

				
# 				LAP_LOG_NORMAL=function(s,a,l,m){               # Laplace transform
# 						  f <- function(x){exp(-(s*exp(a*x*sqrt(2)+l)+x^2))}
# 
# 						  I <- {integrate(f,lower=-Inf,upper=Inf,rel.tol=.flantol)}      # integrate from - to + infinity
# 						  
# 						  I <- I$value                      # get the value
# 						  I <- I/sqrt(pi)			# rescale
# 						  return(I-1/m)                     # return shifted Laplace transform
# 						}
				)
				},globalenv())
}
