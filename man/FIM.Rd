\name{FIM}
\alias{FI.MCN}
\alias{FI.MN}
\alias{FI.MT}
\alias{FI.MSL}
\alias{FI.MSN}
\alias{FI.MSTN}
\alias{FI.MSSL}
\alias{FI.MSCN}
\alias{FI.MSTT}
\alias{FI.MSSL2}
\alias{FI.MSCN2}
\alias{FI.MSNC}
\alias{FI.MSTEC}
\alias{FI.MSSLEC}
\alias{FI.MSCEC}
\title{
Observed Fisher information matrix for distributions in the MSMN, MSMSN, MSSMN and MSMSNC classes.
}
\description{
FI.xxx computes the observed Fisher information (FI) matrix for the distribution xxx, where xxx is any supported
model in the multivariate scale mixtures of normal (MSMN),
multivariate scale mixtures of skew-normal (MSMSN), multivariate skew scale mixtures of normal (MSSMN) or
multivariate scale mixtures of skew-normal-Cauchy (MSMSNC) classes. See details for supported distributions.
}
\usage{
FI.MN(P, y, X)
FI.MT(P, y, X)
FI.MSL(P, y, X)
FI.MCN(P, y, X)
FI.MSN(P, y, X)
FI.MSTN(P, y, X)
FI.MSSL(P, y, X)
FI.MSCN(P, y, X)
FI.MSTT(P, y, X, nu)
FI.MSSL2(P, y, X, nu)
FI.MSCN2(P, y, X, nu, gamma)
FI.MSNC(P, y, X)
FI.MSTEC(P, y, X, nu)
FI.MSSLEC(P, y, X, nu)
FI.MSCEC(P, y, X, nu, gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
the estimated parameters returned by a function of the form estimate.xxx, where xxx is a supported distribution.
}
  \item{y}{
The multivariate vector of responses. The univariate case also is supported.
}
  \item{X}{
The regressor matrix.
}
\item{nu}{
nu parameter. Only for MSTT, MSSL2, MSTEC, MSSLEC and MSCEC distributions.
}
\item{gamma}{
gamma parameter. Only for MSCN2 and MSCEC distributions.
}
}
\details{
Supported models are: 

In MSMN class: multivariate normal (MN), multivariate Student t (MT), multivariate slash (MSL),
multivariate contaminated normal (MCN). See Lange and Sinsheimer (1993) for details.

In MSMSN class: multivariate skew-normal (MSN), multivariate skew-T (MSTT), multivariate skew-slash (MSSL2),
multivariate skew-contaminated normal (MSCN2). See Zeller, Lachos and Vilca-Labra (2011) for details.

In MSSMN class: MSN, multivariate skew-t-normal (MSTN), multivariate skew-slash normal 
(MSSL), multivariate skew-contaminated normal (MSCN). See Louredo, Zeller and Ferreira (2021) for details.

In MSMSNC class: multivariate skew-normal-Cauchy (MSNC), multivariate 
skew-t-Expected-Cauchy (MSTEC),  multivariate skew-slash-Expected-Cauchy 
(MSSLEC), multivariate 
skew-contaminated-Expected-Cauchy (MSCEC). See Kahrari et al. (2020) for details.

Note: the MSN distribution belongs to both, MSMSN and MSSMN classes.
}
\value{
A matrix with the observed FI matrix for the specified model.
}
\references{
Kahrari, F., Arellano-Valle, R.B., Ferreira, C.S., Gallardo, D.I. (2020) Some Simulation/computation
in multivariate linear models of scale mixtures of skew-normal-Cauchy distributions. Communications in Statistics - 
Simulation and Computation. In press. DOI: 10.1080/03610918.2020.1804582

Lange, K., Sinsheimer, J.S. (1993). Normal/independent distributions and their applications in robust
regression. Journal of Computational and Graphical Statistics 2, 175-198.

Louredo, G.M.S., Zeller, C.B., Ferreira, C.S. (2021). Estimation and influence diagnostics for the
multivariate linear regression models with skew scale mixtures of normal distributions. Sankhya B.
In press. DOI: 10.1007/s13571-021-00257-y

Zeller, C.B., Lachos, V.H., Vilca-Labra, F.E. (2011). Local influence analysis for regression models 
with scale mixtures of skew-normal distributions. Journal of Applied Statistics 38, 343-368.
}
\author{
Clecio Ferreira, Diego Gallardo and Camila Zeller.
}
\note{
For MSTEC and MSSLEC and distributions, nu>0 is considered as fixed. For MSCEC distribution, 0<nu<1 and 0<gamma<1 are
considered as fixed.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{estimate.MN}}, \code{\link{estimate.MT}}, \code{\link{estimate.MSL}}, \code{\link{estimate.MCN}}, 
%\code{\link{estimate.MSN}}, \code{\link{estimate.MSTN}}, \code{\link{estimate.MSSL}}, \code{\link{estimate.MSCN}},
%\code{\link{estimate.MSTT}}, \code{\link{estimate.MSSL2}}, \code{\link{estimate.MSCN2}},
%\code{\link{estimate.MSNC}}, \code{\link{estimate.MSTEC}}, \code{\link{estimate.MSSLEC}}, \code{\link{estimate.MSCEC}}, 
%\code{\link{choose.MSMN}}, \code{\link{choose.MSMSN}}, \code{\link{choose.MSMSNC}}, \code{\link{choose.MSSMN}},
%\code{\link{choose2}}
%}
\examples{
set.seed(2020)
n=200   # length of the sample
nv<-3   # number of explanatory variables
p<-nv+1 # nv + intercept
m<-4    # dimension of Y
q0=p*m
X<-array(0,c(q0,m,n)) 
for(i in 1:n) {
    aux=rep(1,p)
    aux[2:p]<-rMN(1,mu=rnorm(nv),Sigma=diag(nv))
    mi=matrix(0,q0,m)
    for (j in 1:m) mi[((j-1)*p+1):(j*p),j]=aux
    X[,,i]<-mi
} #Simulated matrix covariates
betas<-matrix(rnorm(q0),ncol=1) ## True betas
Sigmas <- clusterGeneration::genPositiveDefMat(m,rangeVar=c(1,3), 
lambdaLow=1, ratioLambda=3)$Sigma ##True Sigma
lambda<-rnorm(m) ##True lambda
y=matrix(0,n,m)
for(i in 1:n) {
     mui<-t(X[,,i])\%*\%betas
     y[i,]<-rMSN(n=1,c(mui),Sigmas,lambda)}
\donttest{
fit.MSN=estimate.MSN(y,X) ##Estimate parameters for MSN model
fit.MSN  ## Output of estimate.MSN
summary(fit.MSN)
fit.MSN$se ##Estimated standard errors by the estimate.MSN function 
##Estimated standard errors by minus the square root of 
##the diagonal from the observed FI matrix of the MSN model 
sqrt(diag(solve(-FI.MSN(fit.MSN$coefficients, y, X))))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~model }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
