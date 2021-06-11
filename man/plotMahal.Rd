\name{plotMahal}
\alias{plotMahal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mahalanobis distance for fitted models in the MSMN, MSMSN, MSSMN and MSMSNC classes
}
\description{
Compute and plot the Mahalanobis distance for any supported model in the multivariate scale mixtures of normal (MSMN),
multivariate scale mixtures of skew-normal (MSMSN), multivariate skew scale mixtures of normal (MSSMN) or
multivariate scale mixtures of skew-normal-Cauchy (MSMSNC) classes. See details for supported distributions.
}
\usage{
plotMahal(fit, y, X, dist = "MN", alpha = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{an object returned by the function estimate.xxx, where xxx is any supported model. See details for supported distributions.}
  \item{y}{The multivariate vector of responses. The univariate case also is supported.}
  \item{X}{the regressor matrix.}
  \item{dist}{the distribution used to obtain the fit object. See details for supported distributions.}
  \item{alpha}{significance level (0.05 by default).}
  \item{...}{... can include any of the other arguments of plots}
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
A plot with the Mahalanobis distance for all the observations and a cut-off to detect possible influent observations
based on the specified significance (0.05 by default). 
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
Clecio Ferreira, Diego Gallardo and Camila Zeller
}
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
    aux[2:p]<-rMN(1,mu=rnorm(nv),Sigma=diag(nv)) ##simulating covariates
    mi=matrix(0,q0,m)
    for (j in 1:m) mi[((j-1)*p+1):(j*p),j]=aux
    X[,,i]<-mi
} ##X is the simulated regressor matrix
betas<-matrix(rnorm(q0),ncol=1) ##True betas
Sigmas <- clusterGeneration::genPositiveDefMat(m,rangeVar=c(1,3), 
lambdaLow=1, ratioLambda=3)$Sigma ##True Sigma
y=matrix(0,n,m)
for(i in 1:n) {
     mui<-t(X[,,i])\%*\%betas
     y[i,]<-rMN(n=1,c(mui),Sigmas) ## simulating the response vector 
}
theta.MN=estimate.MN(y,X) #fit the MN model
plotMahal(theta.MN, y, X) #Plot the Mahalanobis distance for MN model
theta.MT=estimate.MT(y,X) #fit the MT model
plotMahal(theta.MT, y, X) #Plot the Mahalanobis distance for MT model
}
