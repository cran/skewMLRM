\name{estimateM}
\alias{estimate.MCN}
\alias{estimate.MN}
\alias{estimate.MT}
\alias{estimate.MSL}
\alias{estimate.MSN}
\alias{estimate.MSTN}
\alias{estimate.MSSL}
\alias{estimate.MSCN}
\alias{estimate.MSTT}
\alias{estimate.MSSL2}
\alias{estimate.MSCN2}
\alias{estimate.MSNC}
\alias{estimate.MSTEC}
\alias{estimate.MSSLEC}
\alias{estimate.MSCEC}
\title{
Fitting a model in the MSMN, MSMSN, MSSMN and MSMSNC classes
}
\description{
estimate.Mxxx computes the maximum likelihood estimates for the distribution xxx, where xxx is any supported
model in the multivariate scale mixtures of normal (MSMN),
multivariate scale mixtures of skew-normal (MSMSN), multivariate skew scale mixtures of normal (MSSMN) or
multivariate scale mixtures of skew-normal-Cauchy (MSMSNC) classes. See details for supported distributions.
}
\usage{
estimate.MN(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE)
estimate.MT(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
     nu.min = 2.0001)
estimate.MSL(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
     nu.min = 2.0001)
estimate.MCN(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE)
estimate.MSN(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE)
estimate.MSTN(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
     nu.min = 2.0001)
estimate.MSSL(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
     nu.min = 2.0001)
estimate.MSCN(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE)
estimate.MSTT(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
     nu.fixed = 3, nu.min = 2.0001)
estimate.MSSL2(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
     nu.fixed = 3, nu.min = 2.0001)
estimate.MSCN2(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
      nu.fixed = 0.1, gamma.fixed = 0.5)
estimate.MSNC(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE)
estimate.MSTEC(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
      nu.fixed = 3, nu.min = 2.0001)
estimate.MSSLEC(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
      nu.fixed = 3, nu.min = 2.0001)
estimate.MSCEC(y, X, max.iter = 1000, prec = 1e-04, est.var = TRUE, 
      nu.fixed = 0.1, gamma.fixed = 0.5)
}
\arguments{
  \item{y}{
The multivariate vector of responses. The univariate case also is supported.
}
  \item{X}{
The regressor matrix.
}
  \item{max.iter}{
The maximum number of iterations.
}
  \item{prec}{
The convergence tolerance for parameters.  
}
  \item{est.var}{
Logical. If TRUE the standard errors are estimated.
}
  \item{nu.fixed}{
If a numerical value is provided, the estimation consider nu as fixed. To estimate nu,
use nu.fixed=FALSE. Avaliable for MSTT, MSSL2, MSCN2, MSTEC, MSSLEC and MSCEC distributions.
For MSTT, MSSL2, MSTEC and MSSLEC, the default value is 3 and nu should be greater than 1. 
For MSCN2 and MSCEC, the default value is 0.1 and nu should be in the (0,1) interval.
}
  \item{gamma.fixed}{
If a numerical value is provided, the estimation consider gamma as fixed. To estimate gamma,
use gamma.fixed=FALSE. Avaliable for MSCN2 and MSCEC distributions.
For MSCN2 and MSCEC, the default value is 0.5 and gamma should be in the (0,1) interval.
}
\item{nu.min}{
Lower value to perform the maximization for nu. For MSTT, MSSL2, MSTEC and MSSLEC is used only when
nu.fixed=FALSE. 
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
an object of class "skewMLRM" is returned. The object returned for this functions is a list containing the following
components:

\item{coefficients}{A named vector of coefficients}
\item{se}{A named vector of the standard errors for the estimated coefficients. Valid if est.var is TRUE and 
the hessian matrix is invertible.}
\item{nu}{The estimated or fixed nu (only for MSTT, MSSL2, MSCN2, MSTEC, MSSLEC and MSCEC models)}
\item{gamma}{The estimated or fixed gamma (only for MSCN2 and MSCEC models)}
\item{logLik}{The log-likelihood function evaluated in the estimated parameters}
\item{AIC}{Akaike's Information Criterion}
\item{BIC}{Bayesian's Information Criterion}
\item{iterations}{the number of iterations until convergence (if attached)}
\item{time}{execution time in seconds}
\item{conv}{An integer code. 0 indicates successful completion. 1 otherwise.}
\item{dist}{The distribution for which was performed the estimation.}
\item{class}{The class for which was performed the estimation.}
\item{n}{The sample size}
\item{y}{The multivariate vector of responses. The univariate case also is supported.}
\item{X}{The regressor matrix (in a list form).}
\item{function}{a string with the name of the used function.}
%% ...
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

\note{
In MT, MSL, MSTN, MSSL, MSTT, MSSL2, MSTEC and MSSLEC distributions, nu>2 guarantees that the mean and variance exist,
nu>1 guarantees the existence only for the mean and for nu<=1, the mean and variance of the distribution is 
not finite.
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
data(ais, package="sn") ##Australian Institute of Sport data set
attach(ais)
##It is considered a bivariate regression model
##with Hg and SSF as response variables and
##Hc, Fe, Bfat and LBM as covariates
y<-cbind(Hg,SSF)
n<-nrow(y); m<-ncol(y)
X.aux=model.matrix(~Hc+Fe+Bfat+LBM)
p<-ncol(X.aux)
X<-array(0,dim=c(2*p,m,n))
for(i in 1:n) {
    X[1:p,1,i]=X.aux[i,,drop=FALSE]
    X[p+1:p,2,i]=X.aux[i,,drop=FALSE]
}
##See the covariate matrix X
##X
fit.MN=estimate.MN(y, X) ##Estimate the parameters for the MN regression model
summary(fit.MN)
fit.MT=estimate.MT(y, X) ##Estimate the parameters for the MT regression model
summary(fit.MT)
\donttest{
##may take some time on some systems
fit.MSSL=estimate.MSSL(y, X)   ##Estimate the parameters for the MSSL regression model
summary(fit.MSSL)
fit.MSTT=estimate.MSTT(y, X)   ##Estimate the parameters for the MSTT regression model
summary(fit.MSTT)
fit.MSNC=estimate.MSNC(y, X)   ##Estimate the parameters for the MSNC regression model
summary(fit.MSNC)
fit.MSCEC=estimate.MSCEC(y, X) ##Estimate the parameters for the MSCEC regression model
summary(fit.MSCEC)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
%\keyword{ ~model selection }
%\keyword{ ~covariates selection }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
