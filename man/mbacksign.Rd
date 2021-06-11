\name{mbacksign}
\alias{mbacksign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Backward Based on Significance
}
\description{
mbacksign implements the covariates selection based on the significance of the covariates
in a specified multivariate model in the multivariate scale mixtures of normal (MSMN),
multivariate scale mixtures of skew-normal (MSMSN), multivariate skew scale mixtures of normal (MSSMN) or
multivariate scale mixtures of skew-normal-Cauchy (MSMSNC) classes. 
See details for avaliable distributions.
}
\usage{
mbacksign(y, X = NULL, max.iter = 1000, prec = 1e-04, dist = "MN", 
     significance = 0.05)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{dist}{
the multivariate distribution in which the covariates selection will be implemented.
}
  \item{significance}{
the level of significance to perform the covariate selection. By default is 0.05.
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
A list with the following components:
%%  If it is a LIST, use
\item{comment}{A comment indicating how many coefficients were eliminated.}
\item{estimate}{A matrix with the estimates and standard errors for the final model.}
\item{logLik}{The log-likelihood function evaluated in the estimated parameters}
\item{AIC}{Akaike's Information Criterion}
\item{BIC}{Bayesian's Information Criterion}
\item{conv}{An integer code. 0 indicates successful completion. 1 otherwise.}
\item{eliminated}{An string vector with the eliminated betas (in order of elimination).}
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
Clecio Ferreira, Diego Gallardo and Camila Zeller.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{mbacksign}}, \code{\link{choose.MSMN}}, \code{\link{choose.MSMSN}}, 
%\code{\link{choose.MSMSNC}}, \code{\link{choose.MSSMN}}, \code{\link{choose2}}
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
##See the regressor matrix X
##X
##Perform covariates selection in the MN distribution
##based on a significance level of 1%, 5% and 10% 
mbacksign(y, X, dist="MN", sign=0.01)
mbacksign(y, X, dist="MN", sign=0.05)
mbacksign(y, X, dist="MN", sign=0.10)
\donttest{
##may take some time on some systems
##identical process in the MCN model with 
##significance level of 5%
mbacksign(y, X, dist="MCN")
##for MSSL model
mbacksign(y, X, dist="MSSL")
##for MSNC model
mbacksign(y, X, dist="MSNC")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
