\name{mbackcrit}
\alias{mbackcrit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate backward based on the AIC or BIC criteria
}
\description{
mbackcrit implements the covariates selection based on backward and the Akaike's information criteria (AIC)
or Schwartz's information criteria (BIC) in a specified multivariate model in the multivariate scale mixtures of normal (MSMN),
multivariate scale mixtures of skew-normal (MSMSN), multivariate skew scale mixtures of normal (MSSMN) or
multivariate scale mixtures of skew-normal-Cauchy (MSMSNC) classes. 
See details for avaliable distributions.
}
\usage{
mbackcrit(y, X = NULL, max.iter = 1000, prec = 1e-04, dist = "MN", 
	criteria = "AIC", est.var=TRUE, cluster = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The multivariate vector of responses. The univariate case also is supported.
}
  \item{X}{
The regressor matrix. It should include intercept term for all the variates.
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
  \item{criteria}{
criteria used to perform the covariates selection. AIC (default) and BIC avaliable.
}
  \item{est.var}{
Logical. If TRUE the standard errors are estimated.
}
  \item{cluster}{
logical. If TRUE, parallel computing is used. FALSE is the default value.
}
\item{...}{
Possible aditional arguments. For instance, for MSTT, MSSL2, MSTEC and MSSLEC 
distributions should be added nu.min and nu.fixed related to specifications
for the nu parameter. 
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
\item{logLik}{The log-likelihood function evaluated in the estimated parameters for the selected model}
\item{AIC}{Akaike's Information Criterion for the selected model}
\item{BIC}{Bayesian's Information Criterion for the selected model}
\item{iterations}{the number of iterations until convergence (if attached)}
\item{conv}{An integer code for the selected model. 0 indicates successful completion. 1 otherwise.}
\item{dist}{The distribution for which was performed the estimation.}
\item{class}{The class for which was performed the estimation.}
\item{choose.crit}{the specified criteria to choose the distribution.}
\item{comment}{A comment indicating how many coefficients were eliminated}
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
\donttest{
##Perform covariates selection in the MN distribution
##based on the AIC criteria
##may take some time on some systems
fit.MN=mbackcrit(y, X, dist="MN")
summary(fit.MN)
##Identical process for MT distribution
fit.MT=mbackcrit(y, X, dist="MT")
summary(fit.MT)
##and for MSN distribution
fit.MSN=mbackcrit(y, X, dist="MSN")
summary(fit.MSN)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~model selection }
% \keyword{ ~covariates selection }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
