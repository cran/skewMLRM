\name{se.est}
\alias{se.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimated standard errors for the estimated parameters in distributions belonging to the MSMN, MSMSN, MSSMN and MSMSNC classes.
}
\description{
se.est computes the estimated standard errors based on the hessian matrix for supported
models in the multivariate scale mixtures of normal (MSMN),
multivariate scale mixtures of skew-normal (MSMSN), multivariate skew scale mixtures of normal (MSSMN) and
multivariate scale mixtures of skew-normal-Cauchy (MSMSNC) classes. See details for supported distributions.
}
\usage{
se.est(P, y, X, dist = "MN", nu = 0.5, gamma = 0.5)
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
  \item{dist}{
a supported multivariate distribution. See details.
}
\item{nu}{
nu parameter. Only for MSTT, MSSL2, MSCN2, MSTEC, MSSLEC and MSCEC distributions.
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
A vector with the estimated standard errors.
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
For MSTEC and MSSLEC and distributions, nu is considered as fixed. For MSCEC distribution, nu and gamma are
considered as fixed. This function is mainly for internal use.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
##See the simulated regressor matrix X
##X
aux<-estimate.MN(y,X) ##Estimate the parameters for the MN regression model
aux$estimate[,2] ##Estimated standard errors (se) by the estimate.MN function 
se.est(aux$estimate[,1], y, X, dist="MN") ##Estimated se for the MN model 
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
