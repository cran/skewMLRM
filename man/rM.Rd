\name{rM}
\alias{rMCN}
\alias{rMN}
\alias{rMT}
\alias{rMSL}
\alias{rMSN}
\alias{rMSTN}
\alias{rMSSL}
\alias{rMSCN}
\alias{rMSTT}
\alias{rMSSL2}
\alias{rMSCN2}
\alias{rMSNC}
\alias{rMSTEC}
\alias{rMSSLEC}
\alias{rMSCEC}
\title{
Random generation for models in the MSMN, MSMSN, MSSMN and MSMSNC classes
}
\description{
rxxx generates random values for the distribution xxx, where xxx is any supported
model in the multivariate scale mixtures of normal (MSMN),
multivariate scale mixtures of skew-normal (MSMSN), multivariate skew scale mixtures of normal (MSSMN) or
multivariate scale mixtures of skew-normal-Cauchy (MSMSNC) classes. See details for supported distributions.
}
\usage{
rMN(n, mu, Sigma)
rMT(n, mu, Sigma, nu = 1)
rMSL(n, mu, Sigma, nu = 1)
rMCN(n, mu, Sigma, nu = 0.5, gamma = 0.5)
rMSN(n, mu, Sigma, lambda)
rMSTN(n, mu, Sigma, lambda, nu = 1)
rMSSL(n, mu, Sigma, lambda, nu = 1)
rMSCN(n, mu, Sigma, lambda, nu = 0.5, gamma = 0.5)
rMSTT(n, mu, Sigma, lambda, nu = 1)
rMSSL2(n, mu, Sigma, lambda, nu = 1)
rMSCN2(n, mu, Sigma, lambda, nu = 0.5, gamma = 0.5)
rMSNC(n, mu, Sigma, lambda)
rMSTEC(n, mu, Sigma, lambda, nu = 1)
rMSSLEC(n, mu, Sigma, lambda, nu = 1)
rMSCEC(n, mu, Sigma, lambda, nu = 0.5, gamma = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of observations to be generated.
}
  \item{mu}{
vector of location parameters. 
}
  \item{Sigma}{
covariance matrix (a positive definite matrix).
}
\item{lambda}{
vector of shape parameters.
}
\item{nu}{
nu parameter. A positive scalar for MT, MSL, MSTN, MSSL, MSTT, MSSL2, MSTEC and MSSLEC models. A value in the interval (0,1) for MCN, MSCN, MSCN2 and MSCEC models.
}
\item{gamma}{
gamma parameter. A value in the interval (0,1) for MCN, MSCN, MSCN2 and MSCEC models.
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

MN used \code{\link{mvrnorm}}. For MT, MSL and MCN, the generation is based on the MSMN class. See Lange and Sinsheimer (1993) for details. 
For MSTN, MSSL and MSCN, the generation is based on the MSSMN class. See Ferreira, Lachos and Bolfarine (2016) for details.
For MSTT, MSSL2 and MSCN2, the generation is based on the multivariate scale mixtures of skew-normal (MSMSN) class.
See Branco and Dey (2001) for details.
For MSNC, the generation is based on the stochastic representation in Proposition 2.1 of
Kahrari et al. (2016). For the MSTEC, MSSLEC and MSCEC models, the generation is based on the 
MSMSNC class. See Kahrari et al. (2017) for details.
}
\value{
A matrix with the generated data.
}
%\seealso{
%\code{\link{estimate.MN}}, \code{\link{estimate.MT}}, \code{\link{estimate.MSL}}, \code{\link{estimate.MCN}}, 
%\code{\link{estimate.MSN}}, \code{\link{estimate.MSTN}}, \code{\link{estimate.MSSL}}, \code{\link{estimate.MSCN}},
%\code{\link{estimate.MSTT}}, \code{\link{estimate.MSSL2}}, \code{\link{estimate.MSCN2}},
%\code{\link{estimate.MSNC}}, \code{\link{estimate.MSTEC}}, \code{\link{estimate.MSSLEC}}, \code{\link{estimate.MSCEC}}
%}
\references{
Branco, M.D., Dey, D.K. (2001). A general class of multivariate skew-elliptical distributions. Journal
of Multivariate Analysis 79, 99-113.

Ferreira, C.S., Lachos, V.H., Bolfarine, H. (2016). Likelihood-based inference for multivariate skew
scale mixtures of normal distributions. AStA Advances in Statistical Analysis 100, 421-441.

Kahrari, F., Rezaei, M., Yousefzadeh, F., Arellano-Valle, R.B. (2016). On the multivariate skew-normal-Cauchy distribution.
Statistics and Probability Letters, 117, 80-88. 

Kahrari, F., Arellano-Valle, R.B., Rezaei, M., Yousefzadeh, F. (2017). Scale mixtures of skew-normal-Cauchy distributions.
Statistics and Probability Letters, 126, 1-6.

Lange, K., Sinsheimer, J.S. (1993). Normal/independent distributions and their applications in robust
regression. Journal of Computational and Graphical Statistics 2, 175-198.
}
\author{
Clecio Ferreira, Diego Gallardo and Camila Zeller.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
rMSN(10, mu=c(0,0), Sigma=diag(2), lambda=c(1,-1)) ##bivariate MSN model
rMSNC(10, mu=0, Sigma=2, lambda=1) ##univariate MSNC model
rMSNC(10, mu=1:3, Sigma=2*diag(3), lambda=c(1,-1,0)) ##trivariate MSN model
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
%\keyword{ ~Data generation }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.