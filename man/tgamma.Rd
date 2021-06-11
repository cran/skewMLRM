\name{tgamma}
\alias{dtgamma}
\alias{qtgamma}
\title{
Truncated gamma distribution
}
\description{
Compute the probability density and quantile functions for the truncated gamma distribution with shape and scale parameters, 
restricted to the interval (a,b).
}
\usage{
dtgamma(x, shape, scale = 1, a = 0, b = Inf)
qtgamma(p, shape, scale = 1, a = 0, b = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles
}
\item{p}{
vector of probabilities
}
  \item{shape}{
shape parameter
}
  \item{scale}{
scale parameter
}
  \item{a}{
lower limit of range
}
  \item{b}{
upper limit of range
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
dtgamma gives the density function for the truncated gamma distribution. 
qtgamma gives the quantile function for the truncated gamma distribution.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Clecio Ferreira, Diego Gallardo and Camila Zeller
}
\note{
Auxiliary function to compute the E step for the Slash and Skew-slash models.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
##probability density and quantile function of the truncated gamma
##model with shape and scale parameters equal to 1
##evaluated in 2 and 0.75, respectively
dtgamma(2, shape=1, a=1)
qtgamma(0.75, shape=1, a=1)
##standard gamma distribution with shape parameter 2 evaluated in 1
dtgamma(1, shape=2)
dgamma(1, shape=2)
}