\name{xpnd}
\alias{xpnd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstitute a symmetric matrix from a vector.
}
\description{
xpnd reconstitutes a symmetric matrix from a vector obtained with the vech function.
}
\usage{
xpnd(x, nrow = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector with the components of the upper diagonal of the matrix
}
  \item{nrow}{
dimension of the matrix to be reconstitute.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A symmetric matrix.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Clecio Ferreira, Diego Gallardo and Camila Zeller.
}
\note{
For internal use.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{vech}}
%}
\examples{
A<-matrix(c(1,2,2,5),nrow=2)
##vectorized A matrix
B<-vech(A)
B
##reconstitute matrix A using B
xpnd(B,2)
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
