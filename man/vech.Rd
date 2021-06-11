\name{vech}
\alias{vech}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vectorize a symmetric matrix
}
\description{
vech takes the upper diagonal from a symmetric matrix and vectorizes it.
}
\usage{
vech(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a symmetric matrix.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector with the components of the upper diagonal from the matrix, listed by row.
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
%\code{\link{xpnd}}
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
