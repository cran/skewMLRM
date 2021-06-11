\name{matrix.sqrt}
\alias{matrix.sqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Square root of a matrix
}
\description{
Compute the square root of a matrix
}
\usage{
matrix.sqrt(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
a symmetric semi-definite positive matrix
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A symmetric matrix, say B, such as B^t*B=A
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
A<-matrix(c(1,2,2,5),nrow=2)
B<-matrix.sqrt(A)
##Recovering A
t(B)\%*\%B
A
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
