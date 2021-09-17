\name{solve2}
\alias{solve2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the inverse of a matrix
}
\description{
Computes the inverse of a matrix using the LU decomposition.
}
\usage{
solve2(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
an invertible square matrix.
}
}
\details{
Use the LU decomposition to compute the inverse of a matrix. In some cases, solve produces error to invert a matrix whereas
this decomposition provide a valid solution.
}
\value{
A square matrix corresponding to the inverse of A.
}
\references{
Bellman, R. (1987). Matrix Analysis, Second edition, Classics in Applied Mathematics, Society for Industrial and Applied Mathematics.

Horn, R. A. and C. R. Johnson (1985). Matrix Analysis, Cambridge University Press.
}
\author{
Clecio Ferreira, Diego Gallardo and Camila Zeller 
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
A=matrix(c(1,2,5,6),ncol=2)
solve2(A)
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
