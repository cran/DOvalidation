\name{K.epa}
\alias{K.epa}
\title{
Epanechnikov Kernel
}
\description{
Evaluation of the Epanechnikov kernel function}
\usage{
K.epa(u)
}
\arguments{
  \item{u}{
A vector (or scalar) with the evaluation point(s). }
}
\value{
The value of the kernel function at \code{u}.}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2016). Double one-sided cross-validation of local linear hazards. Journal of the Royal Statistical Society B, 78, 755-779.
}
\author{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.
}

\seealso{
\code{\link{K.sextic}}
}
\examples{
curve(K.epa,-1.5,1.5,main="Epanechnikov kernel",ylab="K(u)",xlab="u")
# The left onesided 
K.epa.left<-function(u) return(2*K.epa(u)*(u<0))
curve(K.epa.left,-1.5,1.5,main="Left onesided Epanechnikov kernel",ylab="K(u)",xlab="u")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
