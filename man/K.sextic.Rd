\name{K.sextic}
\alias{K.sextic}
\title{Sextic Kernel}
\description{
Evaluation of the Sextic kernel function}
\usage{
K.sextic(u)
}
\arguments{
  \item{u}{
A vector (or scalar) with the evaluation point(s). }
}
\value{
The value of the kernel function at "u".}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2014). Do-Validating Local Linear Hazards. 
Available at SSRN: http://dx.doi.org/10.2139/ssrn.2504497
}
\author{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.
}

\seealso{
\code{\link{K.epa}}
}
\examples{
curve(K.sextic,-1.5,1.5,main="Sextic kernel",ylab="K(u)",xlab="u")
# The left onesided 
K.sextic.left<-function(u) return(2*K.sextic(u)*(u<0))
curve(K.sextic.left,-1.5,1.5,main="Left onesided sextic kernel",ylab="K(u)",xlab="u")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
