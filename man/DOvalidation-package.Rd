\name{DOvalidation-package}
\alias{DOvalidation-package}
\alias{DOvalidation}
\docType{package}
\title{
Kernel Hazard Estimation with Best One-Sided and Double One-Sided Cross-Validation
}
\description{
 Local linear hazard estimator and its multiplicatively bias correction, including three bandwidth selection methods: best one-sided cross-validation, double one-sided cross-validation, and standard cross-validation.
}
\details{
\tabular{ll}{
Package: \tab DOvalidation\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2017-10-20\cr
License: \tab GPL-2\cr
}
}
\author{
M.L. Gamiz, E. Mammen, M.D. Martinez-Miranda and J.P. Nielsen

Maintainer: Maria Dolores Martinez-Miranda <mmiranda@ugr.es>
}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2016). Double one-sided cross-validation of local linear hazards. Journal of the Royal Statistical Society B, 78, 755-779.

Gamiz, M.L., Martinez-Miranda, M.D. and Nielsen, J.P. (2017). Multiplicative local linear hazard estimation and best one-sided cross-validation. Available 	 at http://arxiv.org/abs/1710.05575  

Nielsen, J.P. and Tanggaard, C. (2001). Boundary and bias correction in kernel hazard estimation. Scandinavian Journal of Statistics,28, 675-698.
}
\keyword{ package }
\keyword{ nonparametric }
\keyword{survival}
\examples{
data(Iceland)
Oi<-Iceland$D
Ei<-Iceland$E
ti<-40:110  # time is age and it goes from 40 to 110 years
## Local linear hazard estimator 
## and its multiplicatively bias corrected version (MBC)
## with best one-sided cross-validated bandwidths
## Note: use functions b.BO and b.BO.MBC to get these bandwidths
##       (48.7 and 14.6, respectively)
res.LL<-hazard.LL(xi=ti,Oi=Oi,Ei=Ei,x=ti,b=14.6)
res.MBC<-hazard.MBC(xi=ti,Oi=Oi,Ei=Ei,x=ti,b=48.7)
plot(ti,res.LL$hLL,main='Hazard estimates',xlab='age',ylab='',
    type='l',col=4,lwd=2)
lines(ti,res.MBC$hMBC,col=2,lwd=2)
legend("topleft",bt="n",c("Local linear", "MBC"),col=c(4,2),lwd=2)
}