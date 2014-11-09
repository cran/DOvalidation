\name{DOvalidation-package}
\alias{DOvalidation-package}
\alias{DOvalidation}
\docType{package}
\title{
Local Linear Hazard Estimation with Do-Validated and Cross-Validated Bandwidths
}
\description{
Local linear estimator for the univariate hazard (hazard rate) and bandwidth parameter selection using the do-validation method and the standard least squares cross-validation method.
}
\details{
\tabular{ll}{
Package: \tab DOvalidation\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2014-11-09\cr
License: \tab GPL-2\cr
}
}
\author{
M.L. Gamiz, E. Mammen, M.D. Martinez-Miranda and J.P. Nielsen

Maintainer: Maria Dolores Martinez-Miranda <mmiranda@ugr.es>
}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2014). Do-Validating Local Linear Hazards. Available at SSRN: http://dx.doi.org/10.2139/ssrn.2504497

Nielsen, J.P. and Tanggaard, C. (2001). Boundary and bias correction in kernel hazard estimation. Scandinavian Journal of Statistics,28, 675-698.
}
\keyword{ package }
\keyword{ nonparametric }
\examples{
## Calculation of the local linear hazard estimator with do-validated bandwidth.
## This example is described in Gamiz et al. (2014)
data(UK)
Oi<-UK$D
Ei<-UK$E
ti<-40:110  # time is age and it goes from 40 to 110 years
M<-length(ti)

## First we estimate the bandwidth using do-validation
## (mean of left- and right- onesided cross-validation)
my.bs<-seq(1.5,40,length=100)
## The left onesided cross-validated bandwidth
res.left<-b.OSCV(grid.b=my.bs,K="sextic",Ktype="left",xi=ti,Oi=Oi,Ei=Ei)
bleft<-res.left$boscv
## The right onesided cross-validated bandwidth
res.right<-b.OSCV(grid.b=my.bs,K="sextic",Ktype="right",xi=ti,Oi=Oi,Ei=Ei) 
bright<-res.right$boscv
## The do-validated bandwidth
bdo<-(bleft+bright)/2

## Now we calculate the local linear hazard with the do-validated bandwidth
resLL.do<-hazard.LL(ti,Oi,Ei,x=ti,b=bdo,K="sextic",Ktype="symmetric")

## The local linear hazard estimate is hLL.do below
hLL.do<-resLL.do$hLL
plot(ti,hLL.do,main="Local linear hazard estimator",xlab='age',ylab='',
type='l',lwd=2)
}
