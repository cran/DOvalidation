\name{hazard.LL.RH}
\alias{hazard.LL.RH}
\title{Local Linear Hazard Estimator (Ramlau-Hansen Weighting)
}
\description{
Local linear estimator of the unidimensional hazard (or hazard rate) with Ramlau-Hansen weighting as was defined by Nielsen and Tanggaard (2001).}
\usage{
  hazard.LL.RH(xi , Oi , Ei , x , b , K="epa")
}
\arguments{
  \item{xi}{Vector of time points where the counts data are given.
}
  \item{Oi}{Vector with the number (counts) of occurrences observed at each time point (\code{xi}). 
}
  \item{Ei}{Vector with the observed exposure at each time point (\code{xi}).  
}
  \item{x}{Vector (or scalar) with the (time) grid points where the hazard estimator will be evaluated.
}  
  \item{b}{A positive scalar used as the bandwidth. 
}
  \item{K}{Indicates the kernel function to be considered in the estimator. Choose between values \code{"epa"} (for the Epanechnikov kernel) or \code{"sextic"} (see details for its expression).
}
}

\details{The estimator is calculated assuming that the data are given as count data i.e. number of occurrences and exposures.
The function allows to consider two different kernels using the argument \code{K}. These are: Epanechnikov, K(u)=.75*(1-u^2)*(abs(u)<1), and sextic K(u)=(3003/2048)*(1-(u)^2)^6)*(abs(u)<1).   
}

\value{
%hazard.res<-list(x=x,hLL=hLL)
  \item{x }{Vector (or scalar) with the (time) grid points where the hazard estimator has been evaluated. }
  \item{hLL }{Vector (or scalar) with the resulting hazard estimates at grid points x.}
}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2016). Double one-sided cross-validation of local linear hazards. Journal of the Royal Statistical Society B, 78, 755-779.

Nielsen, J.P. and Tanggaard, C. (2001). Boundary and bias correction in kernel hazard estimation. Scandinavian Journal of Statistics,28, 675-698.
}
\author{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hazard.LL}}
}
\examples{
## This example was described by Gamiz et al. (2016) to analyse the exposure robustness of
## local linear hazards with natural weigthing and Ramlau-Hansen weighting
data(Iceland)
Ei<-Iceland$E
Oi<-Iceland$D
xi<-40:110
n<-length(xi)
x<-seq(xi[1],xi[n],length=100)

## Hazard estimates with the original data
b0<-11.9899
alphaIC17<-hazard.LL.RH(xi,Oi,Ei,x,b=b0,K="sextic")$hLL
alLL17<-hazard.LL(xi,Oi,Ei,x,b=b0,K="sextic",Ktype="symmetric")$hLL
hi<-Oi/Ei;hi[Ei==0]<-0
print(round(hi[60:71],3))
## Hazard estimates with the modified data (one change in the exposure)
Ei2<-Ei; Ei2[67]<-2/365
alphaIC005<-hazard.LL.RH(xi,Oi,Ei2,x,b=b0,K="sextic")$hLL
alLL005<-hazard.LL(xi,Oi,Ei2,x,b=b0,K="sextic",Ktype="symmetric")$hLL

## Figure: Exposure robustness
old.par<-par(mfrow=c(2,2))
plot(x[73:100],alphaIC17[73:100],lwd=2,type='l',main='Exposure: 0.17',
xlab='',ylab='Ramlau-Hansen weighting')
plot(x[73:100],alphaIC005[73:100],lwd=2,type='l',main='Exposure: 0.005',
xlab='',ylab='Ramlau-Hansen weighting')
plot(x[73:100],alLL17[73:100],lwd=2,type='l',main='Exposure: 0.17',
xlab='',ylab='Natural  weighting')
plot(x[73:100],alLL005[73:100],lwd=2,type='l',main='Exposure: 0.005',
xlab='',ylab='Natural weighting')

par(old.par)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
