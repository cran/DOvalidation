\name{hazard.LL}
\alias{hazard.LL}
\title{
Local Linear Hazard Estimator (Natural Weighting)}
\description{
Local linear estimator of the unidimensional hazard (or hazard rate) with natural weighting introduced by Nielsen and Tanggaard (2001).}
\usage{
  hazard.LL(xi, Oi, Ei, x, b, K="epa", Ktype="symmetric" , CI=FALSE)
}
\arguments{
  \item{xi}{Vector of time points where the counts data are given.
}
  \item{Oi}{Vector with the number of occurrences observed at each time point (\code{xi}). 
}
  \item{Ei}{Vector with the observed exposure at each time point (\code{xi}).  
}
  \item{x}{Vector (or scalar) with the (time) grid points where the hazard estimator will be evaluated.
}  
  \item{b}{A positive scalar used as the bandwidth. 
}
  \item{K}{Indicates the kernel function to be considered in the estimator. Choose between values \code{"epa"} (for the Epanechnikov kernel) or \code{"sextic"} (see details for its expression).
}
  \item{Ktype}{Indicates the type of kernel to be used. Choose among \code{"symmetric"} for the usual kernel definition (chosen in the argument \code{K}), \code{"left"} for the left-sided version of the kernel, or \code{"right"} for the right-sided version. See details below.}
  \item{CI}{Logical. If \code{TRUE} then 95\% pointwise confidence intervals are provided for the hazard function.}
}

\details{The estimator is calculated assuming that the data are given as count data i.e. number of occurrences and exposures.
The function allows to consider two different kernels using the argument \code{K}. These are: Epanechnikov, K(u)=.75*(1-u^2)*(abs(u)<1), and sextic K(u)=(3003/2048)*(1-(u)^2)^6)*(abs(u)<1). The argument \code{Ktype} will define the usual estimator with whole support kernel as it is defined by \code{K} or the one-sided versions using left-sided kernel, 2*K(u)*(u<0), or right-sided kernel 2*K(u)*(u>0). See more details in Gamiz et al. (2016).  
}

\value{
%hazard.res<-list(x=x,OLL=OLL,ELL=ELL,hLL=hLL,OLL.norm=OLL.norm,ELL.norm=ELL.norm)
  \item{x }{Vector (or scalar) with the (time) grid points where the hazard estimator has been evaluated. }
  \item{OLL }{Vector with the smoothed occurrences (using the local linear kernel).}
  \item{ELL }{Vector with the smoothed exposures (using the local linear kernel).}
  \item{hLL }{Vector (or scalar) with the resulting hazard estimates at grid points \code{x}.}
  \item{OLL.norm }{Vector with the normalized smoothed occurrences (the smoothing weights are defined as for \code{O.LL} but normalized to add up one.}
  \item{ELL.norm }{Vector with the normalized smoothed exposures (the smoothing weights are defined as for \code{E.LL} but normalized to add up one.}
  \item{CI.inf}{Vector with the lower limits for the 95\% confidence intervals. If \code{CI=FALSE} then \code{NA} values are provided.}
    \item{CI.sup}{Vector with the upper limits for the 95\% confidence intervals. If \code{CI=FALSE} then \code{NA} values are provided.}
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
\code{\link{hazard.LL.RH}},\code{\link{hazard.MBC}}
}
\examples{
## Calculation of the local linear hazard estimator with do-validated bandwidth.
## The hazard estimator is shown and decomposed into smoothed occurrences and exposures.
## This example is described in Gamiz et al. (2016).
data(UK)
Oi<-UK$D
Ei<-UK$E
ti<-40:110  # time is age and it goes from 40 to 110 years
M<-length(ti)
country<-'UK'
bdo<-5.11
resLL.do<-hazard.LL(xi=ti,Oi=Oi,Ei=Ei,x=ti,b=bdo,K="sextic",Ktype="symmetric",CI=TRUE)

## The local linear hazard estimate is hLL.do below
hLL.do<-resLL.do$hLL

## The smoothed occurrences and exposures are:
ELL.norm.do<-resLL.do$ELL.norm
OLL.norm.do<-resLL.do$OLL.norm

## The 95\% pointwise confidence intervals based on the asymptotics are
hLL.do.inf<-resLL.do$CI.inf
hLL.do.sup<-resLL.do$CI.sup

# Now we plot the hazard estimator with confidence intervals
old.par<-par(mar=c(3,1.5,1.5,1.5),oma=c(2,0.5,0.5,0.2),
mgp=c(1.5,0.5,0),cex.axis=1,cex.main=1.5,mfrow=c(3,2))

#hazard estimate
tit<-paste(country,"Hazard estimate",sep= ' - ' )
yy<-range(c(hLL.do.inf,hLL.do.sup),na.rm=TRUE)
plot(ti,hLL.do,main=tit,xlab='age',ylab='',type='l',lwd=2,ylim=yy)
# the confidence bands
x1<-ti;x2<-ti[M:1]   
y1<-hLL.do.sup;y2<-hLL.do.inf[M:1]       
polygon(c(x1,x2,x1[1]),c(y1,y2,y1[1]),col=gray(0.7),border=FALSE)
lines(ti,hLL.do,lty=1,lwd=2,col=1)
  
## Zooming at the old mortality
ind.ages<- -c(1:60)  ## only women with ages 100 or higher
ti2<-ti[ind.ages];M2<-length(ti2)
yy2<-range(c(hLL.do.inf[ind.ages],hLL.do.sup[ind.ages]),na.rm=TRUE)
plot(ti2,hLL.do[ind.ages],main=tit,xlab='age',ylab='',type='l',
lwd=2,ylim=yy2)
# the confidence intervals
x1<-ti2;x2<-ti2[M2:1]   
y1<-hLL.do.sup[ind.ages];hLL.do.inf2<-hLL.do.inf[ind.ages]
y2<-hLL.do.inf2[M2:1]       
polygon(c(x1,x2,x1[1]),c(y1,y2,y1[1]),col=gray(0.7),border=FALSE)
lines(ti2,hLL.do[ind.ages],lty=1,lwd=2,col=1)
  
## We decompose the estimator in the smooth occurrences and exposures
#   The occurrences with a zoom at old-age mortality
yy<-range(OLL.norm.do,na.rm=TRUE)
plot(ti,OLL.norm.do,main="Smoothed occurrences",xlab='age',ylab='',type='l',
lwd=2,ylim=yy)
yy2<-range(OLL.norm.do[ind.ages],na.rm=TRUE)
plot(ti2,OLL.norm.do[ind.ages],main="Smoothed occurrences",xlab='age',
ylab='',type='l',lwd=2,ylim=yy2)
  
#   The exposures with a zoom at old-age mortality
yy<-range(ELL.norm.do,na.rm=TRUE)
plot(ti,ELL.norm.do,main="Smoothed exposures",xlab='age',ylab='',type='l',
lwd=2,ylim=yy)
yy2<-range(ELL.norm.do[ind.ages],na.rm=TRUE)
plot(ti2,ELL.norm.do[ind.ages],main="Smoothed exposures",xlab='age',ylab='',
type='l',lwd=2,ylim=yy2)

# Revert the changes made in the graphics options
par(old.par)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
