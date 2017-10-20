\name{discretise.data}
\alias{discretise.data}
\title{Aggregate data in the form of occurrences and exposures
}
\description{
Aggregate data in the form of occurrences and exposures from individual survival data (possibly right censored and/or left truncated).}
\usage{
  discretise.data(Li, Zi, deltai, xi, M)
}
\arguments{
  \item{Li}{Vector of truncation levels: the datum is registered only if the life time is greater than the truncation level. 
}
  \item{Zi}{Vector of observed life times (Zi=min(Ti,Ci) with Ci censoring value and Ti the true life time). 
}
  \item{deltai}{Vector with non-censoring indicator values (0 if datum is censored, 1 otherwise).  
}
  \item{xi}{Optional. Vector with the grid of time points where the occurrences and exposures should  be calculated. If not provide the grid is calculated automatically.
}  
  \item{M}{Optional. A positive scalar used as the grid size. If not provided it is chosen automatically.
}
}

\details{  The hazard estimators and bandwidth selectors available in the DOvalidation package work from data aggregated in the form of occurrences and exposures. This function can be used to work with individual survival data in the form (Li,Zi,deltai) -- left-truncation level (Li), observed time (Zi) and non-censoring indicator (deltai). If data are not truncated then Li can be chosen as 0.      
}

\value{
%hazard.res<-list(x=x,hLL=hLL)
  \item{xi }{Vector with the time grid points. }
  \item{Oi }{Vector with the calculated occurrences at the time grid points.}
   \item{Ei }{Vector with the calculated exposures at the time grid points.}
}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2016). Double one-sided cross-validation of local linear hazards. Journal of the Royal Statistical Society B, 78, 755-779.
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
## A simple example with right censored data
Zi<-c(3,6,7,7,8,10,11,11,11,12,13,13,14,16,20,20,22,32,34,36)
n<-length(Zi)
Li<-deltai<-rep(0,n) 
ind.cens<-c(1,3,4,8,9,13,14,15,16)
deltai[-ind.cens]<-1
## Obtain the occurrences and exposures for a grid of 6 time points
res<-discretise.data(Li,Zi,deltai,M=6)
## Now calculate the local linear hazard estimator
hazard.LL(res$xi,res$Oi,res$Ei,res$xi,b=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
