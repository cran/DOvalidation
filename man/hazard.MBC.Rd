\name{hazard.MBC}
\alias{hazard.MBC}
\title{
Multiplicative Bias Corrected Hazard Estimator}
\description{
Multiplicatively bias corrected local linear estimator of the unidimensional hazard with natural weighting introduced by Nielsen and Tanggaard (2001).}
\usage{
  hazard.MBC(xi, Oi, Ei, x, b, K="sextic", Ktype="symmetric")
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
}

\details{The estimator is calculated assuming that the data are given as count data i.e. number of occurrences and exposures.
The function allows to consider two different kernels using the argument \code{K}. These are: Epanechnikov, K(u)=.75*(1-u^2)*(abs(u)<1), and sextic K(u)=(3003/2048)*(1-(u)^2)^6)*(abs(u)<1). The argument \code{Ktype} will define the usual estimator with whole support kernel as it is defined by \code{K} or the one-sided versions using left-sided kernel, 2*K(u)*(u<0), or right-sided kernel 2*K(u)*(u>0). See more details in Gamiz et al. (2017).  
}

\value{
  \item{x }{Vector (or scalar) with the (time) grid points where the hazard estimator has been evaluated. }
  \item{hMBC }{Vector (or scalar) with the resulting hazard estimates at grid points \code{x}.}
}
\references{
Gamiz, M.L., Martinez-Miranda, M.D. and Nielsen, J.P. (2017). Multiplicative local linear hazard estimation and best one-sided cross-validation. Available 	 at http://arxiv.org/abs/1710.05575  

Nielsen, J.P. and Tanggaard, C. (2001). Boundary and bias correction in kernel hazard estimation. Scandinavian Journal of Statistics,28, 675-698.
}
\author{
Gamiz, M.L., Martinez-Miranda, M.D. and Nielsen, J.P.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hazard.LL}}
}
\examples{
data(Iceland)
Oi<-Iceland$D
Ei<-Iceland$E
ti<-40:110  # time is age and it goes from 40 to 110 years
res<-hazard.MBC(xi=ti,Oi=Oi,Ei=Ei,x=ti,b=48.7)
plot(ti,res$hMBC,main='Hazard estimate',xlab='age',ylab='',type='l',lwd=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
