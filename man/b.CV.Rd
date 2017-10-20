\name{b.CV}
\alias{b.CV}
\title{
Least Squares Cross-Validation for Local Linear Hazards}
\description{
Bandwidth selection for local linear hazard estimation using least squares cross-validation}
\usage{
b.CV(grid.b, nb , K = "epa", xi, Oi, Ei, wei = "exposure")
}
\arguments{
  \item{grid.b}{
Optional. A vector of bandwidths to minimise the cross-validation score. If not specified it will be considered an equally-spaced grid of \code{nb} bandwidths between "amp/(M+1)" and "amp/2" for "amp" being the range of \code{xi} and "M" its length.}
\item{nb}{
Optional. The number of bandwidths used to minimise the cross-validation score. If \code{grid.b} is provided then the argument \code{nb} will be ignored (if specified). }
  \item{K}{Indicates the kernel function to be considered in the local linear hazard estimator. Choose between values \code{"epa"} (for the epanechnikov kernel) or \code{"sextic"}  (see details of \code{\link{hazard.LL}} for the definition).
}
  \item{xi}{Vector of time points where the count data are given.
}
  \item{Oi}{Vector with the number (counts) of occurrences observed at each time point (\code{xi}). 
}
  \item{Ei}{Vector with the observed exposure at each time point (\code{xi}).  
}
  \item{wei}{Indicates the weights used in the cross-validation score. Choose between the value \code{"exposure"} or \code{"same"}. See details below.}
}
\details{It is assumed that the data are given as count data i.e. number of occurrences and exposures.

If the cross-validation score is strictly increasing or decreasing then a warning will be shown together with the cross-validated bandwidth (in this case one of the extremes in \code{grid.b}).

The cross-validation score is defined with two different weighting functions. This is controlled with the parameter \code{wei}. By default \code{wei="exposure"} that means that only areas where the exposure is significant contribute to the criterion. Specify \code{wei="same"} to allow  all time points contribute the same to the criterion (see Gamiz et al. 2017).}
\value{
%  cv.res<-list(bcv=bcv,ind.cv=ind.cv,cv.values=cv.values,grid.b=grid.b)
  \item{bcv }{The cross-validated bandwidth.}
  \item{ind.cv }{The position of the cross-validated bandwidth into \code{grid.b}.}
  \item{cv.values }{The values of the cross-validation score for each bandwidth in \code{grid.b}.}
  \item{b.grid }{The grid of bandwidths where the cross-validation score has been evaluated.}
}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2016). Double one-sided cross-validation of local linear hazards. Journal of the Royal Statistical Society B, 78, 755-779.

Gamiz, M.L., Martinez-Miranda, M.D. and Nielsen, J.P. (2017). Multiplicative local linear hazard estimation and best one-sided cross-validation. Available 	 at http://arxiv.org/abs/1710.05575 

Nielsen, J.P. and Tanggaard, C. (2001). Boundary and bias correction in kernel hazard estimation. Scandinavian Journal of Statistics, 28, 675-698.
}
\author{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hazard.LL}},\code{\link{b.OSCV}},\code{\link{b.BO}}
}
\examples{
data(UK)
Oi<-UK$D
Ei<-UK$E 
ti<-40:110  # time is age and it goes from 40 to 110 years
M<-length(ti)
my.bs<-seq(1,5,length=50)
res.cv<-b.CV(grid.b=my.bs,K="sextic",xi=ti,Oi=Oi,Ei=Ei)
bcv<-res.cv$bcv
cv.values<-res.cv$cv.values
plot(my.bs,cv.values,main="Cross-validation score",xlab="Bandwidth")
print(paste("The cross-validated bandwidth is:", bcv,sep=" "))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
