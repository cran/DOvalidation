\name{b.CV}
\alias{b.CV}
\title{
Least Squares Cross-Validation for Local Linear Hazards}
\description{
Bandwidth selection for local linear hazard estimation using least squares cross-validation}
\usage{
b.CV(grid.b, nb = 50, K = "epa", xi, Oi, Ei)
}
\arguments{
  \item{grid.b}{
Optional. A vector of bandwidths to minimise the cross-validation score. If not specified it will be considered an equally-spaced grid of "nb" bandwidths between "amp/(M+1)" and "amp/2" for "amp" being the range of "xi" and "M" the lenght of xi.}
  \item{nb}{
The number of bandwidths used to minimise the cross-validation score. If "grid.b" is provided then the argument "nb" will be ignored (if specified). }
  \item{K}{Indicates the kernel function to be considered in the local linear hazard estimator. Choose between values "epa" (for the epanechnikov kernel) or "sextic" (see details of \code{\link{hazard.LL}} for the definition).
}
  \item{xi}{Vector of time points where the count data are given.
}
  \item{Oi}{Vector with the number (counts) of occurrences observed at each time point (xi). 
}
  \item{Ei}{Vector with the observed exposure at each time point (xi).  
}
}
\details{It is assumed that the data are given as count data i.e. number of occurrences and exposures.

If the cross-validation score is strictly increasing or decreasing then a warning will be shown together with the cross-validated bandwidth (in this case one of the extremes in "grid.b".)}
\value{
%  cv.res<-list(bcv=bcv,ind.cv=ind.cv,cv.values=cv.values,grid.b=grid.b)
  \item{bcv }{The cross-validated bandwidth.}
  \item{ind.cv }{The position of the cross-validated bandwidth into "grid.b".}
  \item{cv.values }{The values of the cross-validation score for each bandwidth in "grid.b".}
  \item{b.grid }{The grid of bandwidths where the cross-validation score has been evaluated.}
}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2014). Do-Validating Local Linear Hazards. 
Available at SSRN: http://dx.doi.org/10.2139/ssrn.2504497

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
\code{\link{hazard.LL}},\code{\link{b.OSCV}}
}
\examples{
data(UK)
Oi<-UK$D
Ei<-UK$E 
ti<-40:110  # time is age and it goes from 40 to 110 years
M<-length(ti)
my.bs<-seq(1.5,40,length=100)
res.cv<-b.CV(grid.b=my.bs,K="sextic",xi=ti,Oi=Oi,Ei=Ei)
bcv<-res.cv$bcv
cv.values<-res.cv$cv.values
bs<-my.bs
cvs<-cv.values
plot(bs,cvs,main="Cross-validation score",xlab="Bandwidth")
print(paste("The cross-validated bandwidth is:", bcv,sep=" "))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
