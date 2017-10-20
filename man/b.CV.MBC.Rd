\name{b.CV.MBC}
\alias{b.CV.MBC}
\title{
Least Squares Cross-Validation for Multiplicative Bias Corrected Hazard Estimators}
\description{
Bandwidth selection for multiplicatively bias corrected local linear hazard estimation using least squares cross-validation}
\usage{
b.CV.MBC(grid.b, nb , K = "sextic", xi, Oi, Ei, wei = "same")
}
\arguments{
  \item{grid.b}{
Optional. A vector of bandwidths to minimise the cross-validation score. If not specified it will be considered an equally-spaced grid of \code{nb} bandwidths between "amp/(M+1)" and "amp/2" for "amp" being the range of \code{xi} and "M" its length.}
  \item{nb}{
Optional. The number of bandwidths used to minimise the cross-validation score. If "grid.b" is provided then the argument "nb" will be ignored (if specified). }
  \item{K}{Indicates the kernel function to be considered in the local linear hazard estimator. Choose between values "epa" (for the epanechnikov kernel) or "sextic" (see details of \code{\link{hazard.MBC}} for the definition).
}
  \item{xi}{Vector of time points where the count data are given.
}
  \item{Oi}{Vector with the number (counts) of occurrences observed at each time point (xi). 
}
  \item{Ei}{Vector with the observed exposure at each time point (xi).  
}
  \item{wei}{Indicates the weights used in the cross-validation score. Choose between the value "exposure" or "same". See details below.}
}
\details{It is assumed that the data are given as count data i.e. number of occurrences and exposures.

If the cross-validation score is strictly increasing or decreasing then a warning will be shown together with the cross-validated bandwidth (in this case one of the extremes in "grid.b").

The cross-validation score is defined with two different weighting functions. This is controlled with the parameter \code{wei}. By default \code{wei="exposure"} that means that only areas where the exposure is significant contribute to the criterion. Specify \code{wei="same"} to allow  all time points contribute the same to the criterion (see Gamiz et al. 2017).}
\value{
%  cv.res<-list(bcv=bcv,ind.cv=ind.cv,cv.values=cv.values,grid.b=grid.b)
  \item{bcv }{The cross-validated bandwidth.}
  \item{ind.cv }{The position of the cross-validated bandwidth into "grid.b".}
  \item{cv.values }{The values of the cross-validation score for each bandwidth in "grid.b".}
  \item{b.grid }{The grid of bandwidths where the cross-validation score has been evaluated.}
}
\references{
Gamiz, M.L., Martinez-Miranda, M.D. and Nielsen, J.P. (2017). Multiplicative local linear hazard estimation and best one-sided cross-validation. Available 	 at http://arxiv.org/abs/1710.05575  

Nielsen, J.P. and Tanggaard, C. (2001). Boundary and bias correction in kernel hazard estimation. Scandinavian Journal of Statistics, 28, 675-698.
}
\author{
Gamiz, M.L., Martinez-Miranda, M.D. and Nielsen, J.P.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hazard.MBC}},\code{\link{b.BO.MBC}},\code{\link{b.CV}}
}
\examples{
data(Iceland)
Oi<-Iceland$D
Ei<-Iceland$E
ti<-40:110  # time is age and it goes from 40 to 110 years
my.bs<-seq(50,80,length.out=30)
res.cv<-b.CV.MBC(grid.b=my.bs,K="sextic",xi=ti,Oi=Oi,Ei=Ei,wei="same")
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
