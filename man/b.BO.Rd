\name{b.BO}
\alias{b.BO}
\title{
Best One-Sided Cross-Validation for Local Linear Hazards}
\description{
Bandwidth selection for local linear hazard estimation using best one-sided cross-validation}
\usage{
b.BO(grid.b, nb , K = "sextic", type.bo = "Oi", xi, Oi, Ei, wei = "same")
}
\arguments{
  \item{grid.b}{Optional. A vector of bandwidths to minimise the cross-validation score. If not specified it will be considered an equally-spaced grid of \code{nb} bandwidths between "amp/(M+1)" and "amp/2" for "amp" being the range of \code{xi} and "M" its length.}
  \item{nb}{
Optional. The number of bandwidths used to minimise the cross-validation score. If \code{grid.b} is provided then the argument \code{nb} will be ignored (if specified). }
  \item{K}{Indicates the kernel function to be considered in the local linear hazard estimator. Choose between values \code{"epa"} (for the epanechnikov kernel) or \code{"sextic"} (see details in\code{\link{hazard.LL}} for the definition).
}
  \item{type.bo}{Choose between \code{"Oi"} or \code{"Ei"} to find the best side using the occurrences or the exposures, respectively.
}
  \item{xi}{Vector of time points where the count data are given.
}
  \item{Oi}{Vector with the number (counts) of occurrences observed at each time point (\code{xi}). 
}
  \item{Ei}{Vector with the observed exposure at each time point (\code{xi}).  
}
  \item{wei}{Indicates the weights used in the cross-validation score. Choose between the value \code{"exposure"} or \code{"same"}. See details below.}
}
\details{
It is assumed that the data are given as count data i.e. number of occurrences and exposures.

The BO-validated bandwidth is calculated as the minimizer of a cross-validation score with a indirect kernel. If the score is strictly increasing or decreasing then a warning will be shown together with the selected bandwidth (in this case one of the extremes in \code{grid.b}, adjusted by the rescaling constant, which is 0.5371  for the Epanechnikov kernel, and 0.5874 for the sextic kernel).

The score is defined with two different weighting functions. This is controlled with the parameter \code{wei}. By default \code{wei="exposure"} that means that only areas where the exposure is significant contribute to the criterion. Specify \code{wei="same"} to allow  all time points contribute the same to the criterion (see Gamiz et al. 2017).}
\value{
%  list(bbo=bbo,ind.bo=ind.bo,cvbo.values=cvbo.values,Cval=Cval, grid.b=grid.b)
  \item{bbo }{The onesided cross-validated bandwidth.}
  \item{ind.bo }{The position of the best one-sided cross-validated bandwidth into "grid.b".}
  \item{cvbo.values }{The values of the cross-validation score for each bandwidth in \code{grid.b}.}
  \item{b.grid }{The grid of bandwidths where the score has been evaluated.}
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
\code{\link{hazard.LL}},\code{\link{b.OSCV}},\code{\link{b.CV}}
}
\examples{
data(Iceland)
Oi<-Iceland$D
Ei<-Iceland$E
ti<-40:110  # time is age and it goes from 40 to 110 years
my.bs<-seq(20,40,length.out=30)
res.bo<-b.BO(grid.b=my.bs,K="sextic",type.bo="Oi",xi=ti,Oi=Oi,Ei=Ei,wei="same")
bbo<-res.bo$bbo
cvs<-res.bo$cvbo.values
plot(my.bs,cvs,main="BO-validation score",xlab="Bandwidth")
print(paste("The best one-sided cross-validated bandwidth is:", bbo,sep=" "))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{nonparametric}
\keyword{survival}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
