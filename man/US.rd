\name{US}
\alias{US}
\docType{data}
\title{
US Female Mortality Data
}
\description{
Mortality data of women in the calendar year 2006 from United States. 
The data were obtained from the Humam Mortality Database. Only ages from 40 to 110 have been included.
}
\usage{data(US)}
\format{
 This data frame contains 71 rows and the following 2 columns.
   \describe{
    \item{\code{D}}{Death counts for women of ages between 40 and 110 during the calendar year 2006. Some of these numbers are estimates (of population size or numbers of deaths), not actual counts, and therefore may be expressed as non-integers.
    
    }
    \item{\code{E}}{"Person-years" lived in the female population during the year 2006 for each age-group (from 40 to 110)}  }
}
\source{
Human Mortality Database. University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). Available at www.mortality.org or www.humanmortality.de
}
\references{
Gamiz, M.L., Mammen, E., Martinez-Miranda, M.D. and Nielsen, J.P.(2014). Do-Validating Local Linear Hazards. 
Available at SSRN: http://dx.doi.org/10.2139/ssrn.2504497

Spreeuw, J., Nielsen, J.P. and Jarner, S.F. (2013). A visual test of mixed hazard models, SORT, 37, 149-170.
}
\examples{
data(US)
}
\keyword{datasets}
