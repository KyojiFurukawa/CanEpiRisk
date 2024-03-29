% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCER.r
\name{plotCER}
\alias{plotCER}
\title{Plot the age-specific risk (attributable probability rate) due to radiation exposure}
\usage{
plotCER(
  cer,
  xlab = "Attained age",
  ylab = "Attributable probability rate",
  main = "Attributable probability rate",
  col = 1
)
}
\arguments{
\item{cer}{a CER class object (returned from function CER)}
}
\value{
risk information(vector)
}
\description{
Plot the age-specific risk (attributable probability rate) due to radiation exposure
}
\examples{
   # The following examples use default data provided in the JHPSrisk package
   # for riskmodel (LSS R14 all solid cancer model),
   #     baseline (all solid cancer mortality rates in Japan 2018) and
   #     mortality  (all cause mortality rates in Japan 2018)

   # Cumulated excess risk for male exposed to 0.1 Gy at age 10
   #  followed up to age 90 with err transfer
   cer <- CER( agex=10, doseGy=0.1, sex=1, maxage=90 )
   plot(cer)

#'@export
}
