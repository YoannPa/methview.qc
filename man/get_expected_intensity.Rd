% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_fun.R
\name{get_expected_intensity}
\alias{get_expected_intensity}
\title{Provides expected intensity for a given methylation array probe ID.}
\usage{
get_expected_intensity(
  DT.QC.meta,
  probe.id,
  channel.names = c("Cy3 - Green", "Cy5 - Red")
)
}
\arguments{
\item{DT.QC.meta}{A \code{data.table} with methylation array quality
control metadata obtained using
\link{load_metharray_QC_meta}.}

\item{probe.id}{A \code{character} string matching a methylation array
probe ID listed in 'DT.QC.meta'.}

\item{channel.names}{A \code{character} vector specifying the names of the 2
color channels to use for creating the data.table of
expected intensity (Default:
channel.names = c("Cy3 - Green", "Cy5 - Red")).}
}
\value{
A \code{data.table} with a 'Channel' column and an
        'Expected intensity' column providing the probe intensity expected
        for green and red channels.
}
\description{
Provides expected intensity for a given methylation array probe ID.
}
\examples{
#Create the data.table with quality control metadata from HM450K
dt.meta <- load_metharray_QC_meta(array.meta = "controls450")
#Check expected intensity for QC probe "27630314"
get_expected_intensity(DT.QC.meta = dt.meta, probe.id = "27630314")
}
\author{
Yoann Pageaud.
}
