% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_functions.R
\name{readDataAtFrame}
\alias{readDataAtFrame}
\title{Read all data belonging to a particular entry.}
\usage{
readDataAtFrame(hdf5Identifier, datasetName, indicesName, frameNumbers)
}
\arguments{
\item{hdf5Identifier}{The name of the hdf5 file.}

\item{datasetName}{The dataset that contains the data you wish to load.}

\item{indicesName}{A mapping of factor indices to rows in the data frame to load.}

\item{frameNumbers}{The factors for which you fish to load your data.}
}
\description{
Read all data belonging to a particular entry.
}
