% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_functions.R
\name{makeH5Index}
\alias{makeH5Index}
\title{Create a index table for an HDF5 dataset based on column}
\usage{
makeH5Index(filename, datasetName, indexCol, dataIndexName, rowSize = 1e+06)
}
\arguments{
\item{filename}{String specifying the h5 file.}

\item{datasetName}{The dataset to be indexed.}

\item{indexCol}{The column that stores indexing entries.}

\item{dataIndexName}{The index dataset name.}

\item{rowSize}{A parameter to specify how many rows of data to load at any time.}
}
\value{
The data frame that specifies the row indices.
}
\description{
\code{makeH5Index} takes a specified HDF5 dataset in a given file, and
produces an index of the data based off of entries in a specified column of the dataset.
This function assumes \code{datasetName} is sorted by the column that will be used as an index.
If not, the user should use \code{sortByIndex} to produce a sorted dataset for indexing.
}
\examples{
makeH5Index("foo.h5","/theData",3,"frameIndices")
}
