# Need to update package DESCRIPTION to make these dependencies.
#suppressWarnings(suppressPackageStartupMessages(require(rhdf5)))
#suppressWarnings(suppressPackageStartupMessages(require(dplyr)))


#' Create a index table for an HDF5 dataset based on column
#'
#' \code{makeH5Index} takes a specified HDF5 dataset in a given file, and
#' produces an index of the data based off of entries in a specified column of the dataset.
#' This function assumes \code{datasetName} is sorted by the column that will be used as an index.
#' If not, the user should use \code{sortByIndex} to produce a sorted dataset for indexing.
#'
#' @param filename String specifying the h5 file.
#' @param datasetName The dataset to be indexed.
#' @param dataIndexName The index dataset name.
#' @param indexCol The column that stores indexing entries.
#' @param rowSize A parameter to specify how many rows of data to load at any time.
#'
#' @return The data frame that specifies the row indices.
#' @export
#'
#' @examples
#' makeH5Index("foo.h5","/theData",3,"frameIndices")
makeH5Index = function(filename, datasetName, indexCol, dataIndexName, rowSize = 1000000) {

  h5comp = h5ls(filename)
  maxRows = as.numeric(strsplit(h5comp$dim, " x ")[[1]][1])

  startSeq = seq(from = 1, to = maxRows, by = rowSize)

  countTables = list()
  for (thisStart in startSeq) {
    thisEnd = min(thisStart + rowSize - 1, maxRows)
    h5data = h5read(filename, datasetName, index = list(thisStart:thisEnd, 3))
    thisPoints = as.data.frame(h5data)
    pointCount = as.data.frame(table(thisPoints$V1))
    pointCount$Var1 = as.numeric(as.character(pointCount$Var1))

    if (length(countTables) > 0) {
      # Merge the first entries i
      lastEntry = tail(countTables, 1)[[1]]
      if (tail(lastEntry$Var1, 1) == pointCount$Var1[1]) {
        lastEntry$Freq[length(lastEntry$Freq)] = lastEntry$Freq[length(lastEntry$Freq)] + pointCount$Freq[1]
        pointCount = pointCount[-1,]
        countTables[[length(countTables)]] = lastEntry
      }
    }
    countTables[[length(countTables) + 1]] = pointCount
  }
  countDF = bind_rows(countTables)
  countDF$dataChunkEnd = cumsum(countDF$Freq)
  numFrames = dim(countDF)[1]
  countDF$dataChunkStart = 1 + c(0,countDF$dataChunkEnd[1:(numFrames-1)])
  countDF = countDF[, c(1, 4, 3)]
  colnames(countDF) = c("Frame", "dataChunkStart", "dataChunkEnd")
  H5close()
  h5write(countDF, filename, dataIndexName)
  return(countDF)
}


readDataAtFrame = function(hdf5Identifier,datasetName,indicesName,frameNumbers) {
  # Given the name of the file and the numeric value of the Frame for which you wish to load velocity field data, returns that longform data
  # Args:
  #   hdf5Identifier: the file name of the hdf5 file containing long-format vector field data.
  #   frameNumbers: A numeric vector of integers specifying the frames for which data is to be loaded.

  longformIndices = as.data.frame(h5read(file = hdf5Identifier, name = indicesName))
  colnames(longformIndices) = c("Frame","dataChunkStart","dataChunkEnd")
  theseFrames = longformIndices[longformIndices$Frame %in% frameNumbers,]
  theseFrames = theseFrames[,c("dataChunkStart", "dataChunkEnd")]
  #theseFrames = filter(longformIndices, Frame %in% frameNumbers) %>% select(dataChunkStart,dataChunkEnd)

  dataOfFrames = readLongformAtChunks(hdf5Identifier,datasetName,theseFrames)
  gc()
  return(dataOfFrames)
}


readLongformAtChunks = function(myFileName,datasetName,chunkDefineTable) {
  # Reads chunks from an HDF5 file with the name provided in filename
  # and the chunks delimited by a data frame called chunkDefineTable.
  # Args:
  #   chunkDefineTable: A table consisting of two columns specifying the start and end indices of the data chunk(s) to be read.
  #   filename: A string specifying the name of the HDF5 file.
  # Output:
  #   The data stored in the HDF5 file at the rows specified by fileChunks. This is concatenated as one big data frame.
  indicesList = mapply(FUN = function(x,y) x:y, chunkDefineTable$dataChunkStart, chunkDefineTable$dataChunkEnd)
  indicesVector = unlist(indicesList)
  framesLongForm = tryCatch(
    {as.data.frame(h5read(myFileName,datasetName, index = list(indicesVector,NULL))) #Load all columns
      }, error = function(cond) {return(NA)}
  )
  H5close()
  return(framesLongForm)
}


reshapeArrayLongForm = function(myarr) {
  ### FUNCTION: Reshapes an array into a long-form data frame, while saving the indices of the data into the data frame.
  newarr = myarr
  dim(newarr) = c(prod(dim(myarr)),1,1)
  zcoord = rep(1:dim(myarr)[3],each=dim(myarr)[1]*dim(myarr)[2])
  ycoord = rep(rep(1:dim(myarr)[2],each=dim(myarr)[1]),dim(myarr)[3])
  xcoord = rep(c(1:dim(myarr)[1]),dim(myarr)[2]*dim(myarr)[3])
  longform = data.frame(xc=xcoord,yc=ycoord,zc=zcoord,val=newarr)
}



#' Plot of velocity field
#'
#' @param frameData
#' @param arrowl
#' @param colormapped
#'
#' @return
#' @export
#'
#' @examples
quiverPlot <- function(frameData,arrowl, colormapped=F) {
  require(ggplot2)
  colnames(frameData) = c("x","y","u","v")
  frameData = mutate(frameData,v=v*arrowl, u=u*arrowl)

  if (colormapped) {
    frameData = frameData %>% mutate(speed = sqrt(v^2 + u^2), angle = atan2(x = u, y = v))
    cbPalette = c("magenta","red","yellow","green","cyan","blue","magenta")
    p = ggplot(frameData, aes(x, y, fill=angle, alpha=speed)) + geom_tile() +
      scale_color_gradientn(colours = cbPalette, limits = c(-pi,pi), breaks=c(-pi, 0, pi),
                            labels=c(expression(paste("-",pi,sep="")), 0, expression(paste(pi)))) +
      scale_fill_gradientn(colours = cbPalette, limits = c(-pi,pi), breaks=c(-pi, 0, pi),
                           labels=c(expression(paste("-",pi,sep="")), 0, expression(paste(pi))))

  } else {
    p = ggplot(frameData, aes(y, x, yend=x+v, xend=y+u)) +
      geom_segment(arrow=arrow(angle=20,length=unit(0.2,"cm")))
  }
  p = p + coord_fixed(ratio = 1) + theme(axis.title = element_blank())
  return(p)
}
