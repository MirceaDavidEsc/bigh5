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




#' Read all data belonging to a particular entry.
#'
#' @param hdf5Identifier The name of the hdf5 file.
#' @param datasetName The dataset that contains the data you wish to load.
#' @param indicesName A mapping of factor indices to rows in the data frame to load.
#' @param frameNumbers The factors for which you fish to load your data.
#'
#' @return
#' @export
#'
#' @examples
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


#' Read specified chunks of data in an HDF5 dataset.
#'
#' This function reads chunks from an HDF5 file with the name provided in filename and the chunks delimited by a data
#'  frame called chunkDefineTable. This function should in general not be called, except through the readDataAtFrame function.
#'
#' chunkDefineTable must contain two columns with the names dataChunkStart and dataChunkEnd that specify the blocks
#' of data to load.
#'
#' @param myFileName The HDF5 file from which to load the data
#' @param datasetName The dataset name within the HDF5 file
#' @param chunkDefineTable A table that contains a dataChunkStart and dataChunkEnd column.
#'
#' @return
#' @export
#'
#' @examples
readLongformAtChunks = function(myFileName,datasetName,chunkDefineTable) {
   indicesList = mapply(FUN = function(x,y) x:y, chunkDefineTable$dataChunkStart, chunkDefineTable$dataChunkEnd)
  indicesVector = unlist(indicesList)
  framesLongForm = tryCatch(
    {as.data.frame(h5read(myFileName,datasetName, index = list(indicesVector,NULL))) #Load all columns
      }, error = function(cond) {return(NA)}
  )
  H5close()
  return(framesLongForm)
}



#' Transform matrix format to flatfile format.
#'
#'  Transforms a matrix of values into a long-form data frame.
#'
#' @param myarr A N-dimensional matrix.
#'
#' @return A data frame with (X*Y*Z...) columns rows and N columns.
#' @export
#'
#' @examples
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
#'This function takes in a 4-column data frame that contains the requisite data on a velocity field, where each row is
#'a vector in space which has an x and y position and an x and y velocity component. It returns a ggplot object.
#'
#' @param frameData A 4-column data frame consisting of (from left to right) the x position, y position, x-component, and y-component of the velocity.
#' @param arrowl A scalar multiple to apply to the velocity components, for visualization.
#' @param colormapped If F, plot as arrow line segments. Else, plot as color heatmap.
#'
#' @return A ggplot2 object that is the plot of the velocity field.
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
      scale_color_gradientn(name = expression(theta), colours = cbPalette, limits = c(-pi,pi), breaks=c(-pi, 0, pi),
                            labels=c(expression(paste("-",pi,sep="")), 0, expression(paste(pi)))) +
      scale_fill_gradientn(name = expression(theta), colours = cbPalette, limits = c(-pi,pi), breaks=c(-pi, 0, pi),
                           labels=c(expression(paste("-",pi,sep="")), 0, expression(paste(pi))))

  } else {
    p = ggplot(frameData, aes(x, y, xend=x+u, yend=y+v)) +
      geom_segment(arrow=arrow(angle=20,length=unit(0.2,"cm")))
  }
  p = p + coord_fixed(ratio = 1) + theme(axis.title = element_blank())
  return(p)
}
