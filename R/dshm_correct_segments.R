#' Correcting segments for overlapping features
#'
#' \code{dshm_correct_segments} uses the object \code{int.ids} created with the function \code{dshm_check_segments} to eliminate overlapping fetures between segments.
#'
#' @param data Segments as SpatialPolygonsDataFrame. You can use the function \code{dshm_split_transects} to create segments from transect lines.
#' @param intersections The \code{int.ids} object created with the function \code{dshm_check_segments}. This is a list for all segment intersections.
#' @return Segments as SpatialPolygonsDataFrame corrected for overlapping features.
#' @details For more information about splitting transects into segments as well as checking and correcting segments you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{split_transects.pdf} tutorial.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}

#' @export

dshm_correct_segments <- function(data, intersections){

  mod.seg <- list() #empty list for modified segments
  mod.seg[[length(data)]] <- data[length(data),] #placing in the last slot of 'mod.seg' the last polygon in 'data'. This automatically assign 'NULL' to all slots in 'mod.seg'
  for (i in 1:length(intersections)){
    if(!is.null(intersections[[i]])){

      for(j in 1:length(intersections[[i]])){ #for each polygon within each slot

        s <- rgeos::gDifference(data[intersections[[i]][j],], data[i,]) #difference between [[i]][j] polygons and [[i]] polygon. Overlapping features of [[i]][j] polygons with [[i]] polygon are took out

        if(i>1&&!is.null(mod.seg[[intersections[[i]][j]]])){ #if a polygon has been already modified, it takes it from 'mod.seg' rather than from 'data'
          s <- rgeos::gDifference(mod.seg[[intersections[[i]][j]]], data[i,])
        }

        s$Transect.Label <- data[intersections[[i]][j],]$Transect.Label #transect label
        s$Sample.Label <- data[intersections[[i]][j],]$Sample.Label #sample label
        s$length <- data[intersections[[i]][j],]$length #length
        s$area <- raster::area(s)/10^6 #recalculate area
        mod.seg[[intersections[[i]][j]]] <- s #storing modified segment
      }
    }
  }

  corr.seg <- list() #empty list for corrected segments

  for (i in 1:length(data)){ #if the ith slot in 'mod.seg' is not NULL, it pastes the corrected segments, if NULL it pastes the original segment from 'data'
    if(!is.null(mod.seg[[i]])){

      corr.seg[[i]] <- mod.seg[[i]]

    } else {

      corr.seg[[i]] <- data[i,]

    }
  }

  corr.seg.bind <- do.call(raster::bind, corr.seg) #binding all corrected segments together
  return(corr.seg.bind) #returning the bound corrected segments object as SpatialPolygonsDataFrame

}
