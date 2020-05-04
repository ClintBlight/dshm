#' Checking for overlapping features between segments
#'
#' \code{dshm_check_segments} checks if there are overlapping features between segments created by \code{dshm_split_transects}.
#'
#' @param data Segments as SpatialPolygonsDataFrame. You can use the function \code{dshm_split_transects} to create segments from transect lines.
#' @return Nothing if there are no overlapping features between segments. In case of intersecting polygons the function returns:
#' \itemize{
#'   \item int.ids: a map for all intersections. This is used by the function \code{dshm_correct_segments}.
#'   \item inter.Polys: a list of all intersecting polygons.
#' }
#' @details For more information about splitting transects into segments as well as checking and correcting segments you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{split_transects.pdf} tutorial.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}
#' @export
#'
dshm_check_segments <- function(data){

  mat <- rgeos::gIntersects(rgeos::gBuffer(data, byid = TRUE, width = 0), byid = TRUE, checkValidity = TRUE) #matrix for intersections, 'TRUE' if interaction
  int.polys <- list() #empty list to store intersections as SpatialPolygons
  int.ids <- list() #empty list to store intersecting polygons ids

  for (i in 1:length(mat[1,])){

    colu <- mat[,i][-c(1:i)] #taking out the upper part of the matrix while iterating

    if(sum(colu)==0){ #if all values are FALSE (i.e. no interaction) then go to next iteration

      next

    }

    int.poly.ids <- colu[colu==TRUE] #select the polygons sharing overlapping areas
    int.poly <- rgeos::gIntersection(data[i,],
                                     data[as.numeric(rownames(as.data.frame(int.poly.ids))),],
                                     drop_lower_td = TRUE) #extract just the intersection polygons

    if(!class(int.poly)[1]=="SpatialPolygons"){ #if the object is not a polygon then go the next iteration
      next
    }

    if(round(raster::area(int.poly))==0){ #if the object area is below 1 m^2 then go the next iteration
      next
    }

    int.polys[[i]] <- int.poly #store intersction polygons
    int.ids[[i]] <- as.numeric(rownames(as.data.frame(int.poly.ids))) #store ids of intersecting polygons
  }

  if(length(int.polys)==0){ #if 'int.polys' is empty then nothing to do

    cat(paste("No overlapping features. Segments are OK."))

  } else { #if there are overlapping features

    int.polys.nonulls <- int.polys[!sapply(int.polys, is.null)] #take out NULL slots from 'int.polys' list
    cat(paste(length(int.polys.nonulls)," overlapping features found. Go to 'dshm_correct_segments' function."))
    
    if (length(int.polys.nonulls) > 1 ) {
      inter.Polys <- do.call(raster::bind, int.polys.nonulls)
    } else {
      inter.Polys <- int.polys.nonulls[[1]]
    }
    
    return(list(inter.IDs = int.ids, inter.Polys = inter.Polys))

  }
}
