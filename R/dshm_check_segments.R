#' Check if there are overlapping features between SpatialPolygons
#'
#' @param data Segments or other datasets containing polygons as SpatialPolygonsDataFrame.
#' @return Nothing. It tells if there are intersecting polygons and how many are they. In case of match, it stores in the Global Environent a list 'int.ids' containing the ids of all intersecting polygons. This is used by 'dshm_correct_segments'.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}

#' @export
dshm_check_segments<-function(data){

  mat<-rgeos::gIntersects(rgeos::gBuffer(data, byid=TRUE, width=0),byid=TRUE,checkValidity=TRUE) #matrix for intersections, 'TRUE' if interaction
  int.polys<-list() #empty list to store intersections as SpatialPolygons
  int.ids<-list() #empty list to store intersecting polygons ids
  for (i in 1:length(mat[1,])){
    colu<-mat[,i][-c(1:i)] #taking out the upper part of the matrix while iterating
    if(sum(colu)==0){ #if all values are FALSE (i.e. no interaction) then go to next iteration
      next
    }
    int.poly.ids<-colu[colu==TRUE] #select the polygons sharing overlapping areas
    int.poly<-rgeos::gIntersection(data[i,],data[as.numeric(rownames(as.data.frame(int.poly.ids))),],drop_lower_td = TRUE) #extract just the intersection polygons
    if(!class(int.poly)[1]=="SpatialPolygons"||round(raster::area(int.poly)*10^4)==0){ #if the object is not a polygon or if the area of the polygon is below cm^2 scale then go the next iteration
      next
    }
    int.polys[[i]]<-int.poly #store intersction polygons
    int.ids[[i]]<-as.numeric(rownames(as.data.frame(int.poly.ids))) #store ids of intersecting polygons
  }

  if(length(int.polys)==0){ #if 'int.polys' is empty then nothing to do
    cat(paste("No overlapping features. Segments are OK."))
  } else { #if there are overlapping features
    int.polys.nonulls<-int.polys[!sapply(int.polys,is.null)] #take out NULL slots from 'int.polys' list
    cat(paste(length(int.polys.nonulls)," overlapping features found. Go to 'dshm_correct_segments' function."))

    inter.Polys<-do.call(raster::bind,int.polys.nonulls)

    return(list(inter.IDs = int.ids, inter.Polys = inter.Polys))

  }
}
