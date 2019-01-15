#' Creating a prediction grid
#'
#' \code{dshm_make_grid} creates a prediction grid that can be used for Hurdle model spatial predictions.
#'
#' @param extent Grid extent as specified in \code{\link[raster]{extent}}.
#' @param cell.size Grid cell size in meters.
#' @param projection Grid projection as specified in \code{\link[raster]{crs}}.
#' @return A grid as SpatialPolygonsDataFrame with \code{id} for each grid cell.
#' @details For more information about creating and preparing a prediction grid you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{build_grid.pdf} tutorial.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}
#' @export
#'
dshm_make_grid<-function(extent,cell.size,projection){
  if(sign(extent[1])*sign(extent[2])==1){
    n.cells.x<-floor(abs(abs(extent[2])-abs(extent[1]))/cell.size)
  } else {
    n.cells.x<-floor(abs(abs(extent[2])+abs(extent[1]))/cell.size)
  }

  if(sign(extent[3])*sign(extent[4])==1){
    n.cells.y<-floor(abs(abs(extent[4])-abs(extent[3]))/cell.size)
  } else {
    n.cells.y<-floor(abs(abs(extent[4])+abs(extent[3]))/cell.size)
  }

  topo<-sp::GridTopology(c(extent[1],extent[3])+cell.size/2, c(cell.size,cell.size), c(n.cells.x,n.cells.y))
  grid<-sp::SpatialGrid(topo,proj4string = projection)
  grid.data<-sp::SpatialGridDataFrame(grid,data=data.frame(id=c(1:length(grid))))
  grid.poly<-inlmisc::Grid2Polygons(grid.data,zcol = "id",level = FALSE)
  names(grid.poly)<-c("id")
  return(grid.poly)
}
