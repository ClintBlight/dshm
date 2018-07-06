#' Corrects grid shapefile for land and calculates covariate statistics within each grid cell
#'
#' @param empty.grid Empty grid shapefile.
#' @param land.data Coastline shapefile.
#' @param cov List of covriate rasters.
#' @param fun Function for covariate statistcs within each cell (i.e. mean or median).
#' @param ncores Number of cores (this function is run in parallel only).
#' @param saveRaster If TRUE the corrected and filled grid is saved as a shapefile. Default is FALSE.
#' @param names_simplfied A vector of simplified names for each covariate. Long names might produce problems when saving the grid as a shapefile.
#' @param file_name Name of the saved grid.


dshm_fill_grid<-function(empty.grid,land.data,cov,fun,ncores,saveRaster=FALSE,names_simplfied,file_name){
  grid.cor<-rgeos::gDifference(empty.grid,raster::union(rgeos::gBuffer(land.data, width=0)),byid = TRUE)
  centr.id<-sp::over(empty.grid,raster::aggregate(grid.cor)) #identify cells in empty grid falling within the corrected grid unified, this gives 1's or NA's
  empty.grid$centr.id<-centr.id #adding info to the empty grid dataframe
  empty.grid.cor<-empty.grid[!is.na(empty.grid$centr.id),] #taking only the 1's
  grid.cor$x.coord<-sp::coordinates(empty.grid.cor)[,1] #adding columns for centroid x's
  grid.cor$y.coord<-sp::coordinates(empty.grid.cor)[,2] #adding column for centroid y's
  grid.cor$area<-raster::area(grid.cor)/10^6 #adding information about area

  if(Sys.info()[[1]]=="Windows"){
    cl<-parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    cl<-doMC::registerDoMC(ncores) #register cores
  }

  `%dopar%` <- foreach::`%dopar%`
  grid.cor.cov<-foreach::foreach(i = 1:length(cov),.combine = cbind) %dopar% {
    covar.val<-raster::extract(cov[[i]],grid.cor,fun,na.rm=TRUE)[,1]
    return(covar.val)
  }

  if(Sys.info()[[1]]=="Windows"){
    parallel::stopCluster(cl)
  }

  grid.cor.cov<-as.data.frame(grid.cor.cov)
  names(grid.cor.cov)<-paste(names(cov))
  grid.cor@data[,4:(3+length(cov))]<-grid.cor.cov

  na.id<-data.frame(id=c(1:length(grid.cor)))

  for (i in 1:length(grid.cor)){
    na.id$na[i]<-sum(is.na(grid.cor@data)[i,])
  }

  grid.cor.nona<-grid.cor[subset(na.id,na==0)$id,]
  grid.cor.nona$id<-c(1:length(grid.cor.nona))

  if (saveRaster){
    grid_4save<-grid.cor.nona
    names(grid_4save)<-c("x.coord","y.coord","a",names_simplfied,"id")
    rgdal::writeOGR(obj=grid_4save,dsn="/Users/User/Desktop/data",layer=file_name,driver="ESRI Shapefile",overwrite_layer=TRUE)
  }

  return(grid.cor.nona)
}
