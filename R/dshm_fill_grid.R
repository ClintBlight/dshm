#' Corrects grid shapefile for land and calculates covariate statistics within each grid cell
#'
#' @param empty.grid Empty grid shapefile.
#' @param land.data Coastline shapefile.
#' @param cov List of covriate rasters.
#' @param fun Function for covariate statistcs within each cell (i.e. mean or median).
#' @param ncores Number of cores (this function is run in parallel only).


#' @export
dshm_fill_grid<-function(empty.grid,land.data = NULL,cov,fun,ncores = 2){

  if (!is.null(land.data)){
    grid.cor <- rgeos::gDifference(empty.grid,raster::union(rgeos::gBuffer(land.data, width=0)),byid = TRUE)
    centr.id<-sp::over(empty.grid,raster::aggregate(grid.cor)) #identify cells in empty grid falling within the corrected grid unified, this gives 1's or NA's
    empty.grid$centr.id<-centr.id #adding info to the empty grid dataframe
    empty.grid.cor<-empty.grid[!is.na(empty.grid$centr.id),] #taking only the 1's
    grid.cor$x.coord<-sp::coordinates(empty.grid.cor)[,1] #adding columns for centroid x's
    grid.cor$y.coord<-sp::coordinates(empty.grid.cor)[,2] #adding column for centroid y's
    grid.cor$area<-raster::area(grid.cor)/10^6 #adding information about area
  } else {
    grid.cor <- empty.grid
    grid.cor$x.coord<-sp::coordinates(grid.cor)[,1] #adding columns for centroid x's
    grid.cor$y.coord<-sp::coordinates(grid.cor)[,2] #adding column for centroid y's
    grid.cor$area<-raster::area(grid.cor)/10^6 #adding information about area
    grid.cor@data<-grid.cor@data[,-1]
  }

  if(Sys.info()[[1]]=="Windows"){
    cl<-parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    cl<-doMC::registerDoMC(ncores) #register cores
  }

  #setting up jobs
  cells<-length(grid.cor)
  m.cells.job<-round(cells/ncores)
  jobs<-list()
  for (i in 1:ncores){

    jobs[[i]]<-m.cells.job
    if (i==ncores){
      jobs[[i]]<-m.cells.job+(cells-(m.cells.job*ncores))
    }

    if (i>1){
      jobs[[i]]<-jobs[[i]]+jobs[[i-1]]
    }
  }

  `%dopar%` <- foreach::`%dopar%`
  grid.cor.cov<-foreach::foreach(i = 1:length(jobs)) %dopar% {

    if(i==1){
      covar.val<-matrix(ncol=length(cov),nrow=jobs[[i]])
      for (j in 1:length(cov)){
        covar.val[,j]<-raster::extract(cov[[j]],grid.cor[1:jobs[[i]],],mean,na.rm=TRUE)
      }

    } else {
      covar.val<-matrix(ncol=length(cov),nrow=(jobs[[i]]-jobs[[i-1]]))
      for (j in 1:length(cov)){
        covar.val[,j]<-raster::extract(cov[[j]],grid.cor[(jobs[[i-1]]+1):jobs[[i]],],mean,na.rm=TRUE)
      }

    }
    return(covar.val)
  }

  if(Sys.info()[[1]]=="Windows"){
    parallel::stopCluster(cl)
  }

  grid.cor.cov<-do.call(grid.cor.cov,what = rbind)

  grid.cor.cov<-as.data.frame(grid.cor.cov)
  names(grid.cor.cov)<-paste(names(cov))
  grid.cor@data[,4:(3+length(cov))]<-grid.cor.cov


  na.id<-data.frame(id=c(1:length(grid.cor)))

  for (i in 1:length(grid.cor)){
    na.id$na[i]<-sum(is.na(grid.cor@data)[i,])
  }

  grid.cor.nona<-grid.cor[subset(na.id,na.id$na==0)$id,]
  grid.cor.nona@data<-data.frame(id=c(1:length(grid.cor.nona)),grid.cor.nona@data)
  rownames(grid.cor.nona@data)<-NULL

  return(grid.cor.nona)
}
