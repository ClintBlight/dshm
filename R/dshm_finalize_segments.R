#' Preparing segments for spatial analysis
#'
#' \code{dshm_finalize_segments} calculates covariate statisitcs within each segment. It also corrects segments for land or other obstacles not relevant to the analysis.
#'
#' @param segment.data Segments as SpatialPolygonsDataFrame. You can use the function \code{dshm_split_transects} to create segments from transect lines.
#' @param land.data A SptialPolygonsDataFrame such as land (for marine species) or other obstacles not relevant to predictions.
#' @param covariates List of covariate raters.
#' @param fun Function for covariate statistcs within each segment (i.e. mean or median).
#' @param parallel If \code{TRUE} the function is run on multiple cores. Default is \code{FALSE}.
#' @param ncores Number of cores if parallel is \code{TRUE}.
#' @return Segments as SptialPolygonsDataFrame with following data attached:
#' \itemize{
#'   \item Transect.Label: ID for split transect.
#'   \item Sample.Label: ID for segment.
#'   \item length: segment length.
#'   \item area: segment area.
#'   \item XYZ covariates: different habitat covariates such as depth, distance to coast, etc. specific to each segment.
#' }
#' @details For more information about splitting transects into segments as well as checking and correcting segments you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{split_transects.pdf} tutorial.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}
#' @export
#'
dshm_finalize_segments<-function(segment.data,land.data = NULL,covariates,fun,parallel=FALSE,ncores=NULL){

  t1<-proc.time() #starts recording time
  cat("\n\n Correcting segments and sampling covariates...\n\n ") #message
  if (parallel=="TRUE") { #parallel execution
    if(Sys.info()[[1]]=="Windows"){
      cl<-parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
    } else {
      cl<-doMC::registerDoMC(ncores) #register cores
    }

    `%dopar%` <- foreach::`%dopar%`

    segments<-foreach::foreach(i=1:length(segment.data)) %dopar% { #running the 'dshm_split_segments' on multiple cores
      if (!is.null(land.data)) {
        ext<-rgeos::gDifference(segment.data[i,],rgeos::gBuffer(land.data, byid=TRUE, width=0))
      } else {
        ext <- segment.data[i,]
      }
      ext$Transect.Label<-segment.data[i,]$Transect.Label
      ext$Sample.Label<-segment.data[i,]$Sample.Label
      ext$length<-segment.data[i,]$length
      ext$area<-raster::area(ext)/10^6
      for (i in 1:length(covariates)){
        ext@data[,i+4]<-raster::extract(covariates[[i]],ext,fun=fun,na.rm=TRUE)[1]
        colnames(ext@data)[i+4]<-paste(names(covariates)[i])
      }
      return(ext) #returning segments
    }

    if(Sys.info()[[1]]=="Windows"){
      parallel::stopCluster(cl)
    }

    t2<-(proc.time()-t1) #stopping recording time
    cat(paste(round(t2[3]/60,3)," minutes elapsed.")) #printing elapsed time
    segments.bind<-do.call(raster::bind, segments) #binding all the segments together in one object
  } else { #non-parallel execution
    `%do%` <- foreach::`%do%`
    pb <- utils::txtProgressBar(min = 0, max = length(segment.data), style = 3) #setting progress bar (not available for parallel)

    segments<-foreach::foreach(j=1:length(segment.data)) %do% { #running the 'dshm_split_segments'
      if (!is.null(land.data)) {
        ext<-rgeos::gDifference(segment.data[j,],rgeos::gBuffer(land.data, byid=TRUE, width=0))
      } else {
        ext <- segment.data[j,]
      }
      ext$Transect.Label<-segment.data[j,]$Transect.Label
      ext$Sample.Label<-segment.data[j,]$Sample.Label
      ext$length<-segment.data[j,]$length
      ext$area<-raster::area(ext)/10^6
      for (i in 1:length(covariates)){
        ext@data[,i+4]<-raster::extract(covariates[[i]],ext,fun=fun,na.rm=TRUE)[1]
        colnames(ext@data)[i+4]<-paste(names(covariates)[i])
      }
      utils::setTxtProgressBar(pb, j) #updating progress bar at each iteration
      return(ext) #returning segments
    }

    t2<-(proc.time()-t1) #stopping recording time
    cat(paste("\n\n ",round(t2[3]/60,3)," minutes elapsed.\n\n ")) #printing elapsed time
    segments.bind<-do.call(raster::bind, segments) #binding all the segments together in one object
  }

  return(segments.bind) #returning bound segments

}
