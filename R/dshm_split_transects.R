#' Splitting transect lines into segments
#'
#' Splits transect lines into segments, adds buffer, and attaches data.
#'
#' @param transect.data Transect lines as SpatialLinesDataFrame.
#' @param inter.dist Distance between each segment in meters. Default is 10 cm.
#' @param lwr Lower limit in meters for segments.
#' @param search.time Time in seconds the routine has to search for the solution. Default is 15s. Note this is not the time required to finish calculations. See the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{split_transects.pdf} tutorial for more information.
#' @param w Segment buffer in meters.
#' @param parallel If \code{TRUE} calculations are carried out on multiple cores. Default is \code{FALSE}.
#' @param ncores Number of cores.
#' @param cap If \code{TRUE} a cap is added to the ends of each splitted transect.
#' @return Segments as SpatialPolygonsDataFrame with attached data:
#' \itemize{
#'   \item Transect.Label: ID for split transect.
#'   \item Sample.Label: ID for segment.
#'   \item lenght: segment length.
#'   \item area: segment area.
#' }
#' @details For more information about splitting transects into segments as well as checking and correcting segments you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{split_transects.pdf} tutorial.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}
#' @export
#'
dshm_split_transects<-function(transect.data,inter.dist = 0.01,lwr,search.time = 15,w,parallel=FALSE,ncores=NULL,cap = TRUE){

  t1<-proc.time() #starts recording time
  cat("\n\n Splitting transects...\n\n ") #message
  if (parallel=="TRUE") { #parallel execution
    if(Sys.info()[[1]]=="Windows"){
      cl<-parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
    } else {
      cl<-doMC::registerDoMC(ncores) #register cores
    }

    `%dopar%` <- foreach::`%dopar%`

    segments<-foreach::foreach(i=1:length(transect.data)) %dopar% { #running the 'dshm_split_segments' on multiple cores
      ext<-dshm_split_transect(transect.data[i,],inter.dist,lwr,search.time,w,cap=cap)
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
    pb <- utils::txtProgressBar(min = 0, max = length(transect.data), style = 3) #setting progress bar (not available for parallel)

    segments<-foreach::foreach(i=1:length(transect.data)) %do% { #running the 'dshm_split_segments'
      ext<-dshm_split_transect(transect.data[i,],inter.dist,lwr,search.time,w,cap=cap)
      utils::setTxtProgressBar(pb, i) #updating progress bar at each iteration
      return(ext) #returning segments
    }

    t2<-(proc.time()-t1) #stopping recording time
    cat(paste("\n\n ",round(t2[3]/60,3)," minutes elapsed.\n\n ")) #printing elapsed time
    segments.bind<-do.call(raster::bind, segments) #binding all the segments together in one object
  }

  return(segments.bind) #returning bound segments

}


