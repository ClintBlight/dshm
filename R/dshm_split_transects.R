#' Splits transect lines into segments and adds buffer as well as metadata
#'
#' @param transect.data Transect lines as SpatialLinesDataFrame.
#' @param inter.dist Distance between each segment in meters.
#' @param lwr Lower limit in meters for segments.
#' @param search.time Time the routine has to search for each set of cutting points.
#' @param w Buffer in meters.
#' @param parallel If TRUE the rountine is carried out by multiple cores. Default is FALSE.
#' @param ncores Number of cores.
#' @param cap If TRUE a cap is added to the ends of each splitted transect.
#' @return Splitted transect lines with buffer and information about: (1) 'Transect.Label' the original transect id, (2) 'Sample.Label' the segment id, (3) 'lenght' the length of the segment, and (4) 'area' the area of the segment.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}

#' @export
dshm_split_transects<-function(transect.data,inter.dist,lwr,search.time,w,parallel=FALSE,ncores=NULL,cap){

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
    pb <- txtProgressBar(min = 0, max = length(transect.data), style = 3) #setting progress bar (not available for parallel)

    segments<-foreach::foreach(i=1:length(transect.data)) %do% { #running the 'dshm_split_segments'
      ext<-dshm_split_transect(transect.data[i,],inter.dist,lwr,search.time,w,cap=cap)
      setTxtProgressBar(pb, i) #updating progress bar at each iteration
      return(ext) #returning segments
    }

    t2<-(proc.time()-t1) #stopping recording time
    cat(paste("\n\n ",round(t2[3]/60,3)," minutes elapsed.\n\n ")) #printing elapsed time
    segments.bind<-do.call(raster::bind, segments) #binding all the segments together in one object
  }

  return(segments.bind) #returning bound segments

}


