#' Preparing segments for spatial analysis
#'
#' \code{dshm_finalize_segments} calculates covariate statisitcs within each segment. It also corrects segments for land or other obstacles not relevant to the analysis.
#'
#' @param segment.data Segments as SpatialPolygonsDataFrame. You can use the function \code{dshm_split_transects} to create segments from transect lines.
#' @param land.data A SptialPolygonsDataFrame such as land (for marine species) or other obstacles not relevant to predictions.
#' @param covariates List of covariate rasters.
#' @param fun Function for covariate statistcs within each segment (i.e. mean or median).
#' @param parallel If \code{TRUE} the function is running on multiple cores. Default is \code{FALSE}.
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
  
  #covariates needs to be named list
  if (!is.list(covariates) | is.null(names(covariates))) {
    stop("covariates should be a named list")
  }

  #might get better performance by cropping any larger  
  #land.data and raster covariates by region a bit
  #bigger than the extent of the segments.    
  segment_extent <- raster::extent(segment.data)
    
  #expand this extent's width and height by 20%
  segment_extent_width <- (segment_extent@xmax-segment_extent@xmin)
  segment_extent_height <- (segment_extent@ymax-segment_extent@ymin)
  segment_crs_crop_extent <- raster::extent(
    segment_extent@xmin - segment_extent_width*0.1,
    segment_extent@xmax + segment_extent_width*0.1,
    segment_extent@ymin - segment_extent_height*0.1,
    segment_extent@ymax + segment_extent_height*0.1
  )
  
  #make a polygon with several points per side from 
  #expanded extent rectanagle so sides can curve
  #slightly if transformed to another CRS
  segment_crs_crop_extent_poly<-as(segment_crs_crop_extent, className("SpatialPolygons","sp"))  
  segment_crs_crop_extent_line<-as(segment_crs_crop_extent_poly, className("SpatialLines","sp"))
  segment_crs_crop_extent_points<-as(segment_crs_crop_extent_line, className("SpatialPoints","sp"))
  xs <- approx(segment_crs_crop_extent_points@coords[,1],n=33)
  ys <- approx(segment_crs_crop_extent_points@coords[,2],n=33)
  #can now generate the polygon that will be used for any clipping of  
  segment_crs_crop_poly <- raster::spPolygons(cbind(xs$y,ys$y))
  proj4string(segment_crs_crop_poly) <- segment.data@proj4string
  
  #by default set land.data_gbuffered to be NULL 
  land.datais_gbuffered <- NULL
  
  #pre-process any land.data 
  if (!is.null(land.data)) {
    #if the CRS of the land does not match that of the segments
    #generate an appropriate version by clipping and transforming
    if (!(sp::identicalCRS(segment.data,land.data))) {
      cat("\n\n Trying to clip and transform land.data to make a new version with same CRS as segment.data.\n")
      #transform the cropping polygon to same CRS as the land
      land_crs_crop_extent_poly <- sp::spTransform(segment_crs_crop_poly,land.data@proj4string)
      #now crop the land by just the extent of this polygon as should be fast        
      land.data <- raster::crop(land.data,raster::extent(land_crs_crop_extent_poly))
      #lastly transform newly cropped land to same CRS as the segment.data
      land.data <- sp::spTransform(land.data,segment.data@proj4string)
    }
    
    #at this point any land.data should have same CRS as the
    #segments so can crop by the extent generated above  
    land.data_cropped <- raster::crop(land.data,segment_crs_crop_extent)
    #give R option to free memory
    rm(land.data)
    #if there are still land polygons then buffer just once here
    if (!is.null(land.data_cropped)) {
      land.data_gbuffered <- rgeos::gBuffer(land.data_cropped, byid=TRUE, width=0)
    }
    #give R option to free memory
    rm(land.data_cropped)
  }

  #crop the covariates rasters in place in the list for 
  #pre-cropped covariates shouldn't make much difference
  #but for larger rasters should keep memory footprint 
  #down especially for when parallel=TRUE
  cat("\n\n Pre-processing elements in list of covariates.\n")
  for (covariate_name in names(covariates)) {
    #each covariate in the list could have a different CRS
    #so generate an extent of an appropriately transformed 
    #version of the cropping polygon
    crop_extent_covariate_crs <- raster::extent(
      sp::spTransform(
        segment_crs_crop_poly,
        crs(covariates[[covariate_name]])
        )
     )
    #replace the element in the list with a possibly 
    #smaller raster that has been cropped by the extent 
    covariates[[covariate_name]] <- raster::crop(
      covariates[[covariate_name]],
      crop_extent_covariate_crs
      )
  }
  
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
      if (!is.null(land.data_gbuffered)) {
        ext<-rgeos::gDifference(segment.data[i,],land.data_gbuffered)
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
      if (!is.null(land.data_gbuffered)) {
        ext<-rgeos::gDifference(segment.data[j,],land.data_gbuffered)
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
