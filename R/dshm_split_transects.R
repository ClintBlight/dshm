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


dshm_split_transects<-function(transect.data,inter.dist,lwr,search.time,w,parallel=FALSE,ncores=NULL,cap){

  dshm_split_transect<-function(transect.data,inter.dist,lwr,search.time,w,mute=TRUE,cap){

    x<-floor(sp::SpatialLinesLengths(transect.data)/lwr) #maximum number of segements given the minimum length for each segment
    if(x<=1){ #if the transect is less than 2x minimum segment length, it will not be splitted
      seg<-raster::spLines(sp::coordinates(transect.data)[[1]][[1]],crs=raster::crs(transect.data)) #conversion from SpatialLinesDataFrame to SpatialLines
      seg<-sp::SpatialLinesDataFrame(seg,data.frame(Transect.Label=transect.data$ID,Sample.Label=1)) #re-conversion to SpatialLinesDataFrame adding original transect label and sample (or segment) label
      if(cap){
        seg.buf<-rgeos::gBuffer(seg,width=w,capStyle="ROUND",byid=TRUE) #adding a buffer
      } else {
        seg.buf<-rgeos::gBuffer(seg,width=w,capStyle="FLAT",byid=TRUE)
      }
      seg.buf$length<-sp::SpatialLinesLengths(seg)/10^3 #adding length information in km
      seg.buf$area<-raster::area(seg.buf)/10^6 #adding area information in km^2
    } else { #if the transect is more than 2x minimum segment length, it will be splitted
      repeat{
        t1<-proc.time() #starting recording time
        x<-x-1 #number of cutting points (this is always one order smaller than the number of segments)
        repeat{ #starting routine that splits the transect
          pts<-sp::spsample(transect.data,x,type="stratified") #sampling cutting points on the transect in a 'stratified random approach' (note that sometimes it might happen than one point lays out of the line, this will get an error)
          shape_conv<-sf::st_as_sfc(transect.data) #conversion from SpatialLinesDataFrame to 'sfc' object
          pts_conv<-sf::st_as_sfc(pts) #conversion from SpatialPoints to 'sfc' object
          pts_buff<-sf::st_buffer(pts_conv,dist=inter.dist) #adding a buffer to the points
          diff<-sf::st_difference(shape_conv,sf::st_union(sf::st_combine(pts_buff))) #difference between lines and points gives segments
          diff<-sf::st_cast(diff,"LINESTRING") #conversion from 'MULTILINESTRING' to 'LINESTRING'

          coords<-list() #empty vector to store coordinates
          sls<-list() #empty vector to store segments
          sls.buf<-list() #empty vector to store buffered segments
          d<-c() #empty vector to store segment distances
          for(i in 1:length(diff)){ #for each segment in 'diff' object
            xy<-sf::st_coordinates(diff[i])[,1:2] #get coordinates
            coords[[i]]<-xy #store coordinates
            sls[[i]]<-sp::SpatialLinesDataFrame(raster::spLines(xy,crs=raster::crs(transect.data)),data=data.frame(Transect.Label=transect.data$ID,Sample.Label=i)) #conversion to SpatialLinesDataFrame
            d[i]<-sp::SpatialLinesLengths(sls[[i]]) #calculate and store distance
            sls.buf[[i]]<-rgeos::gBuffer(sls[[i]],width=w,capStyle="FLAT",byid=TRUE) #add buffer with 'flat' end
            sls.buf[[i]]$length<-sp::SpatialLinesLengths(sls[[i]])/10^3 #add info on length in m^2
            sls.buf[[i]]$area<-raster::area(sls.buf[[i]])/10^6 #add info on area in km^2
          }

          if (!mute) { #if true=FALSE, then it prints segment distances at every iteration
            print(d)
          }

          t2<-(proc.time()-t1) #time after each iteration
          if(sum(d>lwr)==length(d)||t2>search.time){ #if the length of all segments meets the input limit OR if the searching time is greater than the input limit, the routine exits and goest to x-1
            break
          }
        }
        if(x==1){ #if the number of cutting points reaches 1, next iteration starts for the maximum again
          x<-floor(sp::SpatialLinesLengths(transect.data)/lwr)
        }
        if(sum(d>lwr)==length(d)){ #if the length of all segments meets the input limit, the routine stops
          break
        }
      }
      seg.buf<-do.call(raster::bind, sls.buf) #all segments is list 'sls.buf' are bound in one single object

      if(cap){

        whole<-rgeos::gBuffer(transect.data,width=w-2,capStyle="ROUND",byid=TRUE) #buffered transect. Note the -2 due to imperfection in difference calculation (below)
        whole_seg<-rgeos::gDifference(whole,seg.buf,byid=TRUE,id=as.character(seg.buf$Transect.Label)) #Difference between whole transect and the segments. This is used to extract the capped ends!

        if(length(seg.buf)==2){ #if the segments are only 2, capped ends are used and bound together to get the final segments
          whole_seg<-sp::disaggregate(whole_seg) #separating each component of the difference between transect and segments
          whole_seg$area<-raster::area(whole_seg) #calculating the area
          caps<-whole_seg[whole_seg$area>lwr*w*2,] #selecting those that have an area > lwr*w*2
        } else if(length(seg.buf)==3){ #if the segments are 3, capped ends are used and bound with the middle segment
          caps<-sp::disaggregate(whole_seg[2,]) #capped ends are equivalent of taking out the middle segment
        } else { #if the segments are > 3
          upr.capped<-sp::disaggregate(whole_seg[2,]) #one capped end
          lwr.capped<-sp::disaggregate(whole_seg[length(whole_seg)-1,]) #opposite capped end
          upr.capped$area<-raster::area(upr.capped) #calculating the area
          lwr.capped$area<-raster::area(lwr.capped)
          upr.capped<-upr.capped[upr.capped$area==min(upr.capped$area),] #extracting the segment with the smaller area
          lwr.capped<-lwr.capped[lwr.capped$area==min(lwr.capped$area),]
          caps<-raster::bind(upr.capped,lwr.capped) #binding capped segments in one object
        }

        caps$Transect.Label<-rep(transect.data$ID,2) #adding transect label information to caps
        caps$Sample.Label<-c(1,length(seg.buf)) #adding segment label (i.e. 1 or max length of segment, order does not matter)
        caps$length<-seg.buf[c(1,length(seg.buf)),]$length #adding length information
        caps$area<-raster::area(caps)/10^6 #recalculating area in km^2

        if(length(seg.buf)==2){ #if the segments are 2, caps are bound together
          seg.buf<-caps
        } else { #if segments are > 2 caps are bound with middle segment/s
          seg.buf<-raster::bind(seg.buf[-c(1,length(seg.buf)),],caps)
        }
      }
    }
    return(seg.buf) #splitted transect with capped ends is returned
  }

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


