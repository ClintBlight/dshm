#' Prepares data for Hurdle model fitting
#'
#' @param det.fn Detection function fitted by ds.
#' @param obsdata Dataframe object containing 4 columns: (1) 'Sample.Label' (i.e. label for segments), (2) 'size' for cluster size, (3) 'distance' (in km) for perpendicular distance of sighting from the transect line, and (4) 'Effort' for segment length (in km).
#' @param segdata Dataframe object with at least 3 columns: (1) 'Transect.Label' (i.e. label for transects), (2) 'Sample.Label' (i.e. label for segments), and (3) 'Effort' for segment length (in km). It may also contain additional columns with relevant habitat covariates specific to each segment that will be fed into the spatial model.
#' @param group If TRUE group abundance is estimated (i.e. 'size' = 1).
#' @param strip.width Transect strip width to calculate segment area if area is not provided.
#' @return Input dataframe 'segdata' with an additional column for animal abundance ('abund') corrected by the detection probability estimated by ds. If they are not included, a column for presence-absence ('pa') is added together with one for segment area (Effort*w, i.e. the strip width estimated by ds).
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}


dshm_prep <- function(det.fn, obsdata, segdata, group=FALSE, strip.width=NULL) {

    if(is.null(det.fn)){
      obsdata$p<-1
      w<-strip.width
    } else {
      obsdata$p <- predict(det.fn$ddf, newdata = obsdata)$fitted  #predicting p
      w <- summary(det.fn)$ddf$meta.data$width  #extracting strip width from ddf object
    }


    if (group) {
      obsdata$size<-1
    }

    #obsdata$size.c <- obsdata$size/obsdata$p  #correcting size with p
    obsdata <- dplyr::group_by(obsdata, Sample.Label)  #grouping by segment
    obsdata <- dplyr::summarize(obsdata, size.c = sum(size, na.rm = TRUE),phat=mean(p,na.rm=TRUE))  #summarizing calculating the sum of individuals for each segment
    obsdata <- as.data.frame(obsdata)  #conversion to dataframe
    #obsdata$size.c <- round(obsdata$size.c)

    segdata$abund <- rep(0)  #adding empty column to segdata for abundance
    segdata$phat <- rep(1)
    for (i in 1:length(obsdata[, 1])) {
        # the loop takes the segment label in obsdata and if it corresponds to any of the segemnt labels in segement data it pastes
        # the corresponding corrected abundance value
        for (j in 1:length(segdata[, 1])) {
            if (obsdata$Sample.Label[i] == segdata$Sample.Label[j]) {
                segdata$abund[j] <- obsdata$size.c[i]
                segdata$phat[j] <- obsdata$phat[i]
            }
        }
    }

    if (sum(colnames(segdata) == "area") == 0) {
        # if there is no area column in the dataset, calculating the area taking the effort and multiplying it by two times the strip
        # width
        possibleError <- tryCatch({segdata$area <- (segdata$Effort * w * 2)}, error = function(e) e)
      if (inherits(possibleError, "error")) {
        stop("Either strip.width not specified (for strip transects) or column for Effort missing (or differently named) in segment.data")
      }
    }

    if (sum(colnames(segdata) == "pa") == 0) {
        # if there is no presence-absence column (pa), calculating it
        segdata$pa <- ifelse(segdata$abund > 0, 1, 0)
    }

    return(segdata)  #the function returns the segdata dataset formatted for modelling

}
