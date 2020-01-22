#wrap function for dshm_split_transects
dshm_split_transect <- function(transect.data, inter.dist, lwr, search.time, w, mute = TRUE, cap){

  x <- floor(sp::SpatialLinesLengths(transect.data)/lwr) #maximum number of segements given the length of the given segment

  #if the transect is less than 2x minimum segment length, it will not be splitted
  if(x<=1){
    #conversion from SpatialLinesDataFrame to SpatialLines
    seg <- raster::spLines(sp::coordinates(transect.data)[[1]][[1]], crs = raster::crs(transect.data))

    #re-conversion to SpatialLinesDataFrame adding original transect label and sample (or segment) label
    seg <- sp::SpatialLinesDataFrame(seg, data.frame(Transect.Label = transect.data$ID, Sample.Label = 1))

    if(cap){

      seg.buf <- rgeos::gBuffer(seg, width = w, capStyle = "ROUND", byid = TRUE) #adding a buffer

    } else {

      seg.buf <- rgeos::gBuffer(seg, width = w, capStyle = "FLAT", byid = TRUE)

    }

    seg.buf$length <- sp::SpatialLinesLengths(seg)/10^3 #adding length information in km
    seg.buf$area <- raster::area(seg.buf)/10^6 #adding area information in km^2

  } else { #if the transect is more than 2x minimum segment length, it will be splitted

    repeat{

      t1 <- proc.time() #starting recording time
      x <- x-1 #number of cutting points (this is always one order smaller than the number of segments)

      repeat{ #starting routine that splits the transect

        #sampling cutting points on the transect in a 'stratified random approach' (note that sometimes it might happen than one point lays out of the line, this will get an error)
        pts <- sp::spsample(transect.data, x, type = "stratified")
        shape_conv <- sf::st_as_sfc(transect.data) #conversion from SpatialLinesDataFrame to 'sfc' object
        pts_conv <- sf::st_as_sfc(pts) #conversion from SpatialPoints to 'sfc' object
        pts_buff <- sf::st_buffer(pts_conv, dist = inter.dist) #adding a buffer to the points
        diff <- sf::st_difference(shape_conv, sf::st_union(sf::st_combine(pts_buff))) #difference between lines and points gives segments
        diff <- sf::st_cast(diff, "LINESTRING") #conversion from 'MULTILINESTRING' to 'LINESTRING'

        coords <- list() #empty vector to store coordinates
        sls <- list() #empty vector to store segments
        sls.buf <- list() #empty vector to store buffered segments
        d <- c() #empty vector to store segment distances

        for(i in 1:length(diff)){ #for each segment in 'diff' object

          xy <- sf::st_coordinates(diff[i])[,1:2] #get coordinates
          coords[[i]] <- xy #store coordinates

          #conversion to SpatialLinesDataFrame
          sls[[i]] <- sp::SpatialLinesDataFrame(raster::spLines(xy, crs = raster::crs(transect.data)),
                                                data = data.frame(Transect.Label = transect.data$ID, Sample.Label = i))

          d[i] <- sp::SpatialLinesLengths(sls[[i]]) #calculate and store distance
          sls.buf[[i]] <- rgeos::gBuffer(sls[[i]], width = w, capStyle = "FLAT", byid = TRUE) #add buffer with 'flat' end

          sls.buf[[i]]$length <- sp::SpatialLinesLengths(sls[[i]])/10^3 #add info on length in m^2
          sls.buf[[i]]$area <- raster::area(sls.buf[[i]])/10^6 #add info on area in km^2

        }

        if (!mute) { #if true=FALSE, then it prints segment distances at every iteration
          print(d)
        }

        t2<-(proc.time() - t1) #time after each iteration
        if(sum(d>lwr) == length(d)||t2 > search.time){ #if the length of all segments meets the input limit OR if the searching time is greater than the input limit, the routine exits and goest to x-1

          break

        }
      }

      if(x==1){ #if the number of cutting points reaches 1, next iteration starts for the maximum again

        x <- floor(sp::SpatialLinesLengths(transect.data)/lwr)

      }

      if(sum(d > lwr)==length(d)){ #if the length of all segments meets the input limit, the routine stops

        break

      }

    }

    seg.buf <- do.call(raster::bind, sls.buf) #all segments is list 'sls.buf' are bound in one single object

    if(cap){

      whole <- rgeos::gBuffer(transect.data, width = w - 1, capStyle = "ROUND", byid = TRUE) #buffered transect. Note the -2 due to imperfection in difference calculation (below)
      whole_seg <- rgeos::gDifference(whole, seg.buf, byid = TRUE, id = as.character(seg.buf$Transect.Label)) #Difference between whole transect and the segments. This is used to extract the capped ends!

      if(length(seg.buf)==2){ #if the segments are only 2, capped ends are used and bound together to get the final segments

        whole_seg <- sp::disaggregate(whole_seg) #separating each component of the difference between transect and segments
        whole_seg$area <- raster::area(whole_seg) #calculating the area
        caps <- whole_seg[whole_seg$area > lwr*w*2,] #selecting those that have an area > lwr*w*2

      } else if(length(seg.buf)==3){ #if the segments are 3, capped ends are used and bound with the middle segment

        caps <- sp::disaggregate(whole_seg[2,]) #capped ends are equivalent of taking out the middle segment

      } else { #if the segments are > 3

        upr.capped <- sp::disaggregate(whole_seg[2,]) #one capped end
        lwr.capped <- sp::disaggregate(whole_seg[length(whole_seg)-1,]) #opposite capped end
        upr.capped$area <- raster::area(upr.capped) #calculating the area
        lwr.capped$area <- raster::area(lwr.capped)
        upr.capped <- upr.capped[upr.capped$area==min(upr.capped$area),] #extracting the segment with the smaller area
        lwr.capped <- lwr.capped[lwr.capped$area==min(lwr.capped$area),]
        caps <- raster::bind(upr.capped, lwr.capped) #binding capped segments in one object

      }

      caps$Transect.Label <- rep(transect.data$ID,2) #adding transect label information to caps
      caps$Sample.Label <- c(1,length(seg.buf)) #adding segment label (i.e. 1 or max length of segment, order does not matter)
      caps$length <- seg.buf[c(1,length(seg.buf)),]$length #adding length information
      caps$area <- raster::area(caps)/10^6 #recalculating area in km^2

      if(length(seg.buf)==2){ #if the segments are 2, caps are bound together

        seg.buf <- caps

      } else { #if segments are > 2 caps are bound with middle segment/s

        seg.buf <- raster::bind(seg.buf[-c(1,length(seg.buf)),],caps)

      }
    }
  }

  return(seg.buf) #splitted transect with capped ends is returned

}

#function that prepares data for dshm_fit
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
    obsdata$p <- stats::predict(det.fn$ddf, newdata = obsdata)$fitted  #predicting p
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

#wrap function for bootstrap
dshm_fit_4boot <- function(det.fn.par, effects.pa,effects.ab, method, lim, distdata, obsdata, segdata, grid, w.pa, w.ab, ID.pa, ID.ab, knots.pa, knots.ab, group, indices,stratification) {

  # for each simulation, random sampling of the rows of original dataset according to the speciefied number
  if(stratification=="none"){
    id<-sample(rownames(distdata), replace = TRUE)
  } else {

    id<-list()
    if(stratification=="stratum"){
      for (j in 1:length(levels(distdata$Region.Label))){
        id[[j]]<-sample(rownames(subset(distdata,distdata$Region.Label==levels(distdata$Region.Label)[j])),replace=TRUE)
      }
    } else if (stratification == "transect") {
      for (j in 1:length(levels(distdata$Transect.Label))){
        id[[j]]<-sample(rownames(subset(distdata,distdata$Transect.Label==levels(distdata$Transect.Label)[j])),replace=TRUE)
      }
    }
    id<-do.call(what = c,id)
  }

  distdata <- distdata[id, ]  #randomly resampling rows in distance dataset with replcement
  distdata$object <- rownames(distdata)  #changing object name according to resampling
  if (group) {
    distdata$size<-1
  }

  det.fn <- Distance::ds(data = distdata, transect = det.fn.par$transect, key = det.fn.par$key, adjustment = det.fn.par$adjustment, truncation = det.fn.par$truncation,
                         formula = det.fn.par$formula, quiet = TRUE)  #fitting detection function

  data <- dshm_prep(det.fn, obsdata, segdata, group)  #calculating detection probabilities, correcting abundance and prepare segment dataset for analysis

  data <- data[indices, ]

  pa <- data$pa  #presence-absence for Cc
  data.ab <- subset(data, pa > 0)  #presence-only abundance dataset for Cc
  ab <- data.ab$abund  #presence-only abundance data for Cc
  ab.full <- data$abund  #full abundance data for Cc

  # HURDLE MODEL - presence absence----

  rowsL <- length(effects.pa)
  mod.sel.pa <- data.frame(mod = rep(0, rowsL), expl_Dev = rep(0, rowsL), Chisq = rep(0, rowsL), logLik = rep(0, rowsL), AICc = rep(0,
                                                                                                                                    rowsL))  #empty dataframe (model selection table for pa data) to fill with model type, % explained deviance, chi-squared test, log-likelihood, and AICc for different covariate combinations

  mod.sel.pa$mod <- effects.pa  #specifying the single main effects as well as their combinations

  mod.list.pa <- list()  #empty list for fitted model in the loop below

  for (i in 1:rowsL) {
    # this loop fits binomial gam for each of the specifyied combinations allowing for different k

    eq <- paste("pa~", effects.pa[i],"+offset(log(area))")  #specifying the equations for the single main effects, the switch function allows for changing the covariates within the loop framework

    gam.pa <- mgcv::gam(stats::as.formula(eq), family = stats::binomial(link = "logit"), data = data, method = method,knots = knots.pa[[i]])  #fitting the gam according to the specified equations

    mod.list.pa[[i]] <- gam.pa  #storing each of the fitted model in the empty list specified before

    mod.sel.pa[i, 2] <- round(((gam.pa$null.deviance - gam.pa$deviance)/gam.pa$null.deviance) * 100, 2)  #calculating the % of deviance explained by each model
    mod.sel.pa[i, 3] <- round((1 - stats::pchisq(gam.pa$deviance, gam.pa$df.residual)), 2)  #Chi-squared statistics (GOF)
    mod.sel.pa[i, 4] <- round(stats::logLik(gam.pa)[1], 3)  #log-likelihood
    mod.sel.pa[i, 5] <- MuMIn::AICc(gam.pa)  #AICc
  }

  mod.sel.pa <- data.frame(mod.sel.pa, deltaAICc = mod.sel.pa$AICc - min(mod.sel.pa$AICc), w = round(MuMIn::Weights(mod.sel.pa$AICc),
                                                                                                           3))  #adding delta AICc and AIC weights to the model selection table
  mod.sel.pa <- mod.sel.pa[order(mod.sel.pa$deltaAICc), ]  #order the model selection table according to increasing delta AICc

  best.pa <- list()  #empty list for the best models in the model selection table
  sub.pa <- subset(mod.sel.pa, mod.sel.pa$w >= lim)  #taking the models which have a delta AIC less or equal to 4

  if (length(ID.pa) > 1) {
    # conditional part that decides to average models if there are more than one model in the list with the best models
    eval.pa <- array(0, dim = c(length(data[, 1]), 1, length(ID.pa)))  #specifying the array with number of rows equal to the length of the full dataset, 1 column, and dimensions equal to the number of models in the list with the best models
    grid.pa <- array(0, dim = c(length(grid[, 1]), 1, length(ID.pa)))
    for (i in 1:length(ID.pa)) {
      # loop going through all models in the list, extracting the fitted values, multiplying them by each model weight, and storing
      # the resulting vectors into the array
      best.pa[i] <- mod.list.pa[ID.pa[i]]
      eval.pa[, , i] <- stats::predict(mod.list.pa[[ID.pa[i]]], newdata = data, type = "response") * (w.pa[ID.pa[i]]/sum(w.pa[ID.pa]))  #multiplying fitted values by weights
      grid.pa[, , i] <- stats::predict(mod.list.pa[[ID.pa[i]]], newdata = grid, type = "response") * (w.pa[ID.pa[i]]/sum(w.pa[ID.pa]))
    }
    fit.w.pa <- apply(eval.pa, 1, sum)  #sum across the array dimension to get final averaged fitted values
    grid.w.pa <- apply(grid.pa, 1, sum)
  } else {
    # conditional part if there is only one model in the list, no averaging
    best.pa <- mod.list.pa[ID.pa[1]]
    fit.w.pa <- stats::predict(mod.list.pa[[ID.pa[1]]], newdata = data, type = "response")  #taking fitted values of the model
    grid.w.pa <- stats::predict(mod.list.pa[[ID.pa[1]]], newdata = grid, type = "response")
  }  #end of conditional part

  # HURDLE MODEL - abundance----

  mod.sel.ab <- data.frame(mod = rep(0, rowsL), expl_Dev = rep(0, rowsL), Chisq = rep(0, rowsL), logLik = rep(0, rowsL), AICc = rep(0,
                                                                                                                                    rowsL))  #empty dataframe (model selection table for ab data) to fill with model type, % explained deviance, chi-squared test, log-likelihood, and AICc for different covariate combinations

  mod.sel.ab$mod <- effects.ab  #specifying the single main effects as well as their combinations

  mod.list.ab <- list()  #empty list for fitted model in the loop below

  for (i in 1:rowsL) {
    # this loop fits gam for each of the specifyied combinations allowing for different distributions and k

    eq <- paste("ab~", effects.ab[i],"+offset(log(phat*area))")  #specifying the equations for the single main effects, the switch function allows for changing the covariates within the loop framework
    gam.ab <- mgcv::gam(stats::as.formula(eq), family = countreg::ztpoisson(), data = data.ab, method = method,knots = knots.ab[[i]])  #fitting the gam according to the specified equations

    mod.list.ab[[i]] <- gam.ab  #storing each of the fitted model in the empty list specified before

    mod.sel.ab[i, 2] <- round(((gam.ab$null.deviance - gam.ab$deviance)/gam.ab$null.deviance) * 100, 2)  #calculating the % of deviance explained by each model
    mod.sel.ab[i, 3] <- round((1 - stats::pchisq(gam.ab$deviance, gam.ab$df.residual)), 2)  #Chi-squared statistics (GOF)
    mod.sel.ab[i, 4] <- round(stats::logLik(gam.ab)[1], 2)  #log-likelihood
    mod.sel.ab[i, 5] <- MuMIn::AICc(gam.ab)  #AICc
  }

  mod.sel.ab <- data.frame(mod.sel.ab, deltaAICc = mod.sel.ab$AICc - min(mod.sel.ab$AICc), w = round(MuMIn::Weights(mod.sel.ab$AICc),
                                                                                                           3))  #adding delta AICc and AIC weights to the model selection table
  mod.sel.ab <- mod.sel.ab[order(mod.sel.ab$deltaAICc), ]  #order the model selection table according to increasing delta AICc

  best.ab <- list()  #empty list for the best models in the model selection table
  sub.ab <- subset(mod.sel.ab, mod.sel.ab$w >= lim)  #taking the models which have a delta AIC less or equal to 4
  if (length(ID.ab) > 1) {
    # conditional part that decides to average models if there are more than one model in the list with the best models
    eval.ab <- array(0, dim = c(length(data.ab[, 1]), 1, length(ID.ab)))  #specifying the array with number of rows equal to the length of the presence-only abundance dataset, 1 column, and dimensions equal to the number of models in the list with the best models
    eval.ab.full <- array(0, dim = c(length(data[, 1]), 1, length(ID.ab)))  #specifying the array with number of rows equal to the length of the full dataset, 1 column, and dimensions equal to the number of models in the list with the best models
    grid.ab <- array(0, dim = c(length(grid[, 1]), 1, length(ID.ab)))
    for (i in 1:length(ID.ab)) {
      # loop going through all models in the list, extracting the fitted values, multiplying them by each model weight, and storing
      # the resulting vectors into the array
      best.ab[[i]] <- mod.list.ab[[ID.ab[i]]]  #storing each best model into a new list
      eval.ab[, , i] <- stats::predict(mod.list.ab[[ID.ab[i]]], newdata = data.ab, type = "response") * (w.ab[ID.ab[i]]/sum(w.ab[ID.ab]))  #multiplying fitted values by weights
      eval.ab.full[, , i] <- stats::predict(mod.list.ab[[ID.ab[i]]], newdata = data, type = "response") * (w.ab[ID.ab[i]]/sum(w.ab[ID.ab]))  #multiplying predicted values on the full dataset by weights
      grid.ab[, , i] <- stats::predict(mod.list.ab[[ID.ab[i]]], newdata = data.frame(grid,phat=1), type = "response") * (w.ab[ID.ab[i]]/sum(w.ab[ID.ab]))
    }
    fit.w.ab <- apply(eval.ab, 1, sum)  #sum across the array dimension to get final averaged fitted values (presence-only abundance)
    fit.w.ab.full <- apply(eval.ab.full, 1, sum)  #sum across the array dimension to get final averaged fitted values (full dataset, Hurdle model evaluation)
    grid.w.ab <- apply(grid.ab, 1, sum)
  } else {
    # conditional part if there is only one model in the list, no averaging
    best.ab <- mod.list.ab[ID.ab[1]]
    fit.w.ab <- stats::predict(mod.list.ab[[ID.ab[1]]], newdata = data.ab, type = "response")  #taking fitted values of the model
    fit.w.ab.full <- stats::predict(mod.list.ab[[ID.ab[1]]], newdata = data, type = "response")  #taking predicted values (full dataset) of the model
    grid.w.ab <- stats::predict(mod.list.ab[[ID.ab[1]]], newdata = data.frame(grid,phat=1), type = "response")
  }  #end of conditional part

  eval.H <- ifelse(ab.full==0,fit.w.pa,fit.w.ab.full * fit.w.pa)  #multiplying pa with ab (full dataset) predictions (HURDLE)
  grid.H <- grid.w.ab * grid.w.pa

  return(c(grid.H, eval.H, ab.full))

}
