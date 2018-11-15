#' Non-parametric bootstrap for Hurdle model uncertainty
#'
#' @param det.fn.par Detection function parameters. For strucuture see the documentation for the 'ds' package.
#' @param effects.pa List of characters defining the binomial gam models to be fitted. For model structure see \code{\link[mgcv]{gam}}.
#' @param effects.ab List of characters defining the zero-truncated Poisson gam models to be fitted. For model structure see \code{\link[mgcv]{gam}}.
#' @param method Character. Fitting method for gam as in the 'mgcv' package. Note that 'REML' is not available since it is not campatible with the zero-truncated poisson in the 'countreg' package. Default is 'GCV.Cp'.
#' @param lim AIC weight (AICw) threshold for model averaging. Models with AICw < lim are not averaged. Default is 0.1.
#' @param distdata Dataframe for distance sampling observations. For strucuture see the documentation for the 'ds' package.
#' @param obsdata Dataframe object containing 4 columns: (1) 'Sample.Label' (i.e. label for segments), (2) 'size' for cluster size, (3) 'distance' (in km) for perpendicular distance of sighting from the transect line, and (4) 'Effort' for segment length (in km).
#' @param segdata Dataframe object with at least 3 columns: (1) 'Transect.Label' (i.e. label for transects), (2) 'Sample.Label' (i.e. label for segments), and (3) 'Effort' for segment length (in km). It may also contain additional columns with relevant habitat covariates specific to each segment that will be fed into the spatial model.
#' @param grid Grid used for model prediction. Column names for habitat covriates should correspond to those in 'segdata'.
#' @param w.pa Presence-absence submodel variant weights obtained from dshm_fit$info$weights.pa.
#' @param w.ab Weights for submodel variants for abundance conditional on presence obtained from dshm_fit$info$weights.ab.
#' @param ID.pa IDs for best presence-absence submodel variants obtained from dshm_fit$models$pa.
#' @param ID.ab IDs for best submodel variants for abundance conditional on presence obtained from dshm_fit$models$ab.
#' @param knots.pa List of knot gam knot positions for each smooth term of the fitted binomial models.
#' @param knots.ab List of knot gam knot positions for each smooth term of the fitted zero-truncated Poisson models.
#' @param group If TRUE group abundance is estimated (i.e. 'size' = 1).
#' @param nsim Number of simulations,
#' @param parallel If TRUE the simulations are performed on multiple cpus. Default is FALSE.
#' @param ncores Number of cpus for parallel execution.
#' @param mute If TRUE all unrelevant messages are suppressed. Default is TRUE.
#' @param stratification 'transect' or 'stratum'
#' @return A list of two arrays: (1) 'sim_grid' containing all simulated grids and (2) 'obs_fit' containing observed and fitted values for each simulation.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}

#' @export
dshm_boot <- function(det.fn.par, effects.pa = NULL,effects.ab = NULL, method = "GCV.Cp", lim = 0.1, distdata, obsdata, segdata, grid, w.pa, w.ab, ID.pa, ID.ab, knots.pa, knots.ab, group=FALSE, nsim, parallel = FALSE, ncores = NULL, mute = TRUE, stratification="transect") {

    # Defining resampling indices (sampling with replacement)----
    ids <- list()  #empty list
    for (i in 1:nsim) {
      # for each simulation, random sampling of the rows of original dataset according to the speciefied number
      id<-list()

      if(stratification=="stratum"){
        for (j in 1:length(levels(cor.seg.def$Region.Label))){
          id[[j]]<-sample(rownames(subset(cor.seg.def@data,Region.Label==levels(cor.seg.def$Region.Label)[j])),replace=TRUE)
        }
      } else if(stratification=="transect") {
        for (j in 1:length(levels(cor.seg.def$Transect.Label))){
          id[[j]]<-sample(rownames(subset(cor.seg.def@data,Transect.Label==levels(cor.seg.def$Transect.Label)[j])),replace=TRUE)
        }
      }
      id<-do.call(what = c,id)

      ids[[i]] <- id  #random sampling and transformation from chr to num
    }

    # Bootstrapping----
    t1 <- proc.time()  #starting time recording
    cat("\n\n Generating simulations...\n\n ")  #printing that the sim routine is starting
    if (parallel == "TRUE") {
        # parallel routine foreach in parallel package alloes for non- and parallel execution. The function feeds the
        # model.selection.simplified with the previously resampled rows and stores the values within a 2-dimension array. The routine
        # in embedded within a tryCatch function looking for possible fitting errors due to more parameters than available data. In
        # case of error 'Fitting error' is pasted within the array.

      if(Sys.info()[[1]]=="Windows"){
        cl<-parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
      } else {
        cl<-doMC::registerDoMC(ncores) #register cores
      }

        `%dopar%` <- foreach::`%dopar%`
        data.array <- foreach::foreach(i = 1:nsim, .combine = cbind) %dopar% {
            possibleError <- tryCatch({mod <- dshm_fit_4boot(det.fn.par, effects.pa,effects.ab, method, lim, distdata, obsdata, segdata, grid, w.pa, w.ab, ID.pa, ID.ab, knots.pa, knots.ab, group, indices = ids[[i]],stratification)}, error = function(e) e,warning=function(w) w)
            if (inherits(possibleError, "error")|inherits(possibleError, "warning")) {
                return("Fitting error")
            }
            if (inherits(possibleError, "error")&inherits(possibleError, "warning")) {
              return("Fitting error")
            }
            return(mod)
        }

        if(Sys.info()[[1]]=="Windows"){
          parallel::stopCluster(cl)
        }

    } else {
        # same thing as before but without parallelization using %do% insted of %dopar%
        `%do%` <- foreach::`%do%`
        pb <- utils::txtProgressBar(min = 0, max = nsim, style = 3)  #adding a progress bar (only possible without parallelization)
        data.array <- foreach::foreach(i = 1:nsim, .combine = cbind) %do% {
            utils::setTxtProgressBar(pb, i)
            possibleError <- tryCatch({mod <- dshm_fit_4boot(det.fn.par, effects.pa,effects.ab, method, lim, distdata, obsdata, segdata, grid, w.pa, w.ab, ID.pa, ID.ab, knots.pa, knots.ab, group, indices = ids[[i]],stratification)}, error = function(e) e,warning=function(w) w)
            if (inherits(possibleError, "error")|inherits(possibleError, "warning")) {
                return("Fitting error")
            }
            if (inherits(possibleError, "error")&inherits(possibleError, "warning")) {
              return("Fitting error")
            }
            return(mod)
        }
        close(pb)
    }
    t2 <- (proc.time() - t1)  #stopping recording time
    cat(paste(round(t2[[3]]/60, 3), " minutes elapsed.\n\n "))  #printing elapsed time in minutes

    # How many simulation are available (i.e. without fitting problems) and in case of errors reducing the array to the successful
    # iterations----
    new.array <- array(rep(0), dim = c(length(grid[, 1]) + length(segdata[, 1]) * 2, 1, nsim))  #empty array
    if (mute) {
        # option to suppress warnings due to coercion to NA when pasting error messages from data.array to new.array (below)
        options(warn = -1)  #suppressing warnings
    }
    for (i in 1:nsim) {
        # pasting simulations to new.array
        new.array[, , i] <- as.numeric(data.array[, i])
    }
    options(warn = 0)  #enabling warnings again
    d.array.id <- data.frame(id = is.na(new.array[1, , ]))  #logical TRUE if NA are present in each simulation
    d.array.id <- as.numeric(rownames(subset(d.array.id, id == "FALSE")))  #extracting the simulations that were successful
    sim.available <- length(d.array.id)  #counting successful simulations
    cat(paste(sim.available, " simulations are available"))  #printing the number of successful simulations
    decision1 <- readline("Would you like to continue (Yes/No)?\n ")  #user can choose to stop or continue
    if (decision1 == "No") {
        # if user decides not to continue, then stop function and print error message
        stop("User decided to stop", call. = FALSE)
    }
    new.array <- new.array[, , d.array.id]  #if user decides to continue then reduce the array to the successful simulations

    # Possibility to choose the maximum number of animals per grid cells (and further reducing the array) and splitting of
    # new.array into a P.array containing grid predictions and R.array containing 2 columns of the length of the reduced dataset
    # and containing fitted and observed values for further R or var calculation----
    P.array <- array(rep(0), dim = c(length(grid[, 1]), 1, sim.available))  #empty P.array
    R.array <- array(rep(0), dim = c(length(segdata[, 1]), 2, sim.available))  #empty R.array
    decision2 <- readline("Would you like to choose a limit for the amount of animals per grid cell (Yes/No)?\n ")  #user can decide to set a limit of animals per grid cells. Sometimes bootstrapping produces unrealistic predictions with some grid cells having a huge abundance.
    if (decision2 == "Yes") {
        # user decides to select a limit
        amount <- readline("Choose a limit for the amount of animals per grid cell: ")  #enter a value for the limit
        for (i in 1:sim.available) {
            # the loop goes into each of the successful simulations and looks if there is at least one grid value above the chosen limit.
            # In case of detection, it goes to the next simulations. In case of non detection it splits the new.array into two arrays P
            # and R.
            if (sum(new.array[1:length(grid[, 1]), i] > as.numeric(amount)) > 1) {
                next
            } else {
                P.array[, , i] <- new.array[1:length(grid[, 1]), i]
                R.array[, 1, i] <- new.array[(length(grid[, 1]) + 1):(length(grid[, 1]) + length(segdata[, 1])), i]
                R.array[, 2, i] <- new.array[(length(grid[, 1]) + length(segdata[, 1]) + 1):(length(grid[, 1]) + length(segdata[, 1]) *
                  2), i]
            }
        }
        P.array.id <- data.frame(id = P.array[1, 1, ] > 0)  #logical TRUE if the first values of each simulation is 0 (i.e. skipped by the previous loop)
        P.array.id <- as.numeric(rownames(subset(P.array.id, id == "TRUE")))  #extracting the row names for TRUE

        cat(paste(length(P.array.id), " simulations are available"))  #printing the number of available simulations
        decision <- readline("Would you like to continue (Yes/No)?\n ")  #user can decide to continue or stop
        if (decision == "No") {
            # if user decides to stop then end routine and print error message
            stop("User decided to stop", call. = FALSE)
        }
        if (length(P.array.id) < sim.available) {
            # if user decides to continue then select a number of simulations
            newSim <- readline("Choose a number of simulations: ")  #select the number of simulations
            P.array <- P.array[, , sample(P.array.id, newSim)]  #reducing P.array randomly sampling simulations among the available ones
            R.array <- R.array[, , sample(P.array.id, newSim)]  #reducing R.array randomly sampling simulations among the available ones
        }
    } else {
        for (i in 1:sim.available) {
            # part if user decides not to select a limit of animals per grid cell. P and R array are created without further reducing
            # simulation number
            P.array[, , i] <- new.array[1:length(grid[, 1]), i]
            R.array[, 1, i] <- new.array[(length(grid[, 1]) + 1):(length(grid[, 1]) + length(segdata[, 1])), i]
            R.array[, 2, i] <- new.array[(length(grid[, 1]) + length(segdata[, 1]) + 1):(length(grid[, 1]) + length(segdata[, 1]) *
                2), i]
        }
    }

    return(list(sim_grids=P.array,obs_fit=R.array))

}




