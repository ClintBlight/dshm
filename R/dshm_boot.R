#' Non-parametric bootstrap for Hurdle model uncertainty
#'
#' \code{dshm_boot} performs a non-parametric bootstrap to generate Hurdle model predictions on a spatial grid. Prediction grids can then be used to calculate confidence intervals. The function bases on a (stratified) sampling process with replacement at the segment level.
#'
#' @param det.fn.par List of detection function parameters. For strucuture see the documentation for \code{\link[Distance]{ds}}.
#' @param effects.pa List of characters defining the binomial gam models to be fitted. For model structure see \code{\link[mgcv]{gam}}.
#' @param effects.ab List of characters defining the zero-truncated Poisson gam models to be fitted. For model structure see \code{\link[mgcv]{gam}}.
#' @param distdata Dataframe for distance sampling observations. For strucuture see the documentation for \code{\link[Distance]{ds}}.
#' @param obsdata Dataframe object with the following structure:
#' \itemize{
#'   \item Region.Label: ID for stratum where the animal was observed.
#'   \item Transect.Label: ID for transect where the animal was observed.
#'   \item Sample.Label: ID for segment where the animal was observed.
#'   \item distance: sighting perpendicular distance from the transect line.
#'   \item size: sighting size, i.e. number of animals.
#'   \item object: sighting ID.
#' }
#' @param segdata Dataframe object with the following strucuture:
#' \itemize{
#'   \item Region.Label: ID for stratum where the transects and segments are located.
#'   \item Transect.Label: ID for split transect.
#'   \item Sample.Label: ID for segment.
#'   \item length: segment length.
#'   \item area: segment area.
#'   \item XYZ covariates: different habitat covariates such as depth, distance to coast, etc. specific to each segment.
#' }
#' You do not have to create segdata manually. You can use the functions in \code{\link{dshm}} to automatically split transects into segments. For more information you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{split_transects.pdf} tutorial.
#' @param grid Grid used for model prediction. Column names for habitat covriates should correspond to those in 'segdata'.
#' @param model_fit Model fitted with the function \code{\link{dshm_fit}}.
#' @param group If \code{TRUE} group abundance is estimated (i.e. sighting size = 1). Default is \code{FALSE}.
#' @param nsim Number of simulations.
#' @param parallel If \code{TRUE} the simulations are performed on multiple cores. Default is \code{FALSE}.
#' @param ncores Number of cores for parallel execution.
#' @param mute If \code{TRUE} all unrelevant messages are suppressed. Default is \code{TRUE}.
#' @param stratification Bootstrap can be executed at the level of the \code{"transect"} or \code{"stratum"}. Default is \code{stratification = "none"}.
#' @return A list of two arrays:
#' \itemize{
#'   \item sim_grid: simulated grids.
#'   \item obs_fit: observed and fitted values for each simulation.
#' }
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}
#' @export
#'
dshm_boot <- function(det.fn.par, effects.pa = NULL,effects.ab = NULL, distdata, obsdata, segdata, model_fit, grid, group=FALSE, nsim, parallel = FALSE, ncores = NULL, mute = TRUE, stratification="none") {

    # Defining resampling indices (sampling with replacement)----
    ids <- list()  #empty list
    for (i in 1:nsim) {
      # for each simulation, random sampling of the rows of original dataset according to the speciefied number
      if (stratification == "none"){
        ids[[i]] <- sample(rownames(segdata), replace = TRUE)
      } else {
        id<-list()
        if(stratification=="stratum"){
          for (j in 1:length(levels(segdata$Region.Label))){
            id[[j]]<-sample(rownames(subset(segdata,segdata$Region.Label==levels(segdata$Region.Label)[j])),replace=TRUE)
          }
        } else if(stratification=="transect") {
          for (j in 1:length(levels(segdata$Transect.Label))){
            id[[j]]<-sample(rownames(subset(segdata,segdata$Transect.Label==levels(segdata$Transect.Label)[j])),replace=TRUE)
          }
        }
        id<-do.call(what = c,id)
        ids[[i]] <- id  #random sampling and transformation from chr to num
      }
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
            possibleError <- tryCatch({mod <- dshm_fit_4boot(det.fn.par, effects.pa,effects.ab, method = model_fit$info$method, lim = model_fit$info$lim, distdata, obsdata, segdata, grid, w.pa = model_fit$info$weights.pa, w.ab = model_fit$info$weights.ab, ID.pa = model_fit$info$ID.pa, ID.ab = model_fit$info$ID.ab, knots.pa = model_fit$info$k.loc.pa, knots.ab = model_fit$info$k.loc.ab, group, indices = ids[[i]],stratification)}, error = function(e) e,warning=function(w) w)
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
            possibleError <- tryCatch({mod <- dshm_fit_4boot(det.fn.par, effects.pa,effects.ab, method = model_fit$info$method, lim = model_fit$info$lim, distdata, obsdata, segdata, grid, w.pa = model_fit$info$weights.pa, w.ab = model_fit$info$weights.ab, ID.pa = model_fit$info$ID.pa, ID.ab = model_fit$info$ID.ab, knots.pa = model_fit$info$k.loc.pa, knots.ab = model_fit$info$k.loc.ab, group, indices = ids[[i]],stratification)}, error = function(e) e,warning=function(w) w)
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




