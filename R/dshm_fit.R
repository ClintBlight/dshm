#' Fits and averages spatial Hurdle models
#'
#' @param det.fn Detection function fitted by \code{\link[Distance]{ds}}.
#' @param effects.pa List of characters defining the binomial gam models to be fitted. For model structure see \code{\link[mgcv]{gam}}.
#' @param effects.ab List of characters defining the zero-truncated Poisson gam models to be fitted. For model structure see \code{\link[mgcv]{gam}}.
#' @param knots.pa List of knot gam knot positions for each smooth term of the fitted binomial models.
#' @param knots.ab List of knot gam knot positions for each smooth term of the fitted zero-truncated Poisson models.
#' @param method Character. GAM fitting method. Note that 'REML' is not available since it is not campatible with \code{\link[countreg]{ztpoisson}}. Default is 'GCV.Cp'.
#' @param lim AIC weight (AICw) threshold for model averaging. Models with AICw < lim are not averaged. Default is 0.1.
#' @param obsdata Dataframe object containing 4 columns: (1) 'Sample.Label' (i.e. label for segments), (2) 'size' for cluster size, (3) 'distance' (in km) for perpendicular distance of sighting from the transect line, and (4) 'Effort' for segment length (in km).
#' @param segdata Dataframe object with at least 3 columns: (1) 'Transect.Label' (i.e. label for transects), (2) 'Sample.Label' (i.e. label for segments), and (3) 'Effort' for segment length (in km). It may also contain additional columns with relevant habitat covariates specific to each segment that will be fed into the spatial model.
#' @param grid Grid used for model prediction. Column names for habitat covriates should correspond to those in 'segdata'.
#' @param SelectionTable If TRUE model selection table is reported for each submodel. Default is TRUE.
#' @param showSelectedModels If TRUE best fitted submodel variants are reported. Default is FALSE.
#' @param group If TRUE group abundance is estimated (i.e. 'size' = 1).
#' @param strip.width Strip width to calculate segment area if there is no "area" column.
#' @return A list of 3 objects: (1) 'models' containing a list for all presence-absence submodel variants ('pa') and all submodel variants for abundance conditional on presence ('ab'); (2) 'info' containing a list of variant weights and best variant ids for both submodels ('weights.pa' and 'weights.ab', and 'ID.pa' and 'ID.ab' respectively); and (3) 'grid_data' containing a list of prediction grids for presence-absence submodel ('pa'), submodel for abundance conditional on presence ('ab'), and Hurdle model ('H').
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}

#' @export
dshm_fit <- function(det.fn, effects.pa = NULL,effects.ab = NULL,knots.pa=NULL ,knots.ab=NULL, method = "GCV.Cp", lim = 0.1, obsdata,
    segdata, grid, SelectionTable = TRUE, showSelectedModels = FALSE, group = FALSE, strip.width = NULL) {

    data <- dshm_prep(det.fn, obsdata, segdata, group, strip.width)

    pa <- data$pa  #presence-absence for Cc
    data.ab <- subset(data, pa > 0)  #presence-only abundance dataset for Cc
    ab <- data.ab$abund  #presence-only abundance data for Cc
    ab.full <- data$abund  #full abundance data for Cc


    # HURDLE MODEL - presence absence----

    rowsL <- length(effects.pa)
    mod.sel.pa <- data.frame(mod = rep(0, rowsL), k = rep(0, rowsL), expl_Dev = rep(0, rowsL), logLik = rep(0, rowsL), AICc = rep(0,
        rowsL))  #empty dataframe (model selection table for pa data) to fill with model type, % explained deviance, chi-squared test, log-likelihood, and AICc for different covariate combinations

    mod.sel.pa$mod <- effects.pa  #specifying the single main effects as well as their combinations

    mod.list.pa <- list()  #empty list for fitted model in the loop below
    knots.pa.list<-list()
    edf.pa.list<-list()

    for (i in 1:rowsL) {
        # this loop fits binomial gam for each of the specifyied combinations allowing for different k

        eq <- paste("pa~", effects.pa[i],"+offset(log(area))")  #specifying the equations for the single main effects, the switch function allows for changing the covariates within the loop framework

        gam.pa <- mgcv::gam(stats::as.formula(eq), family = stats::binomial(link = "logit"), data = data, method = method,knots = knots.pa[[i]])  #fitting the gam according to the specified equations

        mod.list.pa[[i]] <- gam.pa  #storing each of the fitted model in the empty list specified before
        edf.pa.list[[i]] <- summary(gam.pa)$s.table[,1]

        if(!is.null(gam.pa$smooth[[1]]$xp)){
          knots.smooth.pa<-list()
            for (j in 1:length(gam.pa$smooth)){
              knots.smooth.pa[[j]]<-gam.pa$smooth[[j]]$xp
            }
          smooth.names.pa<-c()
            for (j in 1:length(knots.smooth.pa)){
              smooth.names.pa[j]<-gam.pa$smooth[[j]]$term
            }
          names(knots.smooth.pa)<-smooth.names.pa
          knots.pa.list[[i]]<-knots.smooth.pa
        }

        mod.sel.pa[i, 2] <- round(sum(gam.pa$edf)-1,2)
        mod.sel.pa[i, 3] <- round(((gam.pa$null.deviance - gam.pa$deviance)/gam.pa$null.deviance) * 100, 2)  #calculating the % of deviance explained by each model
        mod.sel.pa[i, 4] <- round(stats::logLik(gam.pa)[1], 2)  #log-likelihood
        mod.sel.pa[i, 5] <- MuMIn::AICc(gam.pa)  #AICc
    }

    mod.sel.pa <- data.frame(mod.sel.pa, deltaAICc = mod.sel.pa$AICc - min(mod.sel.pa$AICc), w = round(MuMIn::Weights(mod.sel.pa$AICc),
        3))  #adding delta AICc and AIC weights to the model selection table
    exdev.pa <- mod.sel.pa$expl_Dev
    w.pa <- mod.sel.pa$w
    k.pa <- mod.sel.pa$k
    mod.sel.pa <- mod.sel.pa[order(mod.sel.pa$deltaAICc), ]  #order the model selection table according to increasing delta AICc

    best.pa <- list()  #empty list for the best models in the model selection table
    sub.pa <- subset(mod.sel.pa, w >= lim)  #taking the models which have a delta AIC less or equal to 4
    ID.pa <- as.integer(rownames(sub.pa))  #extracting the row IDs for each of the best models

    if (length(ID.pa) > 1) {
        # conditional part that decides to average models if there are more than one model in the list with the best models
        eval.pa <- array(0, dim = c(length(data[, 1]), 1, length(ID.pa)))  #specifying the array with number of rows equal to the length of the full dataset, 1 column, and dimensions equal to the number of models in the list with the best models
        grid.pa <- array(0, dim = c(length(grid[, 1]), 1, length(ID.pa)))
        for (i in 1:length(ID.pa)) {
            # loop going through all models in the list, extracting the fitted values, multiplying them by each model weight, and storing
            # the resulting vectors into the array
            best.pa[[i]] <- mod.list.pa[[ID.pa[i]]]
            eval.pa[, , i] <- stats::predict(mod.list.pa[[ID.pa[i]]], newdata = data, type = "response") * (w.pa[ID.pa[i]]/sum(w.pa[ID.pa]))  #multiplying fitted values by weights
            grid.pa[, , i] <- stats::predict(mod.list.pa[[ID.pa[i]]], newdata = grid, type = "response") * (w.pa[ID.pa[i]]/sum(w.pa[ID.pa]))
        }
        fit.w.pa <- apply(eval.pa, 1, sum)  #sum across the array dimension to get final averaged fitted values
        grid.w.pa <- apply(grid.pa, 1, sum)
    } else {
        # conditional part if there is only one model in the list, no averaging
        best.pa <- mod.list.pa[[ID.pa[1]]]
        fit.w.pa <- stats::predict(mod.list.pa[[ID.pa[1]]], newdata = data, type = "response")  #taking fitted values of the model
        grid.w.pa <- stats::predict(mod.list.pa[[ID.pa[1]]], newdata = grid, type = "response")
    }  #end of conditional part


    # HURDLE MODEL - abundance----

    mod.sel.ab <- data.frame(mod = rep(0, rowsL), k = rep(0, rowsL), expl_Dev = rep(0, rowsL), logLik = rep(0, rowsL), AICc = rep(0,
        rowsL))  #empty dataframe (model selection table for ab data) to fill with model type, % explained deviance, chi-squared test, log-likelihood, and AICc for different covariate combinations

    mod.sel.ab$mod <- effects.ab  #specifying the single main effects as well as their combinations

    mod.list.ab <- list()  #empty list for fitted model in the loop below
    knots.ab.list<-list()
    edf.ab.list<-list()

    for (i in 1:rowsL) {
        # this loop fits gam for each of the specifyied combinations allowing for different distributions and k

        eq <- paste("ab~", effects.ab[i],"+offset(log(phat*area))")  #specifying the equations for the single main effects, the switch function allows for changing the covariates within the loop framework
        gam.ab <- mgcv::gam(stats::as.formula(eq), family = countreg::ztpoisson(), data = data.ab, method = method,knots=knots.ab[[i]])  #fitting the gam according to the specified equations

        mod.list.ab[[i]] <- gam.ab  #storing each of the fitted model in the empty list specified before
        edf.ab.list[[i]] <- summary(gam.ab)$s.table[,1]

        if(!is.null(gam.ab$smooth[[1]]$xp)){
          knots.smooth.ab<-list()
            for (j in 1:length(gam.ab$smooth)){
              knots.smooth.ab[[j]]<-gam.ab$smooth[[j]]$xp
            }
          smooth.names.ab<-c()
            for (j in 1:length(knots.smooth.ab)){
              smooth.names.ab[j]<-gam.ab$smooth[[j]]$term
            }
          names(knots.smooth.ab)<-smooth.names.ab
          knots.ab.list[[i]]<-knots.smooth.ab
        }

        mod.sel.ab[i, 2] <- round(sum(gam.ab$edf)-1,2)
        mod.sel.ab[i, 3] <- round(((gam.ab$null.deviance - gam.ab$deviance)/gam.ab$null.deviance) * 100, 2)  #calculating the % of deviance explained by each model
        mod.sel.ab[i, 4] <- round(stats::logLik(gam.ab)[1], 2)  #log-likelihood
        mod.sel.ab[i, 5] <- MuMIn::AICc(gam.ab)  #AICc
    }

    mod.sel.ab <- data.frame(mod.sel.ab, deltaAICc = mod.sel.ab$AICc - min(mod.sel.ab$AICc), w = round(MuMIn::Weights(mod.sel.ab$AICc),
        3))  #adding delta AICc and AIC weights to the model selection table
    exdev.ab <- mod.sel.ab$expl_Dev
    w.ab <- mod.sel.ab$w
    k.ab <- mod.sel.ab$k
    mod.sel.ab <- mod.sel.ab[order(mod.sel.ab$deltaAICc), ]  #order the model selection table according to increasing delta AICc

    best.ab <- list()  #empty list for the best models in the model selection table
    sub.ab <- subset(mod.sel.ab, w >= lim)  #taking the models which have a delta AIC less or equal to 4
    ID.ab <- as.integer(rownames(sub.ab))  #extracting the row IDs for each of the best models

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
        best.ab <- mod.list.ab[[ID.ab[1]]]
        fit.w.ab <- stats::predict(mod.list.ab[[ID.ab[1]]], newdata = data.ab, type = "response")  #taking fitted values of the model

        fit.w.ab.full <- stats::predict(mod.list.ab[[ID.ab[1]]], newdata = data, type = "response")  #taking predicted values (full dataset) of the model
        grid.w.ab <- stats::predict(mod.list.ab[[ID.ab[1]]], newdata = data.frame(grid,phat=1), type = "response")
    }  #end of conditional part

    eval.H <- ifelse(ab.full==0,fit.w.pa,fit.w.ab.full * fit.w.pa)  #multiplying pa with ab (full dataset) predictions (HURDLE)
    grid.H <- grid.w.ab * grid.w.pa
    res.H<-ab.full-eval.H

    if (SelectionTable) {
        cat("\n\n Model selection table for presence-absence submodel variants\n\n ")
        print(mod.sel.pa)
        cat("\n\n Selected model indices for presence-absence submodel\n\n ")
        print(ID.pa)
        cat("\n\n\n\n Model selection table for abundance submodel variants\n\n ")
        print(mod.sel.ab)
        cat("\n\n Selected model indices for abundance submodel\n\n ")
        print(ID.ab)
    }

    if (showSelectedModels) {
        cat("\n\n\n\n Selected presence-absence submodel variants\n\n ")
        print(best.pa)
        cat("\n\n\n\n Selected abundance submodel variants\n\n ")
        print(best.ab)  #print list with best models
    }

    return(list(models=list(pa=mod.list.pa,ab=mod.list.ab),
                info=list(ID.pa=ID.pa,k.pa=k.pa,weights.pa=w.pa,
                          edfs.pa=edf.pa.list,ID.ab=ID.ab,k.ab=k.ab,
                          weights.ab=w.ab,edfs.ab=edf.ab.list,
                          k.loc.pa=knots.pa.list,k.loc.ab=knots.ab.list,
                          exdev.pa=exdev.pa,exdev.ab=exdev.ab,method=method,lim=lim),
                grid_data=list(pa = grid.w.pa, ab = grid.w.ab, H = grid.H),
                fitted=list(pa=fit.w.pa,ab.full=fit.w.ab.full,ab=fit.w.ab),
                obs=ab.full,
                residuals=res.H))

}
