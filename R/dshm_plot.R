#' Plotting prediction grid
#'
#' \code{dshm_plot} plots prediction grids, makes easier changing scale colours and intervals, provides the possibility to save the grid as a .shp file.
#'
#' @param prediction Predictions for each grid cell in the prediction grid.
#' @param grid Grid as SpatialPolygonsDataFrame.
#' @param probability If \code{TRUE} colour label is Pr(Presence), if \code{FALSE} is Abundance. Default is \code{FALSE}.
#' @param sightings sighting dataframe with at least three columns:
#' \itemize{
#'   \item x: sighting x coordinate.
#'   \item y: sighting y coordinate.
#'   \item size: sighting size.
#' }
#' @param plot_title Title of the plot.
#' @param scale_col Vector with colour scale values. Length should correspond to the length of 'scale_val'.
#' @param scale_lim Vector for lower and upper limits of the colour scale.
#' @param scale_val Vector for the scale values corresponding to colour transition as defined in 'scale_col'.
#' @param cex Size of the pixels in the map (visualization purpose only, for true scaled spatial map save as raster).
#' @param saveRaster If \code{TRUE} the map is saved as a raster file. Default is \code{TRUE}.
#' @param raster_name Name of the saved raster.
#' @details For more information about fitting Hurdle models, plotting model predictions and much more you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{fitting_Hurldle.pdf} tutorial.
#' @return Prediction map and raster image (if \code{saveRaster = TRUE}).
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}
#' @export
#'
dshm_plot <- function(prediction, grid, probability = FALSE, sightings = NULL, plot_title = NULL, scale_col = c("blue", "yellow",
    "red"), scale_lim = NULL, scale_val = NULL, cex = 1, saveRaster = FALSE, raster_name = NULL) {

    if (probability) {
        lab = "Pr(presence)"
    } else {
        lab = "Abundance"
    }

    if (saveRaster) {
        grid.H.reduced <- data.frame(id = grid$id, P = prediction)
        cat("\n\n Loading grid shapefile and generating raster...\n\n ")
        #grid.shp <- rgdal::readOGR(grid.shp, verbose = FALSE)
        grid.shp <- grid[, as.numeric(rownames(subset(data.frame(col.n=colnames(grid@data),colnames(grid@data)=="id"),col.n=="id")))]
        GRID.pred <- raster::merge(grid, grid.H.reduced, by = "id")

        r <- raster::raster()
        raster::crs(r) <- raster::crs(grid)
        raster::extent(r) <- raster::extent(grid)

        ras <- raster::rasterize(GRID.pred, r, "P")
        raster::writeRaster(ras, paste(raster_name), format = "GTiff", overwrite = TRUE)
    }


    grDevices::dev.new(width = 7, height = 6)
    if (!is.null(scale_lim)) {
        max <- max(scale_lim)
    } else {
        max <- max(prediction)
    }


    breaks<-c(0,0.25,0.5,0.75,1)

    plot_map <- ggplot2::ggplot() + ggplot2::geom_point(data = grid@data, ggplot2::aes(x = grid$x.coord, y = grid$y.coord, colour = prediction/max), cex = cex,
        pch = 15) + ggplot2::scale_colour_gradientn(colours = scale_col, limits = if (!is.null(scale_lim)) {
        limits = scale_lim/max
    }, values = if (!is.null(scale_val)) {
        values = scale_val/max
    }, breaks=breaks, labels = round(max*breaks,2)) + ggplot2::labs(title = paste(plot_title), x = "x coord", y= "y coord", colour = lab)

    if (!is.null(sightings)) {
        plot_map + ggplot2::geom_point(data = sightings, ggplot2::aes(x = sightings$x, y = sightings$y, size = sightings$size)) + ggplot2::scale_size_continuous(range = c(1,
            5)) + ggplot2::labs(size = "Group size")
    } else {
        plot_map
    }
}
