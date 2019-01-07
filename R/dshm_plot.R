#' Plots prediction grid and saves it as a raster image
#'
#' @param prediction Prediction grid returned by dshm_fit$grid_data.
#' @param species Character. Species name (used for plotting and saving the raster).
#' @param grid Grid shapefile (SpatialPolygonsDataFrame object) with specific projection.
#' @param map_type Character. Type of map to plot: 'presence_absence', 'abundance', 'Hurdle', or 'SD'.
#' @param sightings Dataframe with at least three columns: 'x' and 'y' coordinates and 'size' for loaction and size of each sighitng, respectively.
#' @param scale_col Vector with colour scale values. Length should correspond to the length of 'scale_val'.
#' @param scale_lim Vector for lower and upper limits of the colour scale.
#' @param scale_val Vector for the scale values corresponding to colour transition as defined in 'scale_col'.
#' @param cex Size of the pixels in the map (visualization purpose only, for true scaled spatial map save as raster).
#' @param saveRaster If TRUE the map is saved as a raster file. A grid shapefile with projection is needed. Default is FALSE.
#' @return Map and raster. Raster name contains information in 'map_type' and 'species'.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}

#' @export
dshm_plot <- function(prediction, species = NULL, grid, map_type = "Hurdle", sightings = NULL,scale_col = c("blue", "yellow",
    "red"), scale_lim = NULL, scale_val = NULL, cex = 1, saveRaster = FALSE) {

    if (map_type == "presence_absence") {
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
        raster::writeRaster(ras, paste(map_type,species), format = "GTiff", overwrite = TRUE)
    }


    dev.new(width = 7, height = 6)
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
    }, breaks=breaks, labels = round(max*breaks,2)) + ggplot2::labs(title = paste(map_type,species), x = "x coord", y= "y coord", colour = lab)

    if (!is.null(sightings)) {
        plot_map + ggplot2::geom_point(data = sightings, ggplot2::aes(x = sightings$x, y = sightings$y, size = sightings$size)) + ggplot2::scale_size_continuous(range = c(1,
            5)) + ggplot2::labs(size = "Group size")
    } else {
        plot_map
    }
}
