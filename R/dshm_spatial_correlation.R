#' Calculates spatial correlation for dshm
#'
#' @param coord Centroids of each segment.
#' @param z Fitted values of model residuals.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.


dshm_spatial_correlation<-function(coord,z,xlab,ylab){
  data<-data.frame(x=coord[,1],y=coord[,2],z=z)
  surf<-spatial::surf.ls(6,na.omit(data))
  correl<-spatial::correlogram(krig = surf,nint = length(data[,1]),plotit=FALSE)[]
  plot(correl$x/1000,correl$y,ylim=c(-1,1),type="b",xlab=xlab,ylab=ylab,pch=19,col=rgb(0,0,0,0.5))
  abline(h=0)

  info<-length(subset(correl$y,correl$y>=-0.1 & correl$y<=0.1))/length(correl$y)*100 #proportion of correlation within -0.1 and 0.1
  cat(paste(round(info,1),"% of the data ranges between -0.1 and 0.1"))
}
