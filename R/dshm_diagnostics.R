#' Computes CDF and plots CDF vs. EDF relationship for dshm objects
#'
#' @param fit.pa Fitted values from binomial submodel
#' @param fit.ab Fitted values from zero-truncated Poisson submodel
#' @param obs Observation (same length as fit.pa and fit.ab)
#' @param mute If TRUE returns p-value and test statistics for Kolmogorov-Smirnov test. Default is FALSE.
#' @param plot If TRUE plots CDF vs EDF and fitted vs. observed values plots. Default is TRUE.
#' @param plot.n 1 for CDF vs EDF plot and 2 for both CDF vs EDF and fitted vs. observed values plots.

#' @export
dshm_diagnostics<-function(fit.pa,fit.ab,obs,mute=FALSE,plot=TRUE,plot.n=1){

  gof<-data.frame(obs) #create a dataset with all observations
  gof$pa<-rep(0,length(gof[,1]))
  for (i in 1:length(gof[,1])){
    gof$pa[i]<-ifelse(gof$obs[i]==0,0,1)
  }
  gof$pa_hat<-fit.pa
  gof$ab_hat<-fit.ab
  gof$CDF<-rep(0,length(gof[,1])) #adding an empty column for CDF values

  #Hurdle model CDF(n) = (1-p)(1-lambda)+pCDF(P(n>0))lambda, where lambda is 0 or 1 for absence and presence respectively
  for (i in 1:length(gof[,1])){ #calculating the CDF according to presence and absence
    if(gof$obs[i]==0){ #if absence, then CDF is the probability of absence from the binomial model
      gof$CDF[i]<-1-gof$pa_hat[i]
    } else { #if presence then CDF is probability of presence multiplied by the CDF of zero-truncated poisson of fitted value
      gof$CDF[i]<-gof$pa_hat[i]*countreg::pztpois(q = gof$obs[i] ,lambda = gof$ab_hat[i])
    }
  }

  noise <- stats::rnorm(length(gof$CDF),mean=0,sd=0.00001)
  CDF<-gof$CDF+noise
  CDF<-stats::punif(CDF)
  EDF<-rank(CDF)/length(gof$CDF)
  ord<-order(CDF)
  CDF<-CDF[ord]
  EDF<-EDF[ord]

  if(plot){
    if(plot.n==1){
      plot(CDF,EDF,pch=19,col=grDevices::rgb(0,0,0,0.5),cex=1,main="")
      graphics::lines(c(-1,2),c(-1,2))
    } else if (plot.n==2){
      graphics::par(mfrow=c(1,2))
      plot(CDF,EDF,pch=19,col=grDevices::rgb(0,0,0,0.5),cex=1,main="Q-Q Plot")
      lines(c(-1,2),c(-1,2))
      max.val<-max(fit.pa*fit.ab,obs)
      plot(fit.pa*fit.ab,obs,ylab="Response",xlab="Fitted Values",main="Response vs. Fitted Values",pch=19,col=rgb(0,0,0,0.5),ylim=c(0,max.val),xlim=c(0,max.val))
      graphics::lines(c(-100,100),c(-100,100))
    }
  }
  ks<-stats::ks.test(CDF,"punif",0,1,alternative="two.sided")

  if(!mute){
    print(ks)
  } else {
    return(list(D=ks$statistic[[1]],p.value=ks$p.value))
  }
}
