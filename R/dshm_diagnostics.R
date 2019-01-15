#' Dignostics tool for Hurdle models
#'
#' \code{dshm_diagnostics} computes the Hurdle model cumulative distribution function (CDF) and it plots it against the empirical distribution function (EDF). It also calculates Kolmogorov-Smirnov test statistics.
#'
#' @param model Hurdle model fitted through \code{dshm_fit}.
#' @param mute If \code{TRUE} returns p-value and test statistics for Kolmogorov-Smirnov test. Default is \code{FALSE}.
#' @param plot If \code{TRUE} prints CDF vs EDF and fitted values vs. observed values plots. Default is \code{TRUE}.
#' @param plot.n The number of available plots. If \code{plot.n = 1} then the function prints only the CDF vs EDF plot while if \code{plot.n = 2} the function prints both plots for CDF vs EDF plot and fitted values vs. observed values. Default is \code{plot.n = 1}.
#' @return Two plots for CDF vs. EDF and observed values vs. fitted values. Kolmogorov-Smirnov test statistics and p-value.
#' @details The Hurdle model CDF is calculated using the following equation:
#' \deqn{CDF(n) = (1 - p)(1 - \lambda) + p\lambdaCDF(P(n > 0))}
#' Where \eqn{n} are the number of observations, \eqn{p} is the probability of presence, \eqn{\lambda} is 0 for absence and 1 for presence, and \eqn{CDF(P(n > 0))} is the zero-trucated Poisson cumulative distribution function, i.e the probability of observing x given the zero-trucated Poisson parameter.
#'
#' For more information about fitting Hurdle models you can download the \href{http://github.com/FilippoFranchini/dshm/blob/master/vignettes}{fitting_Hurdle.pdf} tutorial.
#' @author Filippo Franchini \email{filippo.franchini@@outlook.com}
#' @export
#'
dshm_diagnostics<-function(model,mute=FALSE,plot=TRUE,plot.n=1){

  obs <- model$obs
  fit.pa <- model$fitted$pa
  fit.ab <- model$fitted$ab.full

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
      graphics::plot(CDF,EDF,pch=19,col=grDevices::rgb(0,0,0,0.5),cex=1,main="Q-Q Plot")
      graphics::lines(c(-1,2),c(-1,2))
      max.val<-max(fit.pa*fit.ab,obs)
      graphics::plot(fit.pa*fit.ab,obs,ylab="Response",xlab="Fitted Values",main="Response vs. Fitted Values",pch=19,col=grDevices::rgb(0,0,0,0.5),ylim=c(0,max.val),xlim=c(0,max.val))
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
