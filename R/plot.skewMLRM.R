plot.skewMLRM <-
function(x, ...)
{
  if(substring(x$"function",1,4)=="dist")
  {
 pc<-x$cut
 plot(x$Mahal,xlab='obs',ylab='Mahalanobis',ylim=c(0,1.1*max(c(x$Mahal,pc))), ...)
 graphics::abline(h=pc,lty=2,lwd=1,col='red')
  }
  else
  {
   cat("\nPlot only can be used with an skewMLRM object created with distMahal function", sep = "")
  }
}
