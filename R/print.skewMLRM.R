print.skewMLRM <-
function(x, digits = max(3L, getOption("digits") - 3L), ...) 
{
cat("\nUsed function: ",x$"function","\n", sep = "")
  if(substring(x$"function",1,7)=="choose."){
    pp<-x$fitted.models[1]
    if(length(x$fitted.models)>1)
    {
     for(j in 2:length(x$fitted.models))
     {pp<-paste(pp,"/",x$fitted.models[j])}
    }
    cat("Fitted models: ",pp,"\nSelected model: ",x$selected,"\nClass: ", x$class,"\n\n", sep = "")
  }
  if(substring(x$"function",1,5)=="mback"){
   cr<-ifelse(x$choose.crit=="sign",paste("significance (",100*x$significance,"%)",sep=""),x$choose.crit)
    cat("Fitted model: ",x$dist,"\nClass: ", x$class,"\nMethod to select covariates: ",cr,"\n\n", sep = "")
  }
  if(substring(x$"function",1,8)=="estimate") cat("Fitted model: ",x$dist,"\nClass: ", x$class,"\n\n", sep = "")
  if(substring(x$"function",1,7)=="choose2"){
    pp<-x$fitted.models[1]
    if(length(x$fitted.models)>1)
    {
     for(j in 2:length(x$fitted.models))
     {pp<-paste(pp,"/",x$fitted.models[j])}
    }
    cat("Fitted models: ",pp,"\nSelected model: ",x$selected,"\nClass: ", x$class,"\n\n", sep = "")
  }
  if(substring(x$"function",1,4)=="dist"){
        cat("Compute Mahalanobis distance\n")
    cat("Fitted model: ",x$dist,"\nClass: ", x$class,"\n", sep = "")
}
    if (length(x$coefficients)) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2L, 
            quote = FALSE)
    }
    cat("\n")
    invisible(x)
}
