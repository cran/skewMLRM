choose.models <-
function(y, X=NULL, max.iter=1000, prec=1e-4, est.var=TRUE, criteria="AIC",cluster=FALSE)
{
  if(cluster!=TRUE) cluster=FALSE
  if(all(criteria!=c("AIC","BIC")))
        stop("criteria must be AIC or BIC")
  if (cluster==FALSE){      
  fit.msmn=try(choose.MSMN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mssmn=try(choose.MSSMN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.msmsn=try(choose.MSMSN(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.msmsnc=try(choose.MSMSNC(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl);j=1
      kname=c("MSMN","MSSMN","MSMSN","MSMSNC")
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('choosefunctions.r')
     require(skewMLRM)
     XX=kname[j]
     choose.XX <- get(paste("choose.", XX, sep = ""), mode = "function")
     fit4= try(choose.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
   }
   fit.msmn=fit4[[1]]
   fit.mssmn=fit4[[2]]
   fit.msmsn=fit4[[3]]
   fit.msmsnc=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.msmn)[1] & grepl("Error",fit.mssmn)[1] & 
     grepl("Error",fit.msmsn)[1] & grepl("Error",fit.msmsnc)[1])
  {stop("estimation problem in all the models in the SMN class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.msmn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.msmn$AIC)}
  if(grepl("Error",fit.mssmn)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssmn$AIC)}
  if(grepl("Error",fit.msmsn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.msmsn$AIC)}
  if(grepl("Error",fit.msmsnc)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msmsnc$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.msmn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.msmn$BIC)}
  if(grepl("Error",fit.mssmn)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssmn$BIC)}
  if(grepl("Error",fit.msmsn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.msmsn$BIC)}
  if(grepl("Error",fit.msmsnc)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msmsnc$BIC)}
  }
  fitted<-c(fit.msmn$fitted.models, fit.mssmn$fitted.models,
fit.msmsn$fitted.models, fit.msmsnc$fitted.models)
  sel<-c(fit.msmn$selected.model, fit.mssmn$selected.model,
fit.msmsn$selected.model, fit.msmsnc$selected.model)[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.msmn, "2"=fit.mssmn, "3"=fit.msmsn, "4"=fit.msmsnc)
  selected<- sel[index]
 if (index == 4) {
     if(selected=="MSNC") RVAL <- list(fitted.models = fitted, selected.model = selected, 
            estimate = fit$estimate, logLik = fit$logLik, 
            AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, conv = fit$conv)
     if(selected=="MSTEC" | selected=="MSSLEC") RVAL <- list(fitted.models = fitted, selected.model = selected, 
            estimate = fit$estimate, nu = fit$nu, logLik = fit$logLik, 
            AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
            conv = fit$conv)
     if(selected=="MSCEC") RVAL <- list(fitted.models = fitted, selected.model = selected, 
            estimate = fit$estimate, nu = fit$nu, gamma=fit$gamma, logLik = fit$logLik, 
            AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, conv = fit$conv)
    }
    else {
        RVAL <- list(fitted.models = fitted, selected.model = selected, 
            estimate = fit$estimate, logLik = fit$logLik, AIC = fit$AIC, 
            BIC = fit$BIC, iterations = fit$iterations, 
            conv = fit$conv)
    }
    if (est.var) {
        if(!any(selected == c("MSTEC", "MSSLEC", "MSCEC"))) se <- try(se.est(fit$estimate[, 1], y, X, dist = selected), silent = TRUE)
        if(any(selected==c("MSTEC", "MSSLEC"))) se <- try(se.est(fit$estimate[, 1], y, X, dist = selected, nu=fit$nu), silent = TRUE)
        if(selected=="MSCEC") se <- try(se.est(fit$estimate[, 1], y, X, dist = selected, nu=fit$nu, gamma=fit$gamma), silent = TRUE)
        if (grepl("Error", se[1])) 
            RVAL$warning = "Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
        else RVAL$estimate <- cbind(RVAL$estimate, se);colnames(RVAL$estimate) <- c("estimate", "s.e.")
    }
 RVAL
}
