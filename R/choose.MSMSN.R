choose.MSMSN <-
function(y, X=NULL, max.iter=1000, prec=1e-4, est.var=TRUE, criteria="AIC",cluster=FALSE)
{
  if(cluster!=TRUE) cluster=FALSE
  if(all(criteria!=c("AIC","BIC")))
        stop("criteria must be AIC or BIC")
  if (cluster==FALSE){      
  fit.msn=try(estimate.MSN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mstt=try(estimate.MSTT(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mssl2=try(estimate.MSSL2(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mscn2=try(estimate.MSCN2(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl);j=1
      kname=c("MSN","MSTT","MSSL2","MSCN2")
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('SMSNfunctions.r')
     require(skewMLRM)
     XX=kname[j]
     estimate.dist.XX <- get(paste("estimate.", XX, sep = ""), mode = "function")
     fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
   }
   fit.msn=fit4[[1]]
   fit.mstt=fit4[[2]]
   fit.mssl2=fit4[[3]]
   fit.mscn2=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.mstt)[1] & grepl("Error",fit.mssl2)[1] & 
grepl("Error",fit.mscn2)[1] & grepl("Error",fit.msn)[1])
  {stop("estimation problem in all the models in the SMSN class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.mstt)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstt$AIC)}
  if(grepl("Error",fit.mssl2)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssl2$AIC)}
  if(grepl("Error",fit.mscn2)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscn2$AIC)}
  if(grepl("Error",fit.msn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msn$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.mstt)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstt$BIC)}
  if(grepl("Error",fit.mssl2)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssl2$BIC)}
  if(grepl("Error",fit.mscn2)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscn2$BIC)}
  if(grepl("Error",fit.msn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msn$BIC)}
  }
  fitted<-c("STT","SSL2","SCN2","SN")[index]
  if(ncol(as.matrix(y))>1) fitted<-c("MSTT","MSSL2","MSCN2","MSN")[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.mstt, "2"=fit.mssl2, "3"=fit.mscn2, "4"=fit.msn)
  selected<- c("STT","SSL2","SCN2","SN")[index]
  if(ncol(as.matrix(y))>1) selected<-c("MSTT","MSSL2","MSCN2","MSN")[index]
  RVAL<- list(fitted.models=fitted, selected.model=selected, estimate=fit$estimate, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv)
  RVAL
}
