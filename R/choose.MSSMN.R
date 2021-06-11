choose.MSSMN <-
function(y, X=NULL, max.iter=1000, prec=1e-4, est.var=TRUE, criteria="AIC",cluster=FALSE)
{
  if(cluster!=TRUE) cluster=FALSE
  if(all(criteria!=c("AIC","BIC")))
        stop("criteria must be AIC or BIC")
  if (cluster==FALSE){      
  fit.msn=try(estimate.MSN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mstn=try(estimate.MSTN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mssl=try(estimate.MSSL(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mscn=try(estimate.MSCN(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl);j=1
      kname=c("MSN","MSTN","MSSL","MSCN")
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('SSMNfunctions.r')
     require(skewMLRM)
     XX=kname[j]
     estimate.dist.XX <- get(paste("estimate.", XX, sep = ""), mode = "function")
     fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
   }
   fit.msn=fit4[[1]]
   fit.mstn=fit4[[2]]
   fit.mssl=fit4[[3]]
   fit.mscn=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.msn)[1] & grepl("Error",fit.mstn)[1] & 
     grepl("Error",fit.mssl)[1] & grepl("Error",fit.mscn)[1])
  {stop("estimation problem in all the models in the SSMN class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.mstn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstn$AIC)}
  if(grepl("Error",fit.mssl)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssl$AIC)}
  if(grepl("Error",fit.mscn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscn$AIC)}
  if(grepl("Error",fit.msn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msn$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.mstn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstn$BIC)}
  if(grepl("Error",fit.mssl)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssl$BIC)}
  if(grepl("Error",fit.mscn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscn$BIC)}
  if(grepl("Error",fit.msn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msn$BIC)}
  }
  fitted<-c("STN","SSL","SCN","SN")[index]
  if(ncol(as.matrix(y))>1) fitted<-c("MSTN","MSSL","MSCN","MSN")[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.mstn, "2"=fit.mssl, "3"=fit.mscn, "4"=fit.msn)
  selected<- c("STN","SSL","SCN","SN")[index]
  if(ncol(as.matrix(y))>1) selected<-c("MSTN","MSSL","MSCN","MSN")[index]
  RVAL<- list(fitted.models=fitted, selected.model=selected, estimate=fit$estimate, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv)
  if(est.var)
  {
   se<-try(se.est(fit$estimate[,1],y,X,dist=selected),silent=TRUE)
   if(grepl("Error",se[1])) RVAL$warning="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
   else RVAL$estimate<-cbind(RVAL$estimate,se); colnames(RVAL$estimate)<-c("estimate","s.e.")
  } 
 RVAL
}
