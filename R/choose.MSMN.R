choose.MSMN <-
function(y, X=NULL, max.iter=1000, prec=1e-4, est.var=TRUE, criteria="AIC",cluster=FALSE)
{
  if(cluster!=TRUE) cluster=FALSE
  if(all(criteria!=c("AIC","BIC")))
        stop("criteria must be AIC or BIC")
  if (cluster==FALSE){      
  fit.mn=try(estimate.MN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mtn=try(estimate.MT(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.msl=try(estimate.MSL(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mcn=try(estimate.MCN(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl)
      kname=c("MN","MT","MSL","MCN");j=1
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('SMNfunctions.r')
     require(skewMLRM)
     XX=kname[j]
     estimate.dist.XX <- get(paste("estimate.", XX, sep = ""), mode = "function")
     fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
   }
   fit.mn=fit4[[1]]
   fit.mtn=fit4[[2]]
   fit.msl=fit4[[3]]
   fit.mcn=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.mn)[1] & grepl("Error",fit.mtn)[1] & 
     grepl("Error",fit.msl)[1] & grepl("Error",fit.mcn)[1])
  {stop("estimation problem in all the models in the SMN class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.mtn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mtn$AIC)}
  if(grepl("Error",fit.msl)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.msl$AIC)}
  if(grepl("Error",fit.mcn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mcn$AIC)}
  if(grepl("Error",fit.mn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.mn$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.mtn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mtn$BIC)}
  if(grepl("Error",fit.msl)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.msl$BIC)}
  if(grepl("Error",fit.mcn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mcn$BIC)}
  if(grepl("Error",fit.mn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.mn$BIC)}
  }
  fitted<-c("T","SL","CN","N")[index]
  if(ncol(as.matrix(y))>1) fitted<-c("MT","MSL","MCN","MN")[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.mtn, "2"=fit.msl, "3"=fit.mcn, "4"=fit.mn)
  selected<- c("T","SL","CN","N")[index]
  if(ncol(as.matrix(y))>1) selected<-c("MT","MSL","MCN","MN")[index]
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
