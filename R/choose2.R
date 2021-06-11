choose2 <-
function(y, X=NULL, max.iter=1000, prec=1e-4, class="MSMN", est.var=TRUE, criteria="AIC",criteria.cov="significance", significance=0.05, cluster=FALSE)
{
 if(cluster!=TRUE) cluster=FALSE
 if (!any(class == c("MSMN","MSSMN","MSMSN","MSMSNC","ALL")))
     stop("class is not recognized")
  if (!any(criteria == c("AIC","BIC")))
     stop("criteria must be AIC or BIC")
  if (!any(criteria.cov == c("AIC","BIC","significance")))
     stop("criteria for covariates must be AIC, BIC or significance")
  if(criteria.cov=="significance")
  {
   if(!is.numeric(significance))
         stop("significance must be numeric")  
   if(significance<=0)
         stop("significance must be a positive value")  
   if(significance>0.10)
         stop("significance is usually less than or equal to 0.10")  
  }
sum.ones<-function(x) sum(x==1)
pos.ones<-function(x) which(apply(x,2,sum)==1)
 y<-as.matrix(y)
if(!is.matrix(y))
        stop("y must have at least one element")
  if(is.null(X)){X<-array(c(diag(ncol(y))),c(ncol(y),ncol(y),nrow(y)))}
  if(is.array(X)==FALSE & is.list(X)==FALSE)
        stop("X must be an array or a list")
  if(is.array(X))
  {Xs<-list()
     for (i in 1:nrow(y)){
    Xs[[i]]<- matrix(t(X[,,i]),nrow=ncol(y))};X<-Xs} 
  m=sum(mapply(sum.ones,X))/(length(X)*nrow(X[[1]]))
  q=ncol(y)
  if((ncol(X[[1]])-m)==q)
       stop("X have only the intercept term(s)")
  if (ncol(y) != nrow(X[[1]]))
        stop("y does not have the same number of columns than X")
  if (nrow(y) != length(X))
        stop("y does not have the same number of observations than X")
 if(class=="MSMN") cluster=FALSE
 if(any(class == c("MSMN","MSSMN","MSMSN","MSMSNC")))
  {
    choose.XX <- get(paste("choose.", class, sep = ""), mode = "function")
    fit=try(choose.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE, cluster=cluster),silent=TRUE)
  }
  if(class=="ALL")
  {
   fit=try(choose.models(y, X, max.iter=max.iter, prec=prec, est.var=FALSE, criteria=criteria, cluster=cluster),silent=TRUE)
  }
  if(grepl("Error",fit)[1])
  {stop(paste("estimation problem in all the models in the ",class,"class"))}
  fitted<-fit$fitted
  selected<-fit$selected
  dist<-selected
  lista.final<- list(fitted.models=fitted, selected.model=selected)
  if(criteria.cov=="significance")
  {
     object<-mbacksign(y, X=X, max.iter=max.iter, prec=prec, dist=dist, significance=significance) 
  }
  else
  {
    object=mbackcrit(y, X=X, max.iter=max.iter, prec=prec, dist=dist, criteria=criteria, cluster=cluster)
  }
   if(selected=="MSTEC" | selected=="MSSLEC") lista.final$nu=object$"nu"
   if(selected=="MSCEC") lista.final$nu=object$"nu"; lista.final$gamma=object$"gamma"
   lista.final$comment=object$"comment"
   lista.final$estimate=object$"estimate"
   lista.final$logLik=object$"logLik"
   lista.final$AIC=object$"AIC"
   lista.final$BIC=object$"BIC"
   lista.final$conv=object$"conv"
   if(!is.null(object$"eliminated"))
   {
     lista.final$"eliminated"=object$"eliminated"
   }
   if(!is.null(object$"warnings"))
   {
     lista.final$warnings=object$"warnings"
   }
 lista.final  
}
