mbackcrit <-
function(y, X=NULL, max.iter=1000, prec=1e-4, dist="MN", criteria="AIC", cluster=FALSE)
{
  if (!any(dist == c("MN","MT","MSL","MCN","MSN","MSNC","MSTEC","MSTN","MSTT",
"MSSLEC","MSSL","MSSL2","MSCN","MSCN2","MSCEC")))
     stop("distribution is not recognized")
  if(cluster!=TRUE) cluster=FALSE
  if (!(criteria == "AIC" | criteria == "BIC"))
        stop("AIC or BIC criterion must be chosen")
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
  m=sum(mapply(sum.ones,X))/length(X)
  q=ncol(y)
  if((ncol(X[[1]])-m)==q)
       stop("X have only the intercept term(s)")
  if (ncol(y) != nrow(X[[1]]))
        stop("y does not have the same number of columns than X")
  if (nrow(y) != length(X))
        stop("y does not have the same number of observations than X")
  estimate.dist <- get(paste("estimate.", dist, sep = ""), mode = "function")
  X.complete <- X
  fit.complete<-estimate.dist(y=y,X=X.complete,max.iter=max.iter,prec=prec,est.var=FALSE)
  names.betas<-rownames(fit.complete$"estimate")[1:m]
  tt<-table(unlist(sapply(X,pos.ones)))
  ind.interc<-as.numeric(names(tt)[which(tt==nrow(y))])
  m=ncol(X[[1]])
  if(criteria=="AIC") criteria.complete<-fit.complete$AIC
  if(criteria=="BIC") criteria.complete<-fit.complete$BIC
       lista.final <- list(comment="The final model considered all the betas", estimate=fit.complete$"estimate",
logLik=fit.complete$"logLik", AIC=fit.complete$"AIC",BIC=fit.complete$"BIC", conv=fit.complete$conv)
   historico <- c()
    flag <- FALSE
   while(!flag) {
   if (cluster==FALSE | any(dist == c("MN","MT","MSL","MCN","MSN","MSTN","MSTT",
"MSSL","MSSL2","MSCN","MSCN2"))){
   criteria.complete<-c()
   for(j in 1:length(c(1:m)[-ind.interc]))
   {
    #X.complete.j<-array(c(diag(ncol(y))),c(ncol(y),ncol(y),nrow(y)))
    #if(j>2) X.complete.j <- lapply(X.complete, function(x) { x[,-j]})
    X.complete.j <- lapply(X.complete, function(x) { x[,-(c(1:m)[-ind.interc][j])]})
     fit.complete.j<-estimate.dist(y=y,X=X.complete.j,max.iter=max.iter,prec=prec,est.var=FALSE)
     if(criteria=="AIC") criteria.complete[j]<-fit.complete.j$AIC
     if(criteria=="BIC") criteria.complete[j]<-fit.complete.j$BIC
   } }
   # paralelo
   #if (cluster!=FALSE & any(dist==c("MSNC","MSTEC","MSSLEC","MSCEC"))){ 
   else{
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl)
    if(criteria=="AIC")
    {
j=1
    criteria.complete <- foreach (j = 1:length(c(1:m)[-ind.interc]), .combine=rbind) %dopar%{
     #source('covariatesfunctions.r')
     require(skewMLRM)
     X.complete.j <- lapply(X.complete, function(x) { x[,-(c(1:m)[-ind.interc][j])]})
     #print(dim(X.complete.j))
     fit.complete.j<-estimate.dist(y=y,X=X.complete.j,max.iter=max.iter,prec=prec,est.var=FALSE)
     criteria.complete<-fit.complete.j$AIC
     }
     criteria.complete<-c(criteria.complete)
    }
    if(criteria=="BIC")
    {
    j=1
    criteria.complete <- foreach (j = 1:length(c(1:m)[-ind.interc]), .combine=rbind) %dopar%{
     #source('covariatesfunctions.r')
     require(skewMLRM)
     X.complete.j <- lapply(X.complete, function(x) { x[,-(c(1:m)[-ind.interc][j])]})
     #print(dim(X.complete.j))
     fit.complete.j<-estimate.dist(y=y,X=X.complete.j,max.iter=max.iter,prec=prec,est.var=FALSE)
     criteria.complete<-fit.complete.j$BIC
     }
     criteria.complete<-c(criteria.complete)
    }
   parallel::stopCluster(cl)
   }
   #
   ind=c(1:m)[-ind.interc][which.min(criteria.complete)]
   crit.complete<-ifelse(criteria=="AIC",fit.complete$AIC,fit.complete$BIC)
   if(min(criteria.complete)>=crit.complete) {flag<-TRUE}
   if(min(criteria.complete)<crit.complete) {X.complete <- lapply(X.complete, function(x) { x[,-ind]});
      historico <- c(historico,rownames(fit.complete$"estimate")[1:m][ind])
names.betas<-rownames(fit.complete$"estimate")[1:m][-ind]
  tt<-table(unlist(sapply(X,pos.ones)))
  ind.interc<-as.numeric(names(tt)[which(tt==nrow(y))])
   m=ncol(X.complete[[1]])
   fit.complete<-estimate.dist(y=y,X=X.complete,max.iter=max.iter,prec=prec,est.var=FALSE)
   rownames(fit.complete$"estimate")[1:m]<-names.betas
   lista.final <- list(comment=paste("The final model eliminated",length(historico),"betas"),
         estimate=fit.complete$"estimate")
if(dist=="MSTEC" | dist=="MSSLEC") lista.final$nu=fit.complete$nu
if(dist=="MSCEC") lista.final$nu=fit.complete$nu; lista.final$gamma=fit.complete$gamma
lista.final$logLik=fit.complete$"logLik" 
lista.final$AIC=fit.complete$"AIC"
lista.final$BIC=fit.complete$"BIC" 
lista.final$conv=fit.complete$"conv"
lista.final$"eliminated"=historico
  rownames(lista.final$"estimate")[1:m]<-names.betas
 }
}
  lista.final$estimate<-fit.complete$estimate
 if(dist!="MSTEC" & dist!="MSSLEC" & dist!="MSCEC") se<-try(se.est(fit.complete$"estimate"[,1],y,X=X.complete,dist=dist),silent=TRUE)
 if(dist=="MSTEC" | dist=="MSSLEC") se<-try(se.est(fit.complete$"estimate"[,1],y,X=X.complete,dist=dist, nu=fit.complete$"nu"),silent=TRUE)
 if(dist=="MSCEC") se<-try(se.est(fit.complete$"estimate"[,1],y,X=X.complete,dist=dist,nu=fit.complete$"nu",gamma=fit.complete$"gamma"),silent=TRUE)
  if(is.numeric(se))
  {
   lista.final$"estimate"<-cbind(lista.final$"estimate"[,1],se)
   colnames(lista.final$estimate)<-c("estimate","s.e.")
  }
  rownames(lista.final$"estimate")[1:m]<-names.betas
  return(lista.final)
}
