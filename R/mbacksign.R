mbacksign <-
function(y, X=NULL, max.iter=1000,prec=1e-4, dist="MN", significance=0.05) 
{
  if (!any(dist == c("MN","MT","MSL","MCN","MSN","MSNC","MSTEC","MSTN","MSTT",
"MSSLEC","MSSL","MSSL2","MSCN","MSCN2","MSCEC")))
     stop("distribution is not recognized")
  if(!is.numeric(significance))
        stop("significance must be numeric")  
  if(significance<=0)
        stop("significance must be a positive value")  
  if(significance>0.10)
        stop("significance is usually less than or equal to 0.10")  
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
  tt<-table(unlist(sapply(X,pos.ones)))
  ind.interc<-as.numeric(names(tt)[which(tt==nrow(y))])
  m=ncol(X[[1]])
  z.critico<-qnorm(1-significance/2)
  estimate.dist <- get(paste("estimate.", dist, sep = ""), mode = "function")
  object.1<-estimate.dist(y=y,X=X,max.iter=max.iter,prec=prec, est.var=FALSE)
  beta.est <- object.1$"estimate"[c(1:m)[-ind.interc],1]
  if(dist!="MSTEC" & dist!="MSSLEC" & dist!="MSCEC") se<-try(se.est(object.1$"estimate"[,1],y,X,dist=dist),silent=TRUE)
  if(dist=="MSTEC" | dist=="MSSLEC") se<-try(se.est(c(object.1$"estimate"[,1],object.1$"nu"),y,X,dist=dist),silent=TRUE)
  if(dist=="MSCEC") se<-try(se.est(c(object.1$"estimate"[,1],object.1$"nu",object.1$"gamma"),y,X,dist=dist),silent=TRUE)
  if(!is.numeric(se)) stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
  ep.beta <- se[c(1:m)[-ind.interc]]
  xmenos.multi <- X
  dim.beta <- length(beta.est)
  test.t <- abs(beta.est/ep.beta)
  #for(i in 1:dim.beta)
  #{
    #if(test.t[i] == Inf) test.t[i] <- NA
    #if(is.nan(test.t[i])) test.t[i] <- NA
  #}
  #aux <- X[[1]]
  #aux2 <- ncol(aux)/nrow(aux)
  #for(i in 1:ncol(y)) {
  #  test.t[aux2*(i-1)+1] <- NA
  #}
  #print(test.t)
      lista.final <- list(comment="The final model considered all the betas", estimate=object.1$"estimate",
logLik=object.1$"logLik", AIC=object.1$"AIC",BIC=object.1$"BIC", conv=object.1$conv)
  criterio <- sum(test.t - z.critico > 0, na.rm = TRUE)
  if(criterio > 0){
    t.min <- c(1:m)[-ind.interc][which.min(test.t)]
    historico <- c()
    cont <- sum(is.na(test.t))
   while(criterio < (length(test.t)- cont) && (dim.beta > 1)) {
      #print(t.min)
      #print(criterio)
      historico <- c(historico,rownames(object.1$"estimate")[1:m][t.min])
names.betas<-rownames(object.1$"estimate")[1:m][-(t.min)]
      xmenos.multi <- lapply(xmenos.multi, function(x) { x[,-(t.min)]})
      m=ncol(xmenos.multi[[1]])
      tt<-table(unlist(sapply(xmenos.multi,pos.ones)))
      ind.interc<-as.numeric(names(tt)[which(tt==nrow(y))])
      object.1 <- estimate.dist(y=y,X=xmenos.multi,max.iter=max.iter,prec=prec,est.var=FALSE)
rownames(object.1$"estimate")[1:m]<-names.betas
      beta.est<- object.1$"estimate"[c(1:m)[-ind.interc],1]
  if(dist!="MSTEC" & dist!="MSSLEC" & dist!="MSCEC") se<-try(se.est(object.1$"estimate"[,1],y,X=xmenos.multi,dist=dist),silent=TRUE)
  if(dist=="MSTEC" | dist=="MSSLEC") se<-try(se.est(c(object.1$"estimate"[,1],object.1$"nu"),y,X=xmenos.multi,dist=dist),silent=TRUE)
  if(dist=="MSCEC") se<-try(se.est(c(object.1$"estimate"[,1],object.1$"nu",object.1$"gamma"),y,X=xmenos.multi,dist=dist),silent=TRUE)
      if(!is.numeric(se)) stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
      ep.beta <- se[c(1:m)[-ind.interc]]
      dim.beta <- length(beta.est)
      test.t <- abs(beta.est/ep.beta)
      #aux <- xmenos.multi[[1]]
      #aux2 <- ncol(aux)/nrow(aux)
      #for(i in 1:ncol(y)) {
      #  test.t[aux2*(i-1)+1] <- NA
      #}
      #for(i in 1:(dim.beta-1))
      #{
        #if(test.t[i] == Inf) test.t[i] <- NA
        #if(is.nan(test.t[i])) test.t[i] <- NA
      #}
      #print(test.t)
      #cont <- sum(is.na(test.t))
      t.min <- c(1:m)[-ind.interc][which.min(test.t)]
      criterio <- sum(test.t - z.critico > 0, na.rm = TRUE)
      lista.final <- list(comment=paste("The final model eliminated",length(historico),"betas"),
                   estimate=object.1$"estimate")
if(dist=="MSTEC" | dist=="MSSLEC") lista.final$nu=object.1$nu
if(dist=="MSCEC") lista.final$nu=object.1$nu; lista.final$gamma=object.1$gamma
lista.final$logLik=object.1$"logLik" 
lista.final$AIC=object.1$"AIC"
lista.final$BIC=object.1$"BIC" 
lista.final$conv=object.1$"conv"
lista.final$"eliminated"=historico
    } } 
 if(dist!="MSTEC" & dist!="MSSLEC" & dist!="MSCEC") se<-try(se.est(object.1$"estimate"[,1],y,X=xmenos.multi,dist=dist),silent=TRUE)
 if(dist=="MSTEC" | dist=="MSSLEC") se<-try(se.est(object.1$"estimate"[,1],y,X=xmenos.multi,dist=dist, nu=object.1$"nu"),silent=TRUE)
 if(dist=="MSCEC") se<-try(se.est(object.1$"estimate"[,1],y,X=xmenos.multi,dist=dist,nu=object.1$"nu",gamma=object.1$"gamma"),silent=TRUE)
  if(is.numeric(se))
  {
   lista.final$"estimate"<-cbind(lista.final$"estimate"[,1],se)
   colnames(lista.final$estimate)<-c("estimate","s.e.")
  }
 return(lista.final)
}
