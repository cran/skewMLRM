summary.skewMLRM <-
function(object, ...)
{
if(object$"function"=="distMahal")
{
name.dist<-switch(object$dist, "MN"="Normal", "MT"="t-Student", "MSL"="Slash", "MCN"="Contaminated Normal",
 "MSN"="Skew-Normal","MSTN"="Skew-t-Normal", "MSSL"="Skew-Slash-Normal", "MSCN"="Skew-Contaminated-Normal",
 "MSTT"="Skew-t", "MSSL2"="Skew-Slash", "MSCN2"="Skew-Contaminated-Normal",
 "MSNC"="Skew-Normal-Cauchy","MSTEC"="Skew-t-Expectec-Cauchy", "MSSLEC"="Skew-slash-Expected-Cauchy", "MSCEC"="Skew-Contaminated-Expected-Cauchy")
class.dist<-object$class
dim.y<-ncol(as.matrix(object$y))
names.dim.y<-ifelse(dim.y<=3,switch(dim.y, "1"="Univariate","2"="Bivariate", "3"="Trivariate"),"Multivariate")
perc<-length(which(object$Mahal>object$cut))/length(object$Mahal)
wh<-"none"
if(perc>0) {wh<- which(object$Mahal>object$cut);
wh1<-wh[1] 
if(length(wh)>1) for(jj in 2:length(wh)) wh1<-paste(wh1,"-",wh[jj],sep="")
}
cat("-------------------------------------------------------------------------\n")
cat("    Used function: ",object$"function","\n", sep = "")
cat("    Mahalanobis distance for the\n")
cat("   ",names.dim.y,name.dist,"Regression model\n")
cat("    Class:",class.dist,"\n")
cat("-------------------------------------------------------------------------\n")
cat("    observations:",length(object$Mahal),"\n")
cat("    cutoff point:",object$cut,"\n")
cat("    alpha: ",100*object$alpha,"% \n",sep="")
cat("    Number of points greater than cutoff: ",round(perc*length(object$Mahal))," (",round(100*perc,1),"%) \n",sep="")
if(perc==0) cat("    Observations greater than cutoff: ",wh,"\n",sep="")
if(perc>0) cat("    Observations greater than cutoff: \n    ",wh1,"\n",sep="")
cat("-------------------------------------------------------------------------\n")
}
else
{
if(is.null(object$se)) stop("Standard errors weren't be estimated: Numerical problems with the inversion of the information matrix\n
Summary can't be performed")
asterisk<-function (x) 
{
    if (x > 0.1) {
        ast = " "}
    else {
        if (x > 0.05) {
            ast = "."}
        else {
            if (x > 0.01) {
                ast = "*"}
            else {
                if (x > 0.001) {
                  ast = "**"}
                else {
                  {
                    ast = "***"}
                }
            }
        }
    }
    return(ast)
}
 name.dist<-switch(object$dist, "MN"="Normal", "MT"="t-Student", "MSL"="Slash", "MCN"="Contaminated Normal",
"MSN"="Skew-Normal","MSTN"="Skew-t-Normal", "MSSL"="Skew-Slash-Normal", "MSCN"="Skew-Contaminated-Normal",
"MSTT"="Skew-t", "MSSL2"="Skew-Slash", "MSCN2"="Skew-Contaminated-Normal",
"MSNC"="Skew-Normal-Cauchy","MSTEC"="Skew-t-Expectec-Cauchy", "MSSLEC"="Skew-slash-Expected-Cauchy", "MSCEC"="Skew-Contaminated-Expected-Cauchy")
       class.dist<-object$class
 aa<-sum(sapply(names(object$coefficients),substring, 1,5)=="alpha")
 dim.y<-(-1/2+sqrt(1/4+2*aa))
 bb<-which(sapply(names(object$coefficients),substring, 1,4)=="beta")
 tt<-cbind(object$coefficients[bb],object$se[bb],object$coefficients[bb]/object$se[bb],
pnorm(abs(object$coefficients[bb]/object$se[bb]),lower.tail=FALSE))
    ast = sapply(tt[,4],FUN=asterisk)
    tt = data.frame(round(tt, 5), ast)
 colnames(tt)<-c("estimate","s.e.","z value","Pr(>|z|)","")
 bb2<-which(sapply(names(object$coefficients),substring, 1,5)=="alpha"|sapply(names(object$coefficients),substring, 1,6)=="lambda")
tt2=round(cbind(object$coefficients[bb2],object$se[bb2]),5);colnames(tt2)<-c("estimate","s.e.")
if(object$dist=="MT" | object$dist=="MSL" | object$dist=="MSTN" | object$dist=="MSSL")
{bb3<-which(sapply(names(object$coefficients),substring, 1,2)=="nu");tt3=round(cbind(object$coefficients[bb3],object$se[bb3]),5);colnames(tt3)<-c("estimate","s.e.");rownames(tt3)<-c("nu")}
if(object$dist=="MSTT" | object$dist=="MSSL2" | object$dist=="MSTEC" | object$dist=="MSSLEC")
{tt3<-matrix(object$nu,ncol=1);colnames(tt3)<-"estimate";rownames(tt3)<-c("nu")}
      if(object$dist=="MCN" | object$dist=="MSCN")
{bb3<-which(sapply(names(object$coefficients),substring, 1,2)=="nu"|sapply(names(object$coefficients),substring, 1,5)=="gamma");tt3=round(cbind(object$coefficients[bb3],object$se[bb3]),5);colnames(tt3)<-c("estimate","s.e.");rownames(tt3)<-c("nu","gamma")}
      if(object$dist=="MSCN2" | object$dist=="MSCEC")
{tt3<-matrix(c(object$nu,object$gamma),ncol=1);colnames(tt3)<-"estimate";rownames(tt3)<-c("nu","gamma")}
tt4<-matrix(c(object$logLik,object$AIC, object$BIC),ncol=1)
rownames(tt4)<-c("logLik","AIC","BIC"); colnames(tt4)<-""
 names.dim.y<-ifelse(dim.y<=3,switch(dim.y, "1"="Univariate","2"="Bivariate", "3"="Trivariate"),"Multivariate")
 cat("\n")
  if(substring(object$"function",1,5)=="mback"){
cat("-------------------------------------------------------------------------\n")
cat("    Used function: ",object$"function","\n", sep = "")
cat("    Selecting covariates in the\n")
cat("   ",names.dim.y,name.dist,"Regression model\n")
cat("    Method:",ifelse(object$choose.crit=="sign",paste("significance (",100*object$significance,"%)",sep=""),object$choose.crit))
cat("\n    Class:",class.dist,"\n")
cat("-------------------------------------------------------------------------\n")
cat("  Summary for the Selected Regression Coefficients\n")
}
 if(substring(object$"function",1,8)=="estimate"){
       cat("-------------------------------------------------------------------------\n")
       cat("         Used function: ",object$"function","\n", sep = "")
       cat("        ",names.dim.y,name.dist,"Regression model\n")
       cat("         Class:",class.dist,"\n")
       cat("-------------------------------------------------------------------------\n")
       cat("\n")
 cat("         Regression Coefficients\n")
}
 if(substring(object$"function",1,7)=="choose."){
   if(object$"function"=="choose.models") class.dist<-c("MSMN, MSSMN, MSMSN and MSMSNC")
       cat("-------------------------------------------------------------------------\n")
       cat("    Used function: ",object$"function","\n", sep = "")
       cat("    Selecting a distribution in the\n")
       cat("   ",names.dim.y,class.dist,"Class of distributions\n")
       cat("    Method:",ifelse(object$choose.crit=="sign",paste("significance (",100*object$significance,"%)\n",sep=""),paste(object$choose.crit,"\n")))
       cat("    Distribution:",name.dist,"\n")
       cat("-------------------------------------------------------------------------\n")
       cat("\n")
cat("  Summary for the Selected Regression Coefficients\n")
}
 if(substring(object$"function",1,7)=="choose2"){
   if(object$fitted.classes=="ALL") class.dist<-c("MSMN, MSSMN, MSMSN and MSMSNC")
       cat("-------------------------------------------------------------------------\n")
       cat("    Used function: ",object$"function","\n", sep = "")
       cat("    Selecting a distribution in the\n")
       cat("   ",names.dim.y,class.dist,"Class of distributions\n")
       cat("    and then, select covariates\n")
       cat("    Method to select a distribution:",ifelse(object$choose.crit=="sign",paste("significance (",100*object$significance,"%)\n",sep=""),paste(object$choose.crit,"\n")))
       cat("    Method to select covariates:",ifelse(object$choose.crit=="sign",paste("significance (",100*object$significance,"%)\n",sep=""),paste(object$choose.crit,"\n")))
       cat("    Class:",class.dist,"\n")
       cat("-------------------------------------------------------------------------\n")
       cat("\n")
cat("  Summary for the Selected Regression Coefficients\n")
}
       cat("\n")
       print(tt)
       cat("---\n")
       cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
       cat("\n")
       cat("-------------------------------------------------------------------------\n")
       if(class.dist=="MSMN") cat("        Scale Coefficients\n") else cat("        Scale and Assymetry Coefficients\n")
       cat("\n")
       print(tt2)
       cat("-------------------------------------------------------------------------\n")
 if(object$dist=="MT" | object$dist=="MSL" | object$dist=="MSTN" | object$dist=="MSSL" |  object$dist=="MSTT" | object$dist=="MSSL2" | object$dist=="MSTEC" | object$dist=="MSSLEC")
{cat("        nu Coefficient\n");cat("\n");print(tt3)}
 if(object$dist=="MCN" | object$dist=="MSCN" | object$dist=="MSCN2" | object$dist=="MSCEC")
{cat("        nu and gamma Coefficients\n");cat("\n");print(tt3)}
if(object$dist!="MN" & object$dist!="MSN" & object$dist!="MSNC")
     {cat("-------------------------------------------------------------------------\n")}
cat("        log-Likelihood function, AIC and BIC criteria\n")
       print(tt4)
       cat("-------------------------------------------------------------------------\n")
}
}
