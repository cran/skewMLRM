### 1.1

* option to maintain fixed nu in the following functions: estimate.MSTT, estimate.MSSL2, 
  estimate.MSCN2, estimate.MSTEC, estimate.MSSLEC and estimate.MSCEC.
* added function distMahal
* fixed bugs in mbacksign
* functions mbacksign and mbackcrit require intercept terms in the covariates matrix
* fixed bugs for the univariate case in all the functions estimate.xxx
* functions of the form estimate.xxx, choose.xxx, choose2, mbackcrit, mbacksign and
  distMahal produce object of the class "skewMLRM". Those functions also provide
  new elements in the returned list such as "dist", "class", "y", "X" and 
  "fitted classes".
* S3 method avaliable for summary, print and plot functions for object in the class 
  "skewMLRM"
* functions se.est and plotMahal were deprecated

