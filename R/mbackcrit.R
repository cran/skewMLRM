mbackcrit <-
function (y, X = NULL, max.iter = 1000, prec = 1e-04, dist = "MN", 
    criteria = "AIC", est.var = TRUE, cluster = FALSE, 
    ...) 
{
    se.est <- function(P, y, X, dist = "MN", nu = 3, gamma = 0.5) {
        if (!any(dist == c("MN", "MT", "MSL", 
            "MCN", "MSN", "MSNC", "MSTEC", 
            "MSTN", "MSTT", "MSSLEC", "MSSL", 
            "MSSL2", "MSCN", "MSCN2", "MSCEC"))) 
            stop("distribution is not recognized")
        y <- as.matrix(y)
        if (!is.matrix(y)) 
            stop("y must have at least one element")
        if (is.null(X)) {
            X <- array(c(diag(ncol(y))), c(ncol(y), ncol(y), 
                nrow(y)))
        }
        if (is.array(X) == FALSE & is.list(X) == FALSE) 
            stop("X must be an array or a list")
        if (is.array(X)) {
            Xs <- list()
            if (ncol(y) > 1 | !is.matrix(X)) {
                for (i in 1:nrow(y)) {
                  Xs[[i]] <- matrix(t(X[, , i]), nrow = ncol(y))
                }
            }
            if (ncol(y) == 1 & is.matrix(X)) {
                for (i in 1:nrow(y)) {
                  Xs[[i]] <- matrix(t(X[i, ]), nrow = 1)
                }
            }
            X <- Xs
        }
        if (ncol(y) != nrow(X[[1]])) 
            stop("y does not have the same number of columns than X")
        if (nrow(y) != length(X)) 
            stop("y does not have the same number of observations than X")
        FI.dist <- get(paste("FI.", dist, sep = ""), 
            mode = "function")
        if (dist == "MSTT" | dist == "MSSL2" | dist == 
            "MSTEC" | dist == "MSSLEC") {
            if (!is.numeric(nu)) 
                stop("nu must be a number greater than 1")
            if (as.numeric(nu) < 1) 
                stop("nu must be a number greater than 1")
        }
        if (dist == "MSCN2" | dist == "MSCEC") {
            if (!is.numeric(nu)) 
                stop("nu must be a number between 0 and 1")
            if ((as.numeric(nu) < 0 | as.numeric(nu) > 1)) 
                stop("nu must be a number between 0 and 1")
            if (!is.numeric(gamma)) 
                stop("gamma must be a number between 0 and 1")
            if (as.numeric(gamma) < 0 | as.numeric(gamma) > 1) 
                stop("gamma must be a number between 0 and 1")
        }
        P <- matrix(P, ncol = 1)
        if (!any(dist == c("MSTEC", "MSSLEC", "MSCEC", 
            "MSTT", "MSSL2", "MSCN2"))) 
            MI.obs <- FI.dist(P = P, y = y, X = X)
        if (dist == "MSTEC" | dist == "MSSLEC" | 
            dist == "MSTT" | dist == "MSSL2") 
            MI.obs <- FI.dist(P = P, y = y, X = X, nu = nu)
        if (dist == "MSCEC" | dist == "MSCN2") 
            MI.obs <- FI.dist(P = P, y = y, X = X, nu = nu, gamma = gamma)
        test = try(solve2(MI.obs), silent = TRUE)
        se = c()
        if (is.numeric(test) & max(diag(test)) < 0) {
            se = sqrt(-diag(test))
        }
        else stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
        se
    }
    if (!any(dist == c("MN", "MT", "MSL", "MCN", 
        "MSN", "MSNC", "MSTEC", "MSTN", 
        "MSTT", "MSSLEC", "MSSL", "MSSL2", 
        "MSCN", "MSCN2", "MSCEC"))) 
        stop("distribution is not recognized")
    if (cluster != TRUE) 
        cluster = FALSE
    if (!(criteria == "AIC" | criteria == "BIC")) 
        stop("AIC or BIC criterion must be chosen")
    sum.ones <- function(x) sum(x == 1)
    pos.ones <- function(x) which(apply(x, 2, sum) == 1)
    y <- as.matrix(y)
    if (!is.matrix(y)) 
        stop("y must have at least one element")
    if (is.null(X)) {
        X <- array(c(diag(ncol(y))), c(ncol(y), ncol(y), nrow(y)))
    }
    if (is.array(X) == FALSE & is.list(X) == FALSE) 
        stop("X must be an array or a list")
    if (is.array(X)) {
        Xs <- list()
        if (ncol(y) > 1 | !is.matrix(X)) {
            for (i in 1:nrow(y)) {
                Xs[[i]] <- matrix(t(X[, , i]), nrow = ncol(y))
            }
        }
        if (ncol(y) == 1 & is.matrix(X)) {
            for (i in 1:nrow(y)) {
                Xs[[i]] <- matrix(t(X[i, ]), nrow = 1)
            }
        }
        X <- Xs
    }
    m = sum(mapply(sum.ones, X))/length(X)
    M=ncol(X[[1]])
    q = ncol(y)
    if ((ncol(X[[1]]) - m) == q) 
        stop("X have only the intercept term(s)")
    if (ncol(y) != nrow(X[[1]])) 
        stop("y does not have the same number of columns than X")
    if (nrow(y) != length(X)) 
        stop("y does not have the same number of observations than X")
    mf <- match.call(expand.dots = FALSE)
    if ("nu.fixed" %in% names(mf$...)) 
        nu.fixed = mf$...$nu.fixed
    if ("gamma.fixed" %in% names(mf$...)) 
        gamma.fixed = mf$...$gamma.fixed
    if ("nu.min" %in% names(mf$...)) 
        nu.min = mf$...$nu.min
    if (dist == "MSTT" | dist == "MSSL2" | dist == 
        "MSTEC" | dist == "MSSLEC") {
        if (!exists("nu.fixed")) 
            nu.fixed = 3
        if (!exists("nu.min")) 
            nu.min = 2.0001
        if (!exists("nu.fixed")) 
            nu.fixed = TRUE
        if (!is.numeric(nu.min) | nu.min <= 0) 
            stop("nu.min should be a positive number")
        if (nu.fixed != FALSE & !is.numeric(nu.fixed)) 
            stop("nu fixed must be a number greater than 1")
        if (is.numeric(nu.fixed) & as.numeric(nu.fixed) < 1) 
            stop("nu fixed must be a number greater than 1")
    }
    if (dist == "MSCN2" | dist == "MSCEC") {
        if (!exists("nu.fixed")) 
            nu.fixed = 0.1
        if (!exists("gamma.fixed")) 
            gamma.fixed = 0.5
        if (nu.fixed != FALSE & !is.numeric(nu.fixed)) 
            stop("nu fixed must be a number between 0 and 1")
        if (is.numeric(nu.fixed) & (as.numeric(nu.fixed) < 0 | 
            as.numeric(nu.fixed) > 1)) 
            stop("nu fixed must be a number between 0 and 1")
        if (gamma.fixed != FALSE & !is.numeric(gamma.fixed)) 
            stop("gamma fixed must be a number between 0 and 1")
        if (is.numeric(gamma.fixed) & (as.numeric(gamma.fixed) < 
            0 | as.numeric(gamma.fixed) > 1)) 
            stop("gamma fixed must be a number between 0 and 1")
    }
    estimate.dist <- get(paste("estimate.", dist, sep = ""), 
        mode = "function")
    X.complete <- X
    m = ncol(X.complete[[1]])
    if (dist != "MSTEC" & dist != "MSSLEC" & dist != 
        "MSCEC" & dist != "MSTT" & dist != "MSSL2" & 
        dist != "MSCN2") 
        fit.complete <- estimate.dist(y = y, X = X, max.iter = max.iter, 
            prec = prec, est.var = FALSE)
    if (dist == "MSTEC" | dist == "MSSLEC" | dist == 
        "MSTT" | dist == "MSSL2") 
        fit.complete <- estimate.dist(y = y, X = X, max.iter = max.iter, 
            prec = prec, est.var = FALSE, nu.fixed = nu.fixed, 
            nu.min = nu.min)
    if (dist == "MSCEC" | dist == "MSCN2") 
        fit.complete <- estimate.dist(y = y, X = X, max.iter = max.iter, 
            prec = prec, est.var = FALSE, nu.fixed = nu.fixed, 
            gamma.fixed = gamma.fixed)
    names.betas <- rownames(fit.complete$coefficients)[1:m]
    tt <- table(unlist(sapply(X, pos.ones)))
    ind.interc <- as.numeric(names(tt)[which(tt == nrow(y))])
    if (length(ind.interc) != q) 
        stop("X should be intercept term for each variable")
    if (criteria == "AIC") 
        criteria.complete <- fit.complete$AIC
    if (criteria == "BIC") 
        criteria.complete <- fit.complete$BIC
    lista.final <- list(coefficients = fit.complete$coefficients, 
        logLik = fit.complete$logLik, AIC = fit.complete$AIC, 
        BIC = fit.complete$BIC, conv = fit.complete$conv, dist = fit.complete$dist, 
        class = fit.complete$class, comment = "The final model considered all the betas")
    if (dist == "MSTEC" | dist == "MSSLEC" | dist == 
        "MSTT" | dist == "MSSL2") 
        lista.final$nu = fit.complete$nu
    if (dist == "MSCEC" | dist == "MSCN2") 
        lista.final$nu = fit.complete$nu
    lista.final$gamma = fit.complete$gamma
    historico <- c()
    XICant <- ifelse(criteria == "AIC", fit.complete$AIC, fit.complete$BIC)
    fit.complete00 <- fit.complete
    ind.interc0 <- ind.interc
    t.min=NULL
    flag <- FALSE
    while (!flag) {
          fit.complete0 <-NULL
        if (cluster == FALSE | any(dist == c("MN", "MT", 
            "MSL", "MCN", "MSN", "MSTN", 
            "MSTT", "MSSL", "MSSL2", "MSCN", 
            "MSCN2"))) {
            criteria.complete <- c()
            for (j in 1:length(c(1:m)[-ind.interc])) {
                X.complete.j <- lapply(X.complete, function(x) {
                  x[, -(c(1:m)[-ind.interc][j]), drop = FALSE]
                })
                if (dist != "MSTEC" & dist != "MSSLEC" & 
                  dist != "MSCEC" & dist != "MSTT" & 
                  dist != "MSSL2" & dist != "MSCN2") 
                  fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                    max.iter = max.iter, prec = prec, est.var = FALSE)
                if (dist == "MSTEC" | dist == "MSSLEC" | 
                  dist == "MSTT" | dist == "MSSL2") 
                  fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                    max.iter = max.iter, prec = prec, est.var = FALSE, 
                    nu.fixed = nu.fixed, nu.min = nu.min)
                if (dist == "MSCEC" | dist == "MSCN2") 
                  fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                    max.iter = max.iter, prec = prec, est.var = FALSE, 
                    nu.fixed = nu.fixed, gamma.fixed = gamma.fixed)
                if (criteria == "AIC") 
                  criteria.complete[j] <- fit.complete.j$AIC
                if (criteria == "BIC") 
                  criteria.complete[j] <- fit.complete.j$BIC
                  fit.complete0[[j]]=fit.complete.j
            }
        }
        else {
            cl0 <- parallel::detectCores() - 1
            cl <- parallel::makeCluster(cl0)
            doParallel::registerDoParallel(cl)
            if (criteria == "AIC") {
                j = 1
                criteria.complete <- foreach(j = 1:length(c(1:m)[-ind.interc]), 
                  .combine = rbind) %dopar% {
                  require(skewMLRM)
                  X.complete.j <- lapply(X.complete, function(x) {
                    x[, -(c(1:m)[-ind.interc][j]), drop = FALSE]
                  })
                  if (dist != "MSTEC" & dist != "MSSLEC" & 
                    dist != "MSCEC" & dist != "MSTT" & 
                    dist != "MSSL2" & dist != "MSCN2") 
                    fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                      max.iter = max.iter, prec = prec, est.var = FALSE)
                  if (dist == "MSTEC" | dist == "MSSLEC" | 
                    dist == "MSTT" | dist == "MSSL2") 
                    fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                      max.iter = max.iter, prec = prec, est.var = FALSE, 
                      nu.fixed = nu.fixed, nu.min = nu.min)
                  if (dist == "MSCEC" | dist == "MSCN2") 
                    fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                      max.iter = max.iter, prec = prec, est.var = FALSE, 
                      nu.fixed = nu.fixed, gamma.fixed = gamma.fixed)
                  criteria.complete <- fit.complete.j$AIC
                }
                criteria.complete <- c(criteria.complete)
            }
            if (criteria == "BIC") {
                j = 1
                criteria.complete <- foreach(j = 1:length(c(1:m)[-ind.interc]), 
                  .combine = rbind) %dopar% {
                  require(skewMLRM)
                  X.complete.j <- lapply(X.complete, function(x) {
                    x[, -(c(1:m)[-ind.interc][j]), drop = FALSE]
                  })
                  if (dist != "MSTEC" & dist != "MSSLEC" & 
                    dist != "MSCEC" & dist != "MSTT" & 
                    dist != "MSSL2" & dist != "MSCN2") 
                    fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                      max.iter = max.iter, prec = prec, est.var = FALSE)
                  if (dist == "MSTEC" | dist == "MSSLEC" | 
                    dist == "MSTT" | dist == "MSSL2") 
                    fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                      max.iter = max.iter, prec = prec, est.var = FALSE, 
                      nu.fixed = nu.fixed, nu.min = nu.min)
                  if (dist == "MSCEC" | dist == "MSCN2") 
                    fit.complete.j <- estimate.dist(y = y, X = X.complete.j, 
                      max.iter = max.iter, prec = prec, est.var = FALSE, 
                      nu.fixed = nu.fixed, gamma.fixed = gamma.fixed)
                  criteria.complete <- fit.complete.j$BIC
                }############# fit.complete0 aqui fica complicado de salvar.
                criteria.complete <- c(criteria.complete)
            }
            parallel::stopCluster(cl)
        }
        ind = c(1:m)[-ind.interc][which.min(criteria.complete)]
        indj=which.min(criteria.complete)
        t.min <- c(t.min,ind)
        #crit.complete <- ifelse(criteria == "AIC", fit.complete$AIC, fit.complete$BIC)
        #crit.complete <- min(criteria.complete)
        #if (min(criteria.complete) >= crit.complete) {
        if (min(criteria.complete) >= XICant) { flag <- TRUE}
        if (min(criteria.complete) < XICant) {
            XICant <- min(criteria.complete)
            X.complete <- lapply(X.complete, function(x) {
                x[, -ind, drop = FALSE]
            })
            fit.complete <- fit.complete0[[indj]]
            names.betas <- names(fit.complete00$coefficients)[1:m][-ind] 
            tt <- table(unlist(sapply(X, pos.ones)))
            ind.interc <- as.numeric(names(tt)[which(tt == nrow(y))]) 
            m = ncol(X.complete[[1]]) 
        }
    }
            if(length(t.min)>1){ t.min<-t.min[-length(t.min)]; historico <- names(fit.complete00$coefficients)[1:M][t.min]}
if(length(t.min)==1){t.min<-NULL; historico<-NULL}
            lista.final <- list(coefficients = fit.complete$coefficients)
            if (dist == "MSTEC" | dist == "MSSLEC" | 
                dist == "MSTT" | dist == "MSSL2") 
            {lista.final$nu = fit.complete$nu}
            if (dist == "MSCEC" | dist == "MSCN2")
            {lista.final$nu = fit.complete$nu; lista.final$gamma = fit.complete$gamma}
            lista.final$logLik = fit.complete$logLik
            lista.final$AIC = fit.complete$AIC
            lista.final$BIC = fit.complete$BIC
            lista.final$conv = fit.complete$conv
            lista.final$dist = fit.complete$dist
            lista.final$class = fit.complete$class
            lista.final$comment = paste("The final model eliminated", 
                length(historico), "betas")
            lista.final$eliminated = historico
            lista.final$X = fit.complete$X
            lista.final$y = fit.complete$y
if(!is.null(historico)) names(lista.final$coefficients)[1:m] <- names.betas

    se <- "Error"
    fit <- lista.final
    if (!is.null(fit$eliminated)) {
        names.betas = paste("beta", (1:(ncol(X[[1]])))[-as.numeric(substring(fit$eliminated, 
            5, 10))], sep = "")
        names(fit$coefficients)[1:length(names.betas)] = names.betas
    }
    if (est.var) {
        if (fit$dist != "MSTEC" & fit$dist != "MSSLEC" & 
            fit$dist != "MSCEC" & fit$dist != "MSTT" & 
            fit$dist != "MSSL2" & dist != "MSCN2") 
            se <- try(se.est(fit$coefficients, y, fit$X, dist = fit$dist), 
                silent = TRUE)
        if (fit$dist == "MSTEC" | fit$dist == "MSSLEC" | 
            fit$dist == "MSTT" | fit$dist == "MSSL2") 
            se <- try(se.est(fit$coefficients, y, lista.final$X, 
                dist = fit$dist, nu = fit$nu), silent = TRUE)
        if (fit$dist == "MSCEC" | fit$dist == "MSCN2") 
            se <- try(se.est(fit$coefficients, y, lista.final$X, 
                dist = fit$dist, nu = fit$nu, gamma = fit$gamma), 
                silent = TRUE)
    }
    if (!grepl("Error", se[1])) {
        names(se) <- names(fit$coefficients)
        if (!is.null(fit$eliminated)) 
            names(se)[1:length(names.betas)] = names.betas
        if (fit$class == "MSMSNC") {
            if (fit$dist != "MSCEC") {
                if (fit$dist != "MSNC") 
                  RVAL <- list(coefficients = fit$coefficients, 
                    se = se, nu = fit$nu, logLik = fit$logLik, 
                    AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                    conv = fit$conv, dist = fit$dist, class = fit$class)
                if (fit$dist == "MSNC") 
                  RVAL <- list(coefficients = fit$coefficients, 
                    se = se, logLik = fit$logLik, AIC = fit$AIC, 
                    BIC = fit$BIC, iterations = fit$iterations, 
                    conv = fit$conv, dist = fit$dist, class = fit$class)
            }
            if (fit$dist == "MSCEC") 
                RVAL <- list(coefficients = fit$coefficients, 
                  se = se, nu = fit$nu, gamma = fit$gamma, logLik = fit$logLik, 
                  AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                  conv = fit$conv, dist = fit$dist, class = fit$class)
        }
        if (fit$class == "MSMN" | fit$class == "MSSMN") {
            RVAL <- list(coefficients = fit$coefficients, se = se, 
                logLik = fit$logLik, AIC = fit$AIC, BIC = fit$BIC, 
                iterations = fit$iterations, conv = fit$conv, 
                dist = fit$dist, class = fit$class)
        }
        if (fit$class == "MSMSN") {
            if (fit$dist != "MSCN2") {
                if (fit$class != "MSN") 
                  RVAL <- list(coefficients = fit$coefficients, 
                    se = se, nu = fit$nu, logLik = fit$logLik, 
                    AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                    conv = fit$conv, dist = fit$dist, class = fit$class)
                else RVAL <- list(coefficients = fit$coefficients, 
                  se = se, logLik = fit$logLik, AIC = fit$AIC, 
                  BIC = fit$BIC, iterations = fit$iterations, 
                  conv = fit$conv, dist = fit$dist, class = fit$class)
            }
            else {
                RVAL <- list(coefficients = fit$coefficients, 
                  se = se, nu = fit$nu, gamma = fit$gamma, logLik = fit$logLik, 
                  AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                  conv = fit$conv, dist = fit$dist, class = fit$class)
            }
        }
    }
    else {
        if (fit$class == "MSMSNC") {
            if (fit$dist != "MSCEC") {
                if (fit$dist != "MSNC") 
                  RVAL <- list(coefficients = fit$coefficients, 
                    nu = fit$nu, logLik = fit$logLik, AIC = fit$AIC, 
                    BIC = fit$BIC, iterations = fit$iterations, 
                    conv = fit$conv, dist = fit$dist, class = fit$class)
                if (fit$dist == "MSNC") 
                  RVAL <- list(coefficients = fit$coefficients, 
                    logLik = fit$logLik, AIC = fit$AIC, BIC = fit$BIC, 
                    iterations = fit$iterations, conv = fit$conv, 
                    dist = fit$dist, class = fit$class)
            }
            if (fit$dist == "MSCEC") 
                RVAL <- list(coefficients = fit$coefficients, 
                  nu = fit$nu, gamma = fit$gamma, logLik = fit$logLik, 
                  AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                  conv = fit$conv, dist = fit$dist, class = fit$class)
        }
        if (fit$class == "MSMN" | fit$class == "MSSMN") {
            RVAL <- list(coefficients = fit$coefficients, logLik = fit$logLik, 
                AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                conv = fit$conv, dist = fit$dist, class = fit$class)
        }
        if (fit$class == "MSMSN") {
            if (fit$dist != "MSCN2") {
                if (fit$dist != "MSN") 
                  RVAL <- list(coefficients = fit$coefficients, 
                    nu = fit$nu, logLik = fit$logLik, AIC = fit$AIC, 
                    BIC = fit$BIC, iterations = fit$iterations, 
                    conv = fit$conv, dist = fit$dist, class = fit$class)
                else RVAL <- list(coefficients = fit$coefficients, 
                  logLik = fit$logLik, AIC = fit$AIC, BIC = fit$BIC, 
                  iterations = fit$iterations, conv = fit$conv, 
                  dist = fit$dist, class = fit$class)
            }
            else {
                RVAL <- list(coefficients = fit$coefficients, 
                  nu = fit$nu, gamma = fit$gamma, logLik = fit$logLik, 
                  AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                  conv = fit$conv, dist = fit$dist, class = fit$class)
            }
        }
    }
    class(RVAL) <- "skewMLRM"
    RVAL$choose.crit <- criteria
    RVAL$comment <- fit$comment
    if (!is.null(fit$eliminated)) 
        RVAL$eliminated <- fit$eliminated
    RVAL$y <- y
    RVAL$X <- fit$X
    RVAL$"function" <- "mbackcrit"
    RVAL
}
