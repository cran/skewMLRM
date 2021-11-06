mbacksign <-
function (y, X = NULL, max.iter = 1000, prec = 1e-04, dist = "MN", 
    significance = 0.05, ...) 
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
    if (!is.numeric(significance)) 
        stop("significance must be numeric")
    if (significance <= 0) 
        stop("significance must be a positive value")
    if (significance > 0.1) 
        stop("significance is usually less than or equal to 0.10")
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
    m = sum(mapply(sum.ones, X))/(length(X) * nrow(X[[1]]))
    q = ncol(y)
    if ((ncol(X[[1]]) - m) == q) 
        stop("X have only the intercept term(s)")
    if (ncol(y) != nrow(X[[1]])) 
        stop("y does not have the same number of columns than X")
    if (nrow(y) != length(X)) 
        stop("y does not have the same number of observations than X")
    tt <- table(unlist(sapply(X, pos.ones)))
    ind.interc <- as.numeric(names(tt)[which(tt == nrow(y))])
    if (length(ind.interc) != q) 
        stop("X should be intercept term for each variable")
    m = ncol(X[[1]])
    z.critico <- qnorm(1 - significance/2)
    estimate.dist <- get(paste("estimate.", dist, sep = ""), 
        mode = "function")
    if (dist != "MSTEC" & dist != "MSSLEC" & dist != 
        "MSCEC" & dist != "MSTT" & dist != "MSSL2" & 
        dist != "MSCN2") 
        object.1 <- estimate.dist(y = y, X = X, max.iter = max.iter, 
            prec = prec, est.var = TRUE)
    if (dist == "MSTEC" | dist == "MSSLEC" | dist == 
        "MSTT" | dist == "MSSL2") 
        object.1 <- estimate.dist(y = y, X = X, max.iter = max.iter, 
            prec = prec, est.var = TRUE, nu.fixed = nu.fixed, 
            nu.min = nu.min)
    if (dist == "MSCEC" | dist == "MSCN2") 
        object.1 <- estimate.dist(y = y, X = X, max.iter = max.iter, 
            prec = prec, est.var = TRUE, nu.fixed = nu.fixed, 
            gamma.fixed = gamma.fixed)
    beta.est <- object.1$coefficients[c(1:m)[-ind.interc]]
    if (dist != "MSTEC" & dist != "MSSLEC" & dist != 
        "MSCEC" & dist != "MSTT" & dist != "MSSL2" & 
        dist != "MSCN2") 
        se <- try(se.est(object.1$coefficients, y, X, dist = dist), 
            silent = TRUE)
    if (dist == "MSTEC" | dist == "MSSLEC" | dist == 
        "MSTT" | dist == "MSSL2") 
        se <- try(se.est(object.1$coefficients, y, X, dist = dist, 
            nu = object.1$nu), silent = TRUE)
    if (dist == "MSCEC" | dist == "MSCN2") 
        se <- try(se.est(object.1$coefficients, y, X, dist = dist, 
            nu = object.1$nu, gamma = object.1$gamma), silent = TRUE)
    if (!is.numeric(se)) 
        stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
    ep.beta <- se[c(1:m)[-ind.interc]]
    xmenos.multi <- X
    dim.beta <- length(beta.est)
    test.t <- abs(beta.est/ep.beta)
    lista.final <- list(coefficients = object.1$coefficients, 
        logLik = object.1$logLik, AIC = object.1$AIC, BIC = object.1$BIC, 
        conv = object.1$conv, dist = object.1$dist, class = object.1$class, 
        comment = "The final model considered all the betas")
    lista.final$dist = object.1$dist
    lista.final$class = object.1$class
    lista.final$X <- X
    if (dist == "MSTEC" | dist == "MSSLEC" | dist == 
        "MSTT" | dist == "MSSL2") 
        lista.final$nu = object.1$nu
    if (dist == "MSCEC" | dist == "MSCN2") 
        lista.final$nu = object.1$nu
    lista.final$gamma = object.1$gamma
    criterio <- sum(test.t - z.critico > 0, na.rm = TRUE)
    if (criterio > 0) {
        t.min <- c(1:m)[-ind.interc][which.min(test.t)]
        historico <- c()
        cont <- sum(is.na(test.t))
        while (criterio < (length(test.t) - cont) && (dim.beta > 
            1)) {
            historico <- c(historico, names(object.1$coefficients)[1:m][t.min])
            names.betas <- rownames(object.1$coefficients)[1:m][-(t.min)]
            xmenos.multi <- lapply(xmenos.multi, function(x) {
                x[, -(t.min), drop = FALSE]
            })
            m = ncol(xmenos.multi[[1]])
            tt <- table(unlist(sapply(xmenos.multi, pos.ones)))
            ind.interc <- as.numeric(names(tt)[which(tt == nrow(y))])
            if (dist != "MSTEC" & dist != "MSSLEC" & 
                dist != "MSCEC" & dist != "MSTT" & 
                dist != "MSSL2" & dist != "MSCN2") 
                object.1 <- estimate.dist(y = y, X = xmenos.multi, 
                  max.iter = max.iter, prec = prec, est.var = TRUE)
            if (dist == "MSTEC" | dist == "MSSLEC" | 
                dist == "MSTT" | dist == "MSSL2") 
                object.1 <- estimate.dist(y = y, X = xmenos.multi, 
                  max.iter = max.iter, prec = prec, est.var = TRUE, 
                  nu.fixed = nu.fixed, nu.min = nu.min)
            if (dist == "MSCEC" | dist == "MSCN2") 
                object.1 <- estimate.dist(y = y, X = xmenos.multi, 
                  max.iter = max.iter, prec = prec, est.var = TRUE, 
                  nu.fixed = nu.fixed, gamma.fixed = gamma.fixed)
            rownames(object.1$coefficients)[1:m] <- names.betas
            beta.est <- object.1$coefficients[c(1:m)[-ind.interc]]
            if (dist != "MSTEC" & dist != "MSSLEC" & 
                dist != "MSCEC" & dist != "MSTT" & 
                dist != "MSSL2" & dist != "MSCN2") 
                se <- try(se.est(object.1$coefficients, y, X = xmenos.multi, 
                  dist = dist), silent = TRUE)
            if (dist == "MSTEC" | dist == "MSSLEC" | 
                dist == "MSTT" | dist == "MSSL2") 
                se <- try(se.est(object.1$coefficients, y, X = xmenos.multi, 
                  dist = dist, nu = object.1$nu), silent = TRUE)
            if (dist == "MSCEC") 
                se <- try(se.est(object.1$coefficients, y, X = xmenos.multi, 
                  dist = dist, nu = object.1$nu, gamma = object.1$gamma), 
                  silent = TRUE)
            if (!is.numeric(se)) 
                stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
            ep.beta <- se[c(1:m)[-ind.interc]]
            dim.beta <- length(beta.est)
            test.t <- abs(beta.est/ep.beta)
            t.min <- c(1:m)[-ind.interc][which.min(test.t)]
            criterio <- sum(test.t - z.critico > 0, na.rm = TRUE)
            lista.final <- list(coefficients = object.1$coefficients)
            if (dist == "MSTEC" | dist == "MSSLEC" | 
                dist == "MSTT" | dist == "MSSL2") 
                lista.final$nu = object.1$nu
            if (dist == "MSCEC" | dist == "MSCN2") 
                lista.final$nu = object.1$nu
            lista.final$gamma = object.1$gamma
            lista.final$logLik = object.1$logLik
            lista.final$AIC = object.1$AIC
            lista.final$BIC = object.1$BIC
            lista.final$dist = object.1$dist
            lista.final$class = object.1$class
            lista.final$conv = object.1$conv
            lista.final$comment = paste("The final model eliminated", 
                length(historico), "betas")
            lista.final$eliminated = historico
            lista.final$X = xmenos.multi
        }
    }
    se <- "Error"
    fit <- lista.final
  if(!is.null(fit$eliminated))
{names.betas=paste("beta",(1:(ncol(X[[1]])))[-as.numeric(substring(fit$eliminated,5,10))],sep="");
 names(fit$coefficients)[1:length(names.betas)]=names.betas}
    if (fit$dist != "MSTEC" & fit$dist != "MSSLEC" & 
        fit$dist != "MSCEC" & fit$dist != "MSTT" & 
        fit$dist != "MSSL2" & dist != "MSCN2") 
        se <- try(se.est(P = fit$coefficients, y = y, X = fit$X, 
            dist = fit$dist), silent = TRUE)
    if (fit$dist == "MSTEC" | fit$dist == "MSSLEC" | 
        fit$dist == "MSTT" | fit$dist == "MSSL2") 
        se <- try(se.est(P = fit$coefficients, y = y, X = fit$X, 
            dist = fit$dist, nu = fit$nu), silent = TRUE)
    if (fit$dist == "MSCEC" | fit$dist == "MSCN2") 
        se <- try(se.est(P = fit$coefficients, y = y, X = fit$X, 
            dist = fit$dist, nu = fit$nu, gamma = fit$gamma), 
            silent = TRUE)
    names(se) <- names(fit$coefficients)
    names(se)[1:length(names.betas)]=names.betas
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
            RVAL <- list(coefficients = fit$coefficients, se = se, 
                nu = fit$nu, gamma = fit$gamma, logLik = fit$logLik, 
                AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                conv = fit$conv, dist = fit$dist, class = fit$class)
    }
    if (fit$class == "MSMN" | fit$class == "MSSMN") {
        RVAL <- list(coefficients = fit$coefficients, se = se, 
            logLik = fit$logLik, AIC = fit$AIC, BIC = fit$BIC, 
            iterations = fit$iterations, conv = fit$conv, dist = fit$dist, 
            class = fit$class)
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
                BIC = fit$BIC, iterations = fit$iterations, conv = fit$conv, 
                dist = fit$dist, class = fit$class)
        }
        else {
            RVAL <- list(coefficients = fit$coefficients, se = se, 
                nu = fit$nu, gamma = fit$gamma, logLik = fit$logLik, 
                AIC = fit$AIC, BIC = fit$BIC, iterations = fit$iterations, 
                conv = fit$conv, dist = fit$dist, class = fit$class)
        }
    }
     class(RVAL) <- "skewMLRM"
    RVAL$choose.crit <- "sign"
    RVAL$comment <- fit$comment
    if (!is.null(fit$eliminated)) 
        RVAL$eliminated <- fit$eliminated
    RVAL$y <- y
    RVAL$X <- fit$X
    RVAL$significance <- significance
    RVAL$"function" <- "mbacksign"
    RVAL
}
