#' Run tests and produce plots required for this analysis
#'
#' @param Model A string or vector of strings, determining which emulator model to run. The entries for this argument should be some combination of of \code{c('lorenz63', 'lorenz96', 'ODEperiodic')}. The defualt value is \code{NULL}, all three models are run.
#' @param Ix The index of the model to be the focus of the emulator uncertainty analysis. For \code{Model = 'lorenz63'}, options for \code{Ix} are 1, 2, 3; for \code{Model = 'lorenz96'}, options are 1, 4, 8, 12, and for \code{Model = 'ODEperiodic'}, the options for \code{Ix} are 1, 2, 3, 4. If \code{NULL}, then all the options are run for the simulation model specified by \code{Model}.
#' @param outFile A string specifying the file path to output the standard error results to. If \code{NULL} (the default), the progress will be printed on the R console.
#' @param plotsFolder A string specifying the folder to save the plots file to. If unspecified, it defaults to the current working directory. 
#' @return Nothing is returned
#'
#' @importFrom deSolve lsode
#' @importFrom MASS rlm
#' @importFrom utils head tail
#' @importFrom stats predict smooth.spline
#'
#' @export


runTests <- function(Model = NULL,
                     Ix = NULL,
                     outFile = NULL,
                     plotsFolder = getwd()){

    ## check that argument Model takes one of the permittted values
    if(!is.null(Model)){
        if(!all(Model %in% c('lorenz63','lorenz96','ODEperiodic'))){
            stop("The argument 'Model' should be a vector with values taken from c('lorenz63','lorenz96','ODEperiodic')")
        }
    }
    ## check that argument Ix takes one of the permitted values
    if(!is.null(Ix)){
        ## only works if Model is specified - if not, override
        if(is.null(Model)){
            Ix <- NULL
            warning('Ix was non-NULL, but Model was NULL.... Overriding Ix!')
        } else {
            ## It only works if we are dealing with a single value of Model
            if(length(Model) != 1){
                stop("If arguments Model and Ix are specified (i.e. non-NULL), Model should be of length 1")
            }
            ## Check that the values given by Model and Ix are compatible
            if(Model == 'lorenz63'){
                if(!all(Ix %in% 1:3)){
                    stop("The 'lorenz63' model has indices 1,2,3 - try a different value of Ix")
                }
            } else if(Model == 'lorenz96'){
                if(!all(Ix %in% c(1,4,8,12))){
                    stop("The 'lorenz96' model has indices 1,4,8.12 - try a different value of Ix")
                }
            } else {
                if(!all(Ix %in% 1:4)){
                    stop("The 'ODEperiodic' model has indices 1,2,3,4 - try a different value of Ix")
                }
            }
        }
    }
    
    
    requireNamespace('deSolve')
    requireNamespace('MASS') 
    requireNamespace('utils') 
    requireNamespace('stats')
   

    if(!is.null(outFile)){
        sink(file = outFile)
    }

    ## {{{ Set up the models

    deterministicModelData <- list()

    ## generate test data for the Lorenz-63 model
    nx <- 3
    xinit = sin(1:nx)
    dt = 0.01
    ntime <- 1
    nparams <- 3
    nstepInit <- nstep <- as.integer(ntime/dt)
    nrealisations <- 1e4
    iorder <- sample(nrealisations)
    lorenz63data <- matrix(nrow = nrealisations,
                           ncol = nparams + nx*2)
    itime <- seq(1,nstep+1,by = 100)
    ## only store a subset of the time-series
    fullTimeSeries <- array(dim = c(nrealisations,nx,length(itime)))
    ##
    set.seed(42)
    ## set up the initial conditions
    xin <- runif(nx,0.5,4)
    xin <- c(1,0,0)
    params <- rgamma(nparams,2,2)
    params <- c(10, 28, 8/3)
    res <- .Fortran('model_lorenz63',
                    x = xin,
                    xout = double(nx*(nstep+1)),
                    dt = dt,
                    nstep = nstep,
                    params = params,
                    PACKAGE = 'EmulatorAnalysis')
    ##
    xin <- res$x
    for(i in 1:nrealisations){
        params <- rgamma(nparams,2,2)
        params <- c(10, 28, 8/3) + rnorm(3)*c(1,1,0.1)
        res <- .Fortran('model_lorenz63',
                        x = xin,
                        xout = double(nx*(nstep+1)),
                        dt = dt,
                        nstep = nstep,
                        params = params,
                        PACKAGE = 'EmulatorAnalysis')
        ## 
        xout <- res[['x']]
        lorenz63data[i,] <- c(params,xin,xout)
        fullTimeSeries[i,,] <- matrix(res$xout,nx,nstep+1)[,itime]
        ##
        res <- .Fortran('model_lorenz63',
                        x = xout,
                        xout = double(nx*(nstepInit+1)),
                        dt = dt,
                        nstep = nstepInit,
                        params = params,
                        PACKAGE = 'EmulatorAnalysis')
        xin <- res$x
    }
    if(any(!is.finite(lorenz63data))) stop('Ouch!')
    ## store the data in a list with contents for all the models
    deterministicModelData[['lorenz63']] <- list(
        params = lorenz63data[,1:nparams],
        xin = lorenz63data[,nparams + 1:nx],
        xout = lorenz63data[,(nparams + nx) + 1:nx],
        ts = fullTimeSeries,
        itime = itime,
        nparams = nparams,
        nx = nx,
        nstep = nstep)
    ## clean up
    rm(lorenz63data,fullTimeSeries)
    gc()

    ## generate test data for the Lorenz-96 model
    params <- c(8.0)
    nparams <- length(params)
    nx <- as.integer(15)
    nstep <- as.integer(5)
    nstepInit <- as.integer(50)
    dt = 0.1
    lorenz96data <- array(dim = c(nrealisations, nparams + nx*2))
    itime <- 1:nstep+1
    fullTimeSeries <- array(dim = c(nrealisations,nx,nstep+1))
    ##
    set.seed(42)
    ## set up the initial conditions
    xin <- rnorm(nx)
    F <- rgamma(nparams,32,4)
    ##
    res <- .Fortran('model_lorenz96',
                    x = xin,
                    xout = double(nx*(nstepInit+1)),
                    n = nx,
                    dt = dt,
                    nstep = nstepInit,
                    F = F,
                    PACKAGE = 'EmulatorAnalysis')
    xin <- res$x
    ## 
    for(i in 1:nrealisations){
        F <- rgamma(nparams,32,4)
        ##
        res <- .Fortran('model_lorenz96',
                        x = xin,
                        xout = double(nx*(nstep+1)),
                        n = nx,
                        dt = dt,
                        nstep = nstep,
                        F = F,
                        PACKAGE = 'EmulatorAnalysis')
        xout <- res$x
        ##
        if(any(!is.finite(xout))) stop('NaNs...')
        lorenz96data[i,] <- c(F,xin,xout)
        fullTimeSeries[i,,] <- matrix(res$xout,nx,nstep+1)
        ##
        res <- .Fortran('model_lorenz96',
                        x = xout,
                        xout = double(nx*(nstepInit+1)),
                        n = nx,
                        dt = dt,
                        nstep = nstepInit,
                        F = F,
                        PACKAGE = 'EmulatorAnalysis')
        xin <- res$x
    }
    if(any(!is.finite(lorenz96data))) stop('Ouch!')
    deterministicModelData[['lorenz96']] <- list(
        params = lorenz96data[,1:nparams,drop=FALSE],
        xin = lorenz96data[,nparams + 1:nx],
        xout = lorenz96data[,(nparams + nx) + 1:nx],
        ts = fullTimeSeries,
        itime = itime,
        nparams = nparams,
        nx = nx,
        nstep = nstep)
    ## clean up
    rm(lorenz96data,fullTimeSeries)
    gc()

    ## generate test data for the bank of seemingly stochastic ODEs
    ##
    ## set up the functional form of the derivatives of the ODEs
    gradient <- function(t,y,params, addSine = FALSE){
        if(all(dim(params) == length(y))){
            dy <- params %*% y
        } else if(all(dim(params) == (length(y) + c(0,1)))){
            dy <- params %*% c(1,y)
        }
        if(addSine){
            nn <- length(dy)
            for(i in 1:nn){
                dy[i] <- dy[i] + 0.3*sin((t+i/nn)*pi)*exp(-t/5) + 0.2*cos((t-i/nn)*pi)*exp(-t/5)
            }
        }
        return(list(dy))
    }
    ##
    jacobian <- function(t,y,params, addSine = FALSE){
        if(all(dim(params) == length(y))){
            ddy <- params
        } else if(all(dim(params) == (length(y) + c(0,1)))){
            ddy <- params[,-1]
        }
        return(ddy)
    }
    ##
    ntimes <- 201
    timeMax <- 2
    n <- 4 ## size of the system
    times <- seq(0,timeMax,length.out = ntimes)
    ensemble <- array(dim = c(nrealisations,ntimes,n))
    isample <- 1
    yinis <- array(dim = c(nrealisations,n))
    Ms <- array(dim = c(nrealisations,n^2))
    set.seed(42)
    for(isample in 1:nrealisations){
        yini <- rnorm(n)
        ##
        vecs <- matrix(rnorm(n^2),n,n)
        vals <- -1*sort(runif(n),decreasing=T)
        M <- vecs %*% diag(vals) %*% solve(vecs)
        ##
        out3  <- deSolve::lsode(yini, times, gradient, parms = M, jactype = "fullusr", jacfunc = jacobian, addSine = TRUE)
        ensemble[isample,,] <- out3[,-1]
        yinis[isample,] <- yini
        Ms[isample,] <- c(M)
    }
    if(any(!is.finite(ensemble))) stop('Ouch!')
    itime <- seq(1,ntimes,by = 10)
    deterministicModelData[['ODEperiodic']] <- list(
        params = Ms,
        xin = yinis,
        xout = ensemble[,ntimes,],
        ts = ensemble[,itime,],
        itime = itime,
        nparams = length(M),
        nx = n,
        nstep = ntimes)
    ## clean up
    rm(Ms,yinis,ensemble)
    gc()

    ## }}}


    ## {{{ run the emulators and the tests

    deterministicModels <- names(deterministicModelData)

    if(is.null(Model)){
        Model <- deterministicModels
    }

    nsets <- 3
    ftrain <- 1/nsets
    ntrain <- round(ftrain * nrealisations)
    set.seed(42)
    iorder <- sample(nrealisations)
    stopifnot(all(Model %in% deterministicModels))

    for(model in Model){
        if(model == 'lorenz96'){
            ixs <- head(round(seq(1,deterministicModelData[[model]]$nx,length.out = 5)),-1)
        } else {
            ixs <- 1:deterministicModelData[[model]]$nx
        }

        if(is.null(Ix)){
            Ix <- ixs
        }
        stopifnot(all(Ix %in% ixs))
        for(ix in Ix){
            cat(sprintf('Model = %s, ix = %3i',model,ix),fill=TRUE)
            nrealisations <- dim(deterministicModelData[[model]]$params)[1]
            isets <- list()
            for(iset in 1:nsets){
                inds <- iorder[(1 + ntrain*(iset-1)):ifelse(iset == nsets,nrealisations,ntrain*iset)]
                isets[[iset]] <- inds
            }
            ##

            df <- cbind(deterministicModelData[[model]]$xout[,ix],
                        deterministicModelData[[model]]$xin,
                        deterministicModelData[[model]]$params)
            colnames(df) <- c('y',
                              paste('x',1:deterministicModelData[[model]]$nx,sep=''), 
                              paste('p',1:deterministicModelData[[model]]$nparams,sep=''))
            df <- as.data.frame(df)

            emulatorNames <- c('lm','gam','svm','nnet','xgb','locfit', 'gp')
            nEmulators <- length(emulatorNames)
            
            ##
            out <- applyEmulators(Df = df,
                                  itrain = isets[[1]],
                                  itest = isets[[3]],
                                  predictorVariables = colnames(df)[-1],
                                  responseVariables = rep('y',nEmulators),
                                  returnSeFit = TRUE,
                                  tuneModels = TRUE,
                                  verbose = TRUE)

            predEns <- out$pred
            emulatorNames <- colnames(predEns)
            tunedParamsBase <- out$tunedParams
            nEmulators <- length(emulatorNames)
            predSe <- out$predSe
            rm(out)
            
            ## assess the skill of the individual emulators
            CORfull <- sapply(c('pearson','spearman'),function(m)cor(predEns, df[isets[[3]],1],method = m)[,1])
            RMSEfull <- apply(predEns - df[isets[[3]],1],2,sd)
            MADfull <- apply(predEns - df[isets[[3]],1],2,mad)

            ## calculate the standard deviation across models
            SDmodels <- apply(predEns,1,sd,na.rm=TRUE)
            
            ntrain1 <- 20 ## number of sub-samples
            ftrain2 <- 2/3/ntrain1 ## fraction of the whole data represented by the individual subsample
            ntrain2 <- round(ftrain2 * nrealisations) ## number of obs per subsample
            predEnsSmallLearners <- array(data = NA,
                                          dim = c(length(isets[[3]]),nEmulators,ntrain1),
                                          dimnames = list(NULL,emulatorNames,NULL))

            i2 <- 0
            for(iset in 1:ntrain1){
                cat('set',iset,'of',ntrain1,fill=TRUE)
                i1 <- i2+1
                if(iset == ntrain1){
                    i2 <- 2*ntrain
                } else {
                    i2 <- i1+ntrain2
                }
                out <- applyEmulators(Df = df,
                                      itrain = iorder[i1:i2],
                                      itest = isets[[3]],
                                      returnSeFit = FALSE,
                                      tuneModels = TRUE,
                                      tunedParams = if(iset == 1) NULL else tunedParamsSubset,
                                      verbose = (iset <= 3))
                tunedParamsSubset <- out$tunedParams
                predEnsSmallLearners[,,iset] <- out$pred
                rm(out)
            }

            ## see if the mean of the sub-sampled learners do any better:
            predEnsSmallLearnersMean <- apply(predEnsSmallLearners,1:2,mean)
            CORsub <- sapply(c('pearson','spearman'),function(m)cor(predEnsSmallLearnersMean, df[isets[[3]],1],method = m)[,1])
            RMSEsub <- apply(predEnsSmallLearnersMean - df[isets[[3]],1],2,sd)
            predEnsSmallLearnersSd <- apply(predEnsSmallLearners,1:2,sd)
            cat('cbind(CORfull[,"pearson"], CORsub[,"pearson"])\n')
            print(cbind(CORfull[,"pearson"], CORsub[,"pearson"]))
            cat('cbind(CORfull[,"spearman"], CORsub[,"spearman"])\n')
            print(cbind(CORfull[,"spearman"], CORsub[,"spearman"]))
            cat('cbind(RMSEfull,RMSEsub)\n')
            print(cbind(RMSEfull,RMSEsub))

            ## Next idea: apply the model to set 2, get errors, train error and test on part 3
            out <- applyEmulators(Df = df,
                                  itrain = isets[[1]],
                                  itest = 1:nrow(df), ## 'test' on the _whole_ dataset (for convenience, mostly - out-of-sample separation handled later
                                  returnSeFit = FALSE,
                                  predictorVariables = colnames(df)[-1],
                                  responseVariables = rep('y',nEmulators),
                                  tuneModels = FALSE,
                                  tunedParams = tunedParamsBase,
                                  verbose = TRUE)

            ## pred1 <- out$pred
            CHECK <- FALSE
            if(CHECK){
                ## check the fit - errors highly heteroscedastic
                plot(df$y,out$pred[,'xgb'])
                ## check that we get similar performance levels as above, and indeed we do
                cor(df$y,out$pred[,'xgb']) ## including the training sample: 0.9247479 better than before
                cor(df$y[c(isets[[2]],isets[[3]])],out$pred[c(isets[[2]],isets[[3]]),'xgb']) ## including the training sample: 0.8898612 same as before

            }
            
            errs <- (out$pred - df$y)
            colnames(errs) <- paste('e',colnames(out$pred),sep='_')
            dfAbsErrs <- cbind(as.data.frame(abs(errs)), df[,-1])
            out <- applyEmulators(Df = dfAbsErrs,
                                  itrain = isets[[2]], ## the original model was trained on set 1, train this model on set 2
                                  itest = isets[[3]], ## test on set 3
                                  predictorVariables = colnames(df)[-1],
                                  responseVariables = colnames(errs),
                                  returnSeFit = FALSE,
                                  tuneModels = TRUE,
                                  verbose = TRUE)


            predAbsErrs <- out$pred
            rm(out)

            ## set up quantiles (1% each)
            Qs <- seq(0,1,0.01)
            Qs2 <- seq(0,1,0.02)
            Qs5 <- seq(0,1,0.05)
            Qs6 <- seq(0,1,1/15)
            Qs10 <- seq(0,1,0.10)
            ## quantile mid-points
            QsMid <- (head(Qs,-1) + tail(Qs,-1))/2
            Qs2Mid <- (head(Qs2,-1) + tail(Qs2,-1))/2
            Qs5Mid <- (head(Qs5,-1) + tail(Qs5,-1))/2
            Qs6Mid <- (head(Qs6,-1) + tail(Qs6,-1))/2
            Qs10Mid <- (head(Qs10,-1) + tail(Qs10,-1))/2
            ##
            rmse <- function(x,y) sqrt(mean((x-y)^2))
            mae <- function(x,y) median(abs(x-y))

            
            pdf(file.path(plotsFolder,sprintf('errPlots_%s_ix%i.pdf',model,ix)))
            for(iEm in 1:nEmulators){
                ## only plot when there are data for this emulator
                if(all(!is.finite(predAbsErrs[,iEm]))) next
                ## avoid a problem of a 'null' emulator (affects nnet in some cases)
                quants <- quantile(predAbsErrs[,iEm], Qs)
                ## if(length(unique(quants)) != length(quants)) next
                ## 
                if(any(is.finite(predSe[,iEm]))){
                    nSources <- 3
                } else {
                    nSources <- 2
                }
                mfrow <- c(nSources,2)
                par(mfrow = mfrow,mar = c(0,0,0,0),oma = c(3.5,3.2,1.8,2))
                iplot <- 1
                for(iround in 1:2){
                    for(iSource in 1:nSources){
                        if(iSource == 1){
                            ## use jitter to deal with ties
                            absErrHat <- abs(jitter(predAbsErrs[,iEm],factor=1e-4))
                            xlabs <- c('Percentile of predicted errors','Mean predicted error in percentile')
                        } else if(iSource == 2) {
                            absErrHat <- abs(jitter(predEnsSmallLearnersSd[,iEm],factor=1e-4))
                            ## xlabs <- c('Percentile of ensemble standard deviation','Mean ensemble standard deviation within each percentile')
                        } else if(iSource == 3) {
                            absErrHat <- abs(jitter(predSe[,iEm],factor=1e-4))
                        }
                        ## break into quantiles
                        quants     <- quantile(absErrHat, Qs)
                        quants2    <- quantile(absErrHat, Qs2)
                        quants5    <- quantile(absErrHat, Qs5)
                        quants6    <- quantile(absErrHat, Qs6)
                        quants10    <- quantile(absErrHat, Qs10)
                        quantsMid  <- quantile(absErrHat, QsMid)
                        quants2Mid <- quantile(absErrHat, Qs2Mid)
                        quants5Mid <- quantile(absErrHat, Qs5Mid)
                        quants6Mid <- quantile(absErrHat, Qs6Mid)
                        quants10Mid <- quantile(absErrHat, Qs10Mid)
                        ## deal with any ties in the quantiles vector
                        if(length(unique(quants)) != length(quants)){
                            idup <- which(duplicated(quants))
                            idupdiff <- diff(idup)
                            if(any(idupdiff) == 1){
                                for(i in idup){
                                    quants[i] <- (quants[i-1] + quants[i+1])/2
                                }
                            } else {
                                j <- 1
                                while(j <= length(idup)){
                                    if(j == length(idup)){
                                        k <- j
                                    } else {
                                        for(k in j:length(idup)){
                                            if(! all(idup[j:k] == (idup[j] + seq_len(k-j+1) - 1))){
                                                k <- k-1
                                                break
                                            }
                                        }
                                    }
                                    quants[idup[j:k]] <- head(tail(seq(quants[idup[j]-1],quants[idup[k]+1],length.out = k-j+3),-1),-1)
                                    j <- k+1
                                }
                            }
                        }
                        quantsMid <- (head(quants,-1) + tail(quants,-1))/2
                        suppressWarnings({
                            ## ilastMid <- tail(which(abs(log(quantsMid/spline(Qs5Mid,quants5Mid,xout = QsMid)$y)) < log(1.2)),1)
                            ## ilastMid <- tail(which(abs(log(quantsMid/predict(smooth.spline(Qs5Mid,quants5Mid),QsMid)$y)) < log(1.2)),1)
                            ilastMid <- tail(which(abs(log(quantsMid/spline(Qs10Mid,quants10Mid,xout = QsMid)$y)) < log(1.2)),1)
                            ilab  <- cut(absErrHat, quants,  labels = FALSE, include.lowest=TRUE)
                            ilab2 <- cut(absErrHat, quants2, labels = FALSE, include.lowest=TRUE)
                            ilab5 <- cut(absErrHat, quants5, labels = FALSE, include.lowest=TRUE)
                            sdPerQuant  <- aggregate(errs[isets[[3]],iEm],by = list(quant=ilab ),FUN = function(x) mae(x,0))
                            sdPerQuant2 <- aggregate(errs[isets[[3]],iEm],by = list(quant=ilab2),FUN = function(x) mae(x,0))
                            sdPerQuant5 <- aggregate(errs[isets[[3]],iEm],by = list(quant=ilab5),FUN = function(x) mae(x,0))
                            ilast <- tail(which(abs(log(sdPerQuant$x/spline(Qs5Mid,sdPerQuant5$x,xout = QsMid)$y)) < log(2)),1)
                        })
                        if(iround == 1){
                            ## print(c(sdPerQuant$x[ilast],quantsMid[ilastMid]))
                            if(iSource == 1){
                                ## lim <- c(min(c(sdPerQuant$x,quantsMid)), max(c(sdPerQuant$x,quantsMid)))
                                lim <- c(min(c(sdPerQuant$x,quantsMid)), sdPerQuant$x[ilast])
                            } else {
                                lim <- c(min(c(sdPerQuant$x,quantsMid,lim[1])), max(c(sdPerQuant$x[ilast],quantsMid[ilastMid],lim[2])))
                            }
                        } else {
                            ## lim <- range(pretty(lim))
                            if(zapsmall(lim)[1] == 0) lim[1] <- lim[2]*1e-4
                            xticks <- pretty(c(0,100))[2:5]
                            yticks <- head(tail(pretty(lim),-1),-1)
                            ##
                            plot(QsMid*100,sdPerQuant$x,
                                 xlab = xlabs[1],
                                 ylab = 'MAE within each percentile',
                                 ylim = lim, yaxs = 'i',
                                 xaxt = 'n',yaxt = 'n',
                                 xaxs = 'i', xlim = c(0,100), 
                                 mgp = c(1.75,0.7,0))
                            fracText(0.95,0.85,LETTERS[iplot],cex=1.5)
                            abline(h = yticks, lty = 3, col = 'grey')
                            abline(v = xticks, lty = 3, col = 'grey')
                            ##
                            spl <- NULL
                            try({
                                spl <- smooth.spline(QsMid*100,sdPerQuant$x,cv = TRUE)
                            }, silent = TRUE)
                            if(is.null(spl)){
                                spl <- smooth.spline(QsMid*100,sdPerQuant$x,cv = FALSE)
                            }
                            ##
                            lines(spl,col = 3, lwd = 2)
                            cat(sprintf('model = %s, ix = %i, em = %s, iplot = %i,',model,ix,emulatorNames[iEm],iplot),sprintf('RSD = %.3g',sd(resid(spl))),fill = TRUE)
                            legend('topleft',
                                   sprintf('RSD = %.3g',sd(resid(spl))),
                                   lwd = NA, col = NA,lty = NA, bty = 'n',
                                   text.col = 3)
                            iside <- ((iplot - 1) %% 2)*2 + 2
                            axis(side = iside,mgp = c(1.75,0.7,0))
                            if(iplot >= (prod(mfrow)-1)){
                                axis(side=1,mgp = c(1.75,0.7,0))
                                mtext(xlabs[1],side=1,line=2)
                            }
                            iplot <- iplot + 1
                            ##
                            plot(quantsMid,sdPerQuant$x,
                                 xlab = xlabs[2],
                                 ylim = lim, yaxs = 'i',
                                 xlim = lim, xaxs = 'i',
                                 xaxt = 'n',yaxt = 'n',
                                 ylab = 'MAE within each quantile',
                                 mgp = c(1.75,0.7,0))
                            fracText(0.95,0.85,LETTERS[iplot],cex=1.5)
                            iside <- ((iplot - 1) %% 2)*2 + 2
                            axis(side = iside)
                            abline(h = yticks, lty = 3, col = 'grey')
                            abline(v = yticks, lty = 3, col = 'grey')
                            ## mtext(sprintf('Model: %s, parameter: %i, emulator: %s',model, ix, emulatorNames[iEm]),
                            ##       line = 0, cex = 1.4)
                            spl <- NULL
                            try({
                                spl <- smooth.spline(quantsMid,sdPerQuant$x,cv = TRUE)
                            }, silent = TRUE)
                            if(is.null(spl)){
                                spl <- smooth.spline(quantsMid,sdPerQuant$x,cv = FALSE)
                            }
                            lines(spl,col = 3, lwd = 2)
                            abline(0,1,lwd=2)
                            sdPerQuant$quantsMid <- quantsMid
                            LM <- MASS::rlm(x ~ quantsMid,data = sdPerQuant)
                            CF <- coef(LM)
                            abline(CF[1],CF[2],lwd=2,lty=2,col=2)
                            cat(sprintf('model = %s, ix = %i, em = %s, iplot = %i,',model,ix,emulatorNames[iEm],iplot),sprintf('y = %.3g + %.3g x',CF[1],CF[2]),fill = TRUE)
                            legend('topleft',
                                   c('y = x', sprintf('y = %.3g + %.3g x',CF[1],CF[2]),'Smoothing spline'),
                                   bg = 'white',
                                   lwd = 2,
                                   col = 1:3,lty = c(1,2,1))
                            cat(sprintf('model = %s, ix = %i, em = %s, iplot = %i,',model,ix,emulatorNames[iEm],iplot),
                                sprintf('R[p] = %.3g',cor(quantsMid,sdPerQuant$x)),fill = TRUE)
                            cat(sprintf('model = %s, ix = %i, em = %s, iplot = %i,',model,ix,emulatorNames[iEm],iplot),
                                sprintf('R[s] = %.3g',cor(quantsMid,sdPerQuant$x,method = 'spearman')),fill = TRUE)
                            cat(sprintf('model = %s, ix = %i, em = %s, iplot = %i,',model,ix,emulatorNames[iEm],iplot),
                                sprintf('RSD = %.3g',sd(resid(LM))),fill = TRUE)
                            legend('bottomright',
                                   c(parse(text = sprintf('R[p] == %.3g',cor(quantsMid,sdPerQuant$x))),
                                     parse(text = sprintf('R[s] == %.3g',cor(quantsMid,sdPerQuant$x,method = 'spearman'))),
                                     parse(text = sprintf('RSD == %.3g',sd(resid(LM))))),
                                   lwd = NA, col = NA,lty = NA, bty = 'n', text.col = c(1,1,2))
                            ##
                            if(iplot >= (prod(mfrow)-1)){
                                axis(side=1,mgp = c(1.75,0.7,0))
                                mtext(xlabs[2],side=1,line=2)
                            }
                            iplot <- iplot + 1
                        }
                    }
                }
                mtext(sprintf('Model: %s, parameter: %i, emulator: %s',model, ix, emulatorNames[iEm]),
                      line = 0.2, cex = 1.4,outer=TRUE)
                mtext('MAE within each percentile',side=2,outer=TRUE,line=2)
            }
            dev.off()

            if(model == 'lorenz63' && ix == 1){
                pdf(file.path(plotsFolder,'ErrorDistributionWrtParams-Lorenz63-XGB-ix1.pdf'))
                ## here we show that the parameters have a very strong impact on the spread of the errors
                par(mfrow = c(3,2), oma = c(3,3,0.5,3), mar = c(0,0,0,0))
                for(inputIndex in 1:with(deterministicModelData[['lorenz63']],nparams + nx)){
                    xvals <- errs[isets[[3]],'e_xgb']
                    if(inputIndex <= deterministicModelData[['lorenz63']]$nparams){
                        yvals <- deterministicModelData[['lorenz63']]$params[isets[[3]],inputIndex]
                    } else {
                        yvals <- deterministicModelData[['lorenz63']]$xin[isets[[3]],inputIndex - deterministicModelData[['lorenz63']]$nparams]
                    }
                    xlim <- quantile(abs(xvals),0.999) * c(-1,1)
                    ylim <- quantile(yvals,c(0.001,0.999))
                    smoothScatter(xvals, yvals,
                                  xlim = xlim, ylim = ylim, xlab = '',
                                  ylab = '', xaxt = 'n', yaxt = 'n')
                    yaxisSide <- ((inputIndex-1) %% 2)*2 + 2
                    axis(side = yaxisSide,mgp = c(1.75,0.7,0))
                    if(inputIndex >= 5) axis(side=1,mgp = c(1.75,0.7,0))
                    if(inputIndex <= deterministicModelData[['lorenz63']]$nparams){
                        mtext(sprintf('Parameter %i',inputIndex),line=1.7,side=yaxisSide,cex = 1)
                    } else {
                        mtext(sprintf('Variable %i',inputIndex - deterministicModelData[['lorenz63']]$nparams),line=1.7,side=yaxisSide,cex = 1)
                    }
                }
                mtext('Errors from xgb',outer = TRUE,line=2,side=1,cex = 1)
                dev.off()

                ## Does distance around the 'centre' help?
                nparams <- deterministicModelData[['lorenz63']]$nparams
                nx <- deterministicModelData[['lorenz63']]$nx
                lorenz63data <- with(deterministicModelData[['lorenz63']],cbind(xin,params))
                Colmeans <- colMeans(lorenz63data[,1:(nx + nparams)])
                Cov <- cov(lorenz63data[,1:(nx + nparams)])
                Sds <- apply(lorenz63data[,1:(nx + nparams)],2,sd)
                d <- t(lorenz63data[isets[[3]],1:(nx + nparams)]) - colMeans(lorenz63data[,1:(nx + nparams)])
                ## put together a bunch of different distance metrics
                distances <- data.frame(
                    ## the Mahalonobis distance
                    mdist = mahalanobis(lorenz63data[isets[[3]],1:(nx + nparams)], center = Colmeans, cov = Cov),
                    ## first without scaling the input distances
                    L1 = apply(d,2,function(x)sum(abs(x))),
                    L2 = apply(d,2,function(x)sqrt(sum(x^2))),
                    Linf = apply(d,2,function(x)max(abs(x))),
                    ## then scaling the inputs by the corresponding standard deviations
                    L1n = apply(d,2,function(x)sum(abs(x/Sds))),
                    L2n = apply(d,2,function(x)sqrt(sum((x/Sds)^2))),
                    Linfn = apply(d,2,function(x)max(abs(x/Sds))))
                rm(d)

                ## Mahalonobis distance provides the greatest discrimination, but way less powerful
                ## than a model for the errors
                sapply(distances,function(d)cor(d,abs(errs[isets[[3]],'e_xgb'])))
                
                pdf(file.path(plotsFolder,'ErrorDistributionWrtDistanceMetrics-Lorenz63-XGB-ix1.pdf'))
                par(mfrow = c(2,2), oma = c(1,3,0.5,2), mar = c(1.7,0,0,0))
                iplot <- 1
                for(idist in c(1,5,6,7)){
                    ##
                    quants <- quantile(distances[,idist], Qs)
                    quantsMid <- (head(quants,-1) + tail(quants,-1))/2
                    ilab <- cut(distances[,idist], quants, labels = FALSE, include.lowest=TRUE)
                    ## find the uncertainty within each quantile (measured as the RMSE within each quantile)
                    sdPerQuant <- aggregate(errs[isets[[3]],iEm],by = list(quant=ilab),FUN = function(x) mae(x,0))
                    ##
                    plot(quantsMid,sdPerQuant$x,
                         xlab = '',
                         ## xlab = xlabs[2],
                         ## ylim = lim, yaxs = 'i',
                         ## xlim = lim, xaxs = 'i',
                         yaxt = 'n',
                         ylab = 'MAE within each quantile',
                         mgp = c(1.75,0.7,0))
                    fracText(0.95,0.85,LETTERS[iplot],cex=1.5)
                    iside <- ((iplot - 1) %% 2)*2 + 2
                    axis(side = iside)
                    sdPerQuant$quantsMid <- quantsMid
                    LM <- lm(x ~ quantsMid,data = sdPerQuant)
                    CF <- coef(LM)
                    abline(CF[1],CF[2],lwd=2,lty=2,col=2)
                    legend('bottomright',
                           c(sprintf('R = %.3g',cor(quantsMid,sdPerQuant$x)),
                             sprintf('RSD = %.3g',sd(resid(LM)))),
                           lwd = NA, col = NA,lty = NA, bty = 'n', text.col = c(1,2))
                    ##
                    if(iplot >= (prod(mfrow)-1)){
                        axis(side=1,mgp = c(1.75,0.7,0))
                        mtext(xlabs[2],side=1,line=2)
                    }
                    iplot <- iplot + 1
                }
                mtext('Distance metric', side = 1, line = 0,outer=TRUE)
                mtext('MAE within each percentile',side=2,outer=TRUE,line=2)
                dev.off()

            }
        }
    }
    ## }}}

    if(!is.null(outFile)){
        sink()
    }

    return(invisible(NULL))
}
