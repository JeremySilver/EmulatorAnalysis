#' Apply the emulators to a dataset (e.g. as produced by a deterministic simulation model)
#'
#' @param Df A \code{data.frame} object to apply the emulators to
#' @param itrain Indices of the training set (i.e. row numbers of \code{Df}).
#' @param itest Indices of the test set (i.e. row numbers of \code{Df}).
#' @param responseVariables A string, or vector or strings, giving the name(s) of the response variable(s) to be used in the emulators. These should appear as variable names in the \code{Df} argument. If one string is provided, all emulators will use the same response variable; if this argument is a vector of strings, it should have the same length as the \code{emulatorNames} argument (below). The default value is \code{'y'}.
#' @param predictorVariables A vector of strings giving the names of the predictor variables to be used by the emulators. These should appear as variable names in the \code{Df} argument. The default value is \code{NULL}, which results in all variable names from \code{Df} other than those specified by \code{responseVariables}.
#' @param returnSeFit Calculate and return the `native' estimate of uncertainty from the emulators, where available (i.e., the standard error of the fit, the width of the confidence interval, or the posterior standard deviation of the estimate). This adds considerably to the time and memory requirements for the Gaussian process regression, in particular. The default value is \code{FALSE}.
#' @param tuneModels Boolean: should hyperparameter estimation be undertaken? The default value is \code{TRUE}. If \code{TRUE}, then the resulting hyperparameters are returned.
#' @param tunedParams If the hyperparameter estimation has already been undertaken for this dataset, this can be provided as a list. The default value is \code{NULL}, in which case models are
#' @param emulatorNames A vector of strings giving the names of the emulators, defaults to \code{c('lm','gam','svm','nnet','xgb','locfit', 'gp').}
#' @param verbose Boolean: should extra output be printed to the R console? Defaults to \code{FALSE}.

#' @return A \code{list} object with the following 2 or 3 entries (there will be 3 entries if argument \code{returnSeFit} is \code{TRUE}, and 2 otherwise). \code{pred}: a \code{data.frame} object with column names given by \code{emulatorNames}, and the rows corresponding to \code{itest}. \code{tunedParams}: the list of tuned hyperparameters. \code{predSe}: the native estimate of uncertainty from each emulator for which it is available (this return value is only included if \code{returnSeFit} is \code{TRUE}); this has the same format as the return value of \code{pred} (a \code{data.frame} object with column names given by \code{emulatorNames}, and the rows corresponding to \code{itest}).
#'
#' @importFrom dials grid_latin_hypercube tree_depth learn_rate
#' @importFrom e1071 svm tune.control
#' @importFrom kernlab gausspr sigest
#' @importFrom locfit locfit
#' @importFrom mgcv gam predict.gam
#' @importFrom parsnip boost_tree set_engine set_mode
#' @importFrom rsample vfold_cv
#' @importFrom nnet nnet nnet.formula
#' @importFrom stats lm predict.lm
#' @importFrom tune control_grid select_best tune_grid
#' @importFrom workflows add_formula add_model workflow
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis legend lines mtext par plot smoothScatter text
#' @importFrom stats aggregate as.formula coef cor cov mad mahalanobis median optim quantile resid rgamma rnorm runif sd smooth.spline spline
#' @importFrom magrittr %>%
#'
#' @export

applyEmulators <- function(Df,
                           itrain,
                           itest,
                           responseVariables = 'y',
                           predictorVariables = NULL,
                           returnSeFit = FALSE,
                           tuneModels = TRUE,
                           tunedParams = NULL,
                           emulatorNames = c('lm','gam','svm','nnet','xgb','locfit', 'gp'),
                           verbose = FALSE){
    requireNamespace('xgboost')  ## xgb.train
    requireNamespace('mgcv')     ## gam
    requireNamespace('e1071')    ## svm
    requireNamespace('nnet')     ## nnet
    requireNamespace('locfit')   ## locfit
    requireNamespace('kernlab')  ## gausspr
    ## requireNamespace('tidymodels') ## for tuning xgboost
    requireNamespace('parsnip')
    requireNamespace('dials')
    requireNamespace('rsample')
    requireNamespace('workflows')
    requireNamespace('tune')
    requireNamespace('stats')
    requireNamespace('graphics')
    requireNamespace('grDevices')

    nEmulators <- length(emulatorNames)
    
    ##
    predEns <- matrix(NA,length(itest),nEmulators,
                      dimnames = list(NULL,emulatorNames))
    if(returnSeFit){
        predSe <- matrix(NA,length(itest),nEmulators,
                         dimnames = list(NULL,emulatorNames))
    }
    ##
    stopifnot(all(responseVariables %in% names(Df)))
    if(length(responseVariables) == 1){
        responseVariables <- rep(responseVariables,nEmulators)
    } else if(length(responseVariables) != nEmulators){
        stop('length of responseVariables vector does not match the number of emulators')
    }
    ##
    if(is.null(predictorVariables)){
        predictorVariables <- setdiff(names(Df),responseVariables)
    }
    if(tuneModels){
        tunedParams = list()
    } else {
        stopifnot(class(tunedParams) == 'list')
        if(length(predictorVariables) <= 6){
            stopifnot(all(sort(names(tunedParams)) == sort(c('svm','nnet','xgb','locfit','gp'))))
        } else {
            stopifnot(all(sort(names(tunedParams)) == sort(c('svm','nnet','xgb','gp'))))
        }
    }
    
    ## LM
    iEm <- 1
    res <- try({
        tm <- system.time({
            Formula <- as.formula(sprintf('%s~%s',responseVariables[iEm],paste(predictorVariables,collapse='+')))
            if(verbose) print(Formula,showEnv = FALSE)
            pred <- predict.lm(lm(Formula, data = Df[itrain,]),
                                      newdata = Df[itest,predictorVariables],
                                      se.fit = returnSeFit)
            if(returnSeFit){
                predEns[,iEm] <- pred$fit
                predSe[,iEm] <- pred$se.fit
            } else {
                predEns[,iEm] <- pred
            }
        })
        if(verbose) cat(sprintf('%s: %g seconds\n',emulatorNames[iEm],tm[3]))
        1
    })
    if(class(res) != 'numeric') cat('Issue with',emulatorNames[iEm],fill=TRUE)
    ## GAM
    iEm <- iEm + 1
    res <- try({
        tm <- system.time({
            Formula <- as.formula(paste(responseVariables[iEm],'~',paste(paste('s(',predictorVariables,')'),collapse='+')))
            if(verbose) print(Formula,showEnv = FALSE)
            pred <- predict.gam(gam(Formula, data = Df[itrain,]),
                            newdata = Df[itest,predictorVariables],
                            se.fit = returnSeFit)
            if(returnSeFit){
                predEns[,iEm] <- pred$fit
                predSe[,iEm] <- pred$se.fit
            } else {
                predEns[,iEm] <- pred
            }
            ## rm(pred)
        })
        if(verbose) cat(sprintf('%s: %g seconds\n',emulatorNames[iEm],tm[3]))
        1
    })
    if(class(res) != 'numeric') cat('Issue with',emulatorNames[iEm],fill=TRUE)
    ## SVM
    iEm <- iEm + 1
    res <- try({
        Formula <- as.formula(sprintf('%s~%s',responseVariables[iEm],paste(predictorVariables,collapse='+')))
        if(verbose) print(Formula,showEnv = FALSE)
        if(tuneModels){
            if(verbose) cat(sprintf('tuning %s...\n',emulatorNames[iEm]))
            tm <- system.time({
                FN <- function(par){
                    set.seed(42)
                    if(min(par) < 1e-6){
                        return(Inf)
                    } else {                   
                        suppressWarnings({
                            obj <- svm(Formula, data = Df[itrain,], 
                                       gamma = max(par[1],1e-6),
                                       cost = max(par[2],1e-6),
                                       epsilon = max(par[3],1e-6),
                                       cross =3)
                        })
                        return(mean(obj$MSE))
                    }
                }
                ##
                par0 <- c(gamma = 1/length(predictorVariables), cost = 1, epsilon = .1)
                ##
                out <- optim(par = par0, fn = FN,
                             method = "Nelder-Mead", 
                             control = list(maxit = 100, reltol = 5e-3))
                parms <- out$par
                names(parms) <- names(par0)
                tunedParams[[emulatorNames[iEm]]] <- parms
            })
            if(verbose) cat(sprintf('tune %s: %g seconds\n',emulatorNames[iEm],tm[3]))
        } else {
            parms <- tunedParams[[emulatorNames[iEm]]]
        }
        tm <- system.time({
            predEns[,iEm] <- predict(svm(Formula, data = Df[itrain,],
                                         gamma = parms['gamma'], cost = parms['cost'], epsilon = parms['epsilon'],
                                         fitted = FALSE),
                                     newdata = Df[itest,predictorVariables])
        })
        if(verbose) cat(sprintf('run %s: %g seconds\n',emulatorNames[iEm],tm[3]))
        1
    })
    if(class(res) != 'numeric') cat('Issue with',emulatorNames[iEm],fill=TRUE)
    ## n-net
    iEm <- iEm + 1
    res <- try({
        ## create a new variable for the formula
        yRescaled <- sprintf('%s_rescaled',responseVariables[iEm])
        Formula <- as.formula(sprintf('%s ~ %s',yRescaled,paste(predictorVariables,collapse='+')))
        ## transform the data to a 0-1 scale
        ymin <- min(Df[,responseVariables[iEm]])
        yrange <- diff(range(Df[,responseVariables[iEm]]))
        Df[,yRescaled] <- (Df[,responseVariables[iEm]] - ymin)/yrange
        ##
        if(verbose) print(Formula,showEnv = FALSE)
        if(tuneModels){
            if(verbose) cat(sprintf('tuning %s...\n',emulatorNames[iEm]))
            tm <- system.time({
                obj <- e1071::tune(nnet,
                                   Formula, data = Df[itrain,], 
                                   ranges = list(size = unique(round(seq(1/3*length(predictorVariables), 3*length(predictorVariables),length.out = 7))),
                                                 decay = c(0,5e-4,5e-3)),
                                   maxit = 200,
                                   trace = FALSE,
                                   MaxNWts = 2000,
                                   tunecontrol = tune.control(sampling = "cross", cross =3))
                parms <- unlist(obj$best.parameters)
                tunedParams[[emulatorNames[iEm]]] <- parms
            })
            if(verbose) cat(sprintf('tune %s: %g seconds\n',emulatorNames[iEm],tm[3]))
        } else {
            parms <- tunedParams[[emulatorNames[iEm]]]
        }
        ##
        tm <- system.time({
            predEns[,iEm] <- predict(nnet.formula(formula = Formula,
                                                  data = Df[itrain,],
                                                  size = parms['size'],
                                                  trace = FALSE,
                                                  decay = parms['decay'],
                                                  MaxNWts = 2000,
                                                  maxit = 200),
                                     newdata = Df[itest,predictorVariables ])
            ## transform back to the original scale
            predEns[,iEm] <- predEns[,iEm]*yrange + ymin
        })
        if(verbose) cat(sprintf('run %s: %g seconds\n',emulatorNames[iEm],tm[3]))
        Df[,yRescaled] <- NULL
        1
    })
    if(class(res) != 'numeric') cat('Issue with',emulatorNames[iEm],fill=TRUE)
    ## xgboost with black box regressors - 6 seconds
    iEm <- iEm + 1
    if(tuneModels){
        if(verbose) cat(sprintf('tuning %s...\n',emulatorNames[iEm]))
        tm <- system.time({

            set.seed(123)

            ## following steps suggested in https://juliasilge.com/blog/xgboost-tune-volleyball/

            dtrain <- xgb.DMatrix(data = as.matrix(Df[itrain,predictorVariables]), label = Df[itrain,responseVariables[iEm]])
            xgb_spec <- boost_tree(
                trees = 100, 
                tree_depth = tune::tune(), ## maps to 'max_depth'
                learn_rate = tune::tune(), ##  maps to 'eta'. step size
                ) %>% 
                set_engine("xgboost") %>% 
                set_mode("regression") 

                xgb_grid <- grid_latin_hypercube(
                    tree_depth(range = c(3,12)),
                    learn_rate(range = c(-3,-0.5)),
                    size = 30)

                Formula <- as.formula(sprintf('%s~%s',responseVariables[iEm],paste(predictorVariables,collapse='+')))
                xgb_wf <- workflow() %>%
                    add_formula(Formula) %>%
                    add_model(xgb_spec)
                
                set.seed(123)
                folds <- vfold_cv(Df[itrain,c(predictorVariables,responseVariables[iEm])], v = 3)

                set.seed(234)
                xgb_res <- tune_grid(
                    xgb_wf,
                    resamples = folds,
                    grid = xgb_grid,
                    control = control_grid(save_pred = TRUE))
                
                ## show_best(xgb_res, "rmse")
                best_rmse <- select_best(xgb_res, "rmse")
                parms <- unlist(best_rmse[,1:(ncol(best_rmse)-1)])
                
                tunedParams[[emulatorNames[iEm]]] <- parms
        })
        if(verbose) cat(sprintf('tune %s: %g seconds\n',emulatorNames[iEm],tm[3]))
    } else {
        parms <- tunedParams[[emulatorNames[iEm]]]
    }
    ##
    xgbParams <- list(max_depth = parms['tree_depth'],
                      eta = parms['learn_rate'], 
                      nthread = 2,
                      objective = "reg:squarederror")
    
    res <- try({
        tm <- system.time({
            ## set up the xgb dataset
            dtrain <- xgb.DMatrix(data = as.matrix(Df[itrain,predictorVariables]), label = Df[itrain,responseVariables[iEm]])
            dtest <- xgb.DMatrix(data = as.matrix(Df[itest,predictorVariables]))
            predEns[,iEm] <- predict(xgb.train(params = xgbParams, data = dtrain, nrounds = 100, print_every_n = 10, verbose = 10), newdata = dtest)
            
            ## prevent chronic under-dispersion
            iter <- 1
            while(mad(predEns[,iEm]) < mad(Df[itrain,responseVariables[iEm]])*0.2 && iter <= 20){
                xgbParams$eta <- min(xgbParams$eta*5,0.9)
                cat('\tRerunning with increased learning rate:',xgbParams$eta,fill=TRUE)
                predEns[,iEm] <- predict(xgb.train(params = xgbParams, data = dtrain, nrounds = 100, print_every_n = 10, verbose = 10), newdata = dtest)
                iter <- iter + 1
            }

            
        })
        if(verbose) cat(sprintf('run %s: %g seconds\n',emulatorNames[iEm],tm[3]))
        1
    })
    if(class(res) != 'numeric') cat('Issue with',emulatorNames[iEm],fill=TRUE)
    ## local regression model - 14.890 seconds
    iEm <- iEm + 1
    if(length(predictorVariables) <= 6){
        if(tuneModels){
            if(verbose) cat(sprintf('tuning %s...\n',emulatorNames[iEm]))
            tm <- system.time({
                FN <- function(par){
                    set.seed(42)
                    Formula <- as.formula(sprintf('%s~lp(%s,nn=%g)',responseVariables[iEm],paste(predictorVariables,collapse=','),min(max(par[1],0.01),0.99)))
                    res <- try({
                        obj <- e1071::tune(locfit, Formula, data = Df[itrain,],
                                           maxk = 200,
                                           tunecontrol = tune.control(sampling = "cross", cross =3))
                    }, silent = TRUE)
                    if(class(res) == 'try-error'){
                        obj <- list(best.performance = Inf)
                    }
                    return(obj$best.performance)
                }
                ##
                par0 <- c(nn = 0.7)
                ##
                out <- optim(par = par0, fn = FN,
                             method = "Brent", lower = c(0.1),upper = c(0.99),
                             control = list(maxit = 100, reltol = 1e-2))
                parms <- out$par
                names(parms) <- names(par0)
                tunedParams[[emulatorNames[iEm]]] <- parms
            })
            if(verbose) cat(sprintf('tune %s: %g seconds\n',emulatorNames[iEm],tm[3]))
        } else {
            parms <- tunedParams[[emulatorNames[iEm]]]
        }
        res <- try({
            tm <- system.time({
                Formula <- as.formula(sprintf('%s~lp(%s,nn=%g)',responseVariables[iEm],paste(predictorVariables,collapse=','),parms[1]))
                if(verbose) print(Formula,showEnv = FALSE)
                pred <- predict(locfit(Formula, data = Df[itrain,], maxk = 200),
                                newdata = Df[itest,predictorVariables],
                                se.fit = returnSeFit)
                if(returnSeFit){
                    predEns[,iEm] <- pred$fit
                    predSe[,iEm] <- pred$se.fit
                } else {
                    predEns[,iEm] <- pred
                }
                rm(pred)                    
            })
            if(verbose) cat(sprintf('run %s: %g seconds\n',emulatorNames[iEm],tm[3]))
            1
        })
        if(class(res) != 'numeric') cat('Issue with',emulatorNames[iEm],fill=TRUE)
    }

    ## Gaussian process
    iEm <- iEm + 1
    res <- try({
        if(tuneModels){
            if(verbose) cat(sprintf('tuning %s...\n',emulatorNames[iEm]))
            tm <- system.time({
                FN <- function(par){
                    set.seed(42)
                    if(min(par) < 1e-6){
                        return(Inf)
                    } else {                   
                        gp <- gausspr(x = as.matrix(Df[itrain,predictorVariables]),
                                      y = as.matrix(Df[itrain,responseVariables[iEm]]),
                                      type= 'regression',
                                      kernel="rbfdot",
                                      kpar = list(sigma = max(par[1],1e-3)),
                                      var = max(par[2],1e-3),
                                      cross = 3,
                                      variance.model = FALSE)
                        return(gp@cross)
                    }
                }
                ##
                par0 <- c(sigma = mean(sigest(x = as.matrix(Df[itrain,predictorVariables]), scaled = FALSE)[c(1, 3)]), var = 1)
                ##
                out <- optim(par = par0, fn = FN,
                             method = "Nelder-Mead", 
                             control = list(maxit = 100, reltol = 1e-2))
                parms <- out$par
                names(parms) <- names(par0)
                tunedParams[[emulatorNames[iEm]]] <- parms
            })
            if(verbose) cat(sprintf('tune %s: %g seconds\n',emulatorNames[iEm],tm[3]))
        } else {
            parms <- tunedParams[[emulatorNames[iEm]]]
        }
        ##
        tm <- system.time({
            gp <- gausspr(x = as.matrix(Df[itrain,predictorVariables]),
                          y = as.matrix(Df[itrain,responseVariables[iEm]]),
                          type= 'regression',
                          kernel="rbfdot",
                          kpar = list(sigma = parms['sigma']),
                          var = parms['var'],
                          variance.model = returnSeFit)
            predEns[,iEm] <- predict(gp,newdata = Df[itest,predictorVariables], type = "response")
        })
        if(verbose) cat(sprintf('%s: %g seconds\n',emulatorNames[iEm],tm[3]))
        1
    })
    if(class(res) != 'numeric') cat('Issue with',emulatorNames[iEm],fill=TRUE)

    out <- list(pred = predEns,
                tunedParams = tunedParams)
    if(returnSeFit){
        tm <- system.time({
            predSe[,iEm] <- predict(gp,newdata = Df[itest,predictorVariables], type="sdeviation")
        })
        out$predSe <- predSe
        if(verbose) cat(sprintf('predict,gausspr,sdeviation: %g seconds\n',tm[3]))
    }
    return(out)
}
