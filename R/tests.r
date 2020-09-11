runEpiExamples <- FALSE
if(runEpiExamples){
    library(EpiModel)

    ## Example 1: SI Model (One-Group)
                                        # Set parameters
    param <- param.dcm(inf.prob = 0.2, act.rate = 0.25)
    init <- init.dcm(s.num = 500, i.num = 1)
    control <- control.dcm(type = "SI", nsteps = 500)
    mod1 <- dcm(param, init, control)
    pdf('mod1_dcm.pdf')
    plot(mod1)
    dev.off()
    
    ## Example 2: SIR Model with Vital Dynamics (One-Group)
    param <- param.dcm(inf.prob = 0.2, act.rate = 5,
                       rec.rate = 1/3, a.rate = 1/90, ds.rate = 1/100,
                       di.rate = 1/35, dr.rate = 1/100)
    init <- init.dcm(s.num = 500, i.num = 1, r.num = 0)
    control <- control.dcm(type = "SIR", nsteps = 500)
    mod2 <- dcm(param, init, control)
    pdf('mod2_dcm.pdf')
    plot(mod2)
    dev.off()
    
    ## Example 3: SIS Model with act.rate Sensitivity Parameter
    param <- param.dcm(inf.prob = 0.2, act.rate = seq(0.1, 0.5, 0.1),
                       rec.rate = 1/50)
    init <- init.dcm(s.num = 500, i.num = 1)
    control <- control.dcm(type = "SIS", nsteps = 500)
    mod3 <- dcm(param, init, control)
    pdf('mod3_dcm.pdf')
    plot(mod3)
    dev.off()
    
    ## Example 4: SI Model with Vital Dynamics (Two-Group)
    param <- param.dcm(inf.prob = 0.4,  inf.prob.g2 = 0.1,
                       act.rate = 0.25, balance = "g1",
                       a.rate = 1/100, a.rate.g2 = NA,
                       ds.rate = 1/100, ds.rate.g2 = 1/100,
                       di.rate = 1/50, di.rate.g2 = 1/50)
    init <- init.dcm(s.num = 500, i.num = 1,
                     s.num.g2 = 500, i.num.g2 = 0)
    control <- control.dcm(type = "SI", nsteps = 500)
    mod4 <- dcm(param, init, control)
    pdf('mod4_dcm.pdf')
    plot(mod4)
    dev.off()

    
    ## Example 1: SI Model
    param <- param.icm(inf.prob = 0.2, act.rate = 0.25)
    init <- init.icm(s.num = 500, i.num = 1)
    control <- control.icm(type = "SI", nsteps = 500, nsims = 10)
    mod1 <- icm(param, init, control)
    pdf('mod1_icm.pdf')
    plot(mod1)
    dev.off()
    
    ## Example 2: SIR Model
    param <- param.icm(inf.prob = 0.2, act.rate = 0.25, rec.rate = 1/50)
    init <- init.icm(s.num = 500, i.num = 1, r.num = 0)
    control <- control.icm(type = "SIR", nsteps = 500, nsims = 10)
    mod2 <- icm(param, init, control)
    pdf('mod2_icm.pdf')
    plot(mod2)
    dev.off()
    
    ## Example 3: SIS Model
    param <- param.icm(inf.prob = 0.2, act.rate = 0.25, rec.rate = 1/50)
    init <- init.icm(s.num = 500, i.num = 1)
    control <- control.icm(type = "SIS", nsteps = 500, nsims = 10)
    mod3 <- icm(param, init, control)
    pdf('mod3_icm.pdf')
    plot(mod3)
    dev.off()
    
    ## Example 4: SI Model with Vital Dynamics (Two-Group)
    param <- param.icm(inf.prob = 0.4,  inf.prob.g2 = 0.1,
                       act.rate = 0.25, balance = "g1",
                       a.rate = 1/100, a.rate.g2 = NA,
                       ds.rate = 1/100, ds.rate.g2 = 1/100,
                       di.rate = 1/50, di.rate.g2 = 1/50)
    init <- init.icm(s.num = 500, i.num = 1,
                     s.num.g2 = 500, i.num.g2 = 0)
    control <- control.icm(type = "SI", nsteps = 500, nsims = 10)
    mod4 <- icm(param, init, control)
    pdf('mod4_icm.pdf')
    plot(mod4)
    dev.off()


    ## Example 1: Independent SI Model
    ## Network model estimation
    nw <- network.initialize(n = 100, bipartite = 50, directed = FALSE)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    est1 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15)
    init <- init.net(i.num = 10, i.num.m2 = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, verbose.int = 0)
    mod1 <- netsim(est1, param, init, control)
    
    ## Print, plot, and summarize the results
    pdf('mod1_net.pdf')
    plot(mod1)
    dev.off()
    summary(mod1, at = 50)
    
    ## Example 2: Dependent SIR Model
    ## Recalculate dissolution coefficient with departure rate
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20,
                                   d.rate = 0.0021)
    
    ## Reestimate the model with new coefficient
    est2 <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    ## Reset parameters to include demographic rates
    param <- param.net(inf.prob = 0.3, inf.prob.m2 = 0.15,
                       rec.rate = 0.02, rec.rate.m2 = 0.02,
                       a.rate = 0.002, a.rate.m2 = NA,
                       ds.rate = 0.001, ds.rate.m2 = 0.001,
                       di.rate = 0.001, di.rate.m2 = 0.001,
                       dr.rate = 0.001, dr.rate.m2 = 0.001)
    init <- init.net(i.num = 10, i.num.m2 = 10,
                     r.num = 0, r.num.m2 = 0)
    control <- control.net(type = "SIR", nsteps = 100, nsims = 5)
    
    ## Simulate the model with new network fit
    mod2 <- netsim(est2, param, init, control)
    
    ## Print, plot, and summarize the results
    pdf('mod2_net.pdf')
    plot(mod2)
    dev.off()
    summary(mod2, at = 100)
    
}

runLorenzTlAdTests <- FALSE
if(runLorenzTlAdTests){
    pwd <- getwd()
    ## setwd('~/projects/DSI_emulators/EmulatorProject/src')
    ## system('R CMD SHLIB -o lorenz.so *.f90')
    ## dyn.load('lorenz.so')

    xinit = sin(1:3)
    x1 = xinit
    nstep <- as.integer(100)
    params <- c(1.3,0.6,2.0)
    dt = 0.1

    {
        x1 = xinit
        print(x1, digits=17)
        x1 <- .Fortran('model_lorenz63',
                       x = x1,
                       dt = dt,
                       nstep = nstep,
                       params = params)$x
        print(x1, digits=17)
        ##
        x1 = xinit
        print(x1, digits=17)
        x1 <- .Fortran('model_lorenz63',
                       x = x1,
                       dt = dt,
                       nstep = nstep,
                       params = params)$x
        print(x1, digits=17)
    }
    

    ## call model(x1,dt,n)
    dxinit = sin(1:3 + 123)
    dx = dxinit
    for(m in 0:11){
        x2 = xinit + dx
        x3 = xinit
        x_tl = dx
        ## cat('dx  '); print(dx,digits=17)
        ## cat('x2.0 '); print(x2,digits=17)
        x2 <- .Fortran('model_lorenz63',
                       x = x2,
                       dt = dt,
                       nstep = nstep,
                       params = params)$x
        ## cat('x2.1 '); print(x2,digits=17)
        ## cat('',fill=TRUE)
        out <- .Fortran('model_tl_lorenz63',
                       x = x3,
                       x_tl = x_tl,
                       dt = dt,
                       nstep = nstep,
                       params = params)
        x3 <- out$x
        x_tl <- out$x_tl
        print(c(m, 1 - mean(abs((x2 - x1)/x_tl))),digits=15)
        dx = dx*0.1
    }
    print(x2)
    print(x3)

    dyinit = log(cos(sin(1:3))**2)

    x = xinit
    x_tl = dxinit
    x_tl <- .Fortran('model_tl_lorenz63',
                    x = x,
                    x_tl = x_tl,
                    dt = dt,
                    nstep = nstep,
                    params = params)$x_tl 
    ##
    x = xinit
    x_ad = x_tl
    x_ad <- .Fortran('model_ad_lorenz63',
                     x = x,
                     x_ad = x_ad,
                     dt = dt,
                     nstep = nstep,
                     params = params)$x_ad
    tl_disp = sum( x_tl * x_tl )
    ad_disp = sum( dxinit * x_ad )
    cat('test adj 1',tl_disp, ad_disp, 1.0 - tl_disp / ad_disp, fill = TRUE)

    x = xinit
    x_tl = dxinit
    x_tl <- .Fortran('model_tl_lorenz63',
                    x = x,
                    x_tl = x_tl,
                    dt = dt,
                    nstep = nstep,
                    params = params)$x_tl 
    x = xinit
    x_ad = dyinit
    x_ad <- .Fortran('model_ad_lorenz63',
                     x = x,
                     x_ad = x_ad,
                     dt = dt,
                     nstep = nstep,
                     params = params)$x_ad
    cat( 'test adj 2',sum(x_tl*dyinit), sum(dxinit*x_ad), log(sum(x_tl*dyinit)/sum(dxinit*x_ad)), fill = TRUE)

    ##

    n = as.integer(40)
    dt = 0.1
    nstep <- as.integer(100)
    F = 8.0
  
    xinit = sin(1:n)
    x1 = xinit
    print(x1, digits=17)
    x1 <- .Fortran('model_lorenz96',
                   x = x1,
                   n = n,
                   dt = dt,
                   nstep = nstep,
                   F = F)$x
    print(x1, digits=17)

    dxinit = sin(1:n + 123)
    dx = dxinit
    for(m in 0:15){
        x2 = xinit + dx
        x3 = xinit
        x_tl = dx
        x2 <- .Fortran('model_lorenz96',
                       x = x2,
                       n = n,
                       dt = dt,
                       nstep = nstep,
                       F = F)$x
        out <- .Fortran('model_tl_lorenz96',
                        x = x3,
                        x_tl = x_tl ,
                        n = n,
                        dt = dt,
                        nstep = nstep,
                        F = F)
        x3 <- out$x
        x_tl <- out$x_tl
        cat(m,x2[1],x_tl[1],  1 - sum(abs((x2 - x1)/x_tl))/n,fill=TRUE)
        dx = dx*0.1
    }

    dyinit = log(cos(sin(1:n))**2)

    ## <M \delta x, M \delta x> = <\delta x, M* M \delta x>
    x = xinit
    x_tl = dxinit
    x_tl <- .Fortran('model_tl_lorenz96',
                     x = x,
                     x_tl = x_tl ,
                     n = n,
                     dt = dt,
                     nstep = nstep,
                     F = F)$x_tl
    x = xinit
    x_ad = x_tl
    x_ad <- .Fortran('model_ad_lorenz96',
                     x = x,
                     x_ad = x_ad ,
                     n = n,
                     dt = dt,
                     nstep = nstep,
                     F = F)$x_ad
    cat('test adj 1',sum(x_tl*x_tl), sum(dxinit*x_ad), 1.0 - sum(x_tl*x_tl)/sum(dxinit*x_ad),fill=TRUE)

    ## <M \delta x, \delta y> = <\delta x, M* \delta y>
    x = xinit
    x_tl = dxinit
    x_tl <- .Fortran('model_tl_lorenz96',
                     x = x,
                     x_tl = x_tl ,
                     n = n,
                     dt = dt,
                     nstep = nstep,
                     F = F)$x_tl
    x = xinit
    x_ad = dyinit
    x_ad <- .Fortran('model_ad_lorenz96',
                     x = x,
                     x_ad = x_ad ,
                     n = n,
                     dt = dt,
                     nstep = nstep,
                     F = F)$x_ad
    cat('test adj 2',sum(x_tl*dyinit), sum(dxinit*x_ad), 1.0 - sum(x_tl*dyinit)/sum(dxinit*x_ad),fill = TRUE)
    
    
}    

    
generateLorenzData <- FALSE
if(generateLorenzData){
    params <- c(10, 28, 8/3)
    params <- c(1.3,0.6,2.0)
    nparams <- length(params)
    nx <- 3
    xinit = sin(1:nx)
    dt = 0.1
    ntime <- 500
    nstep <- as.integer(ntime/dt)
    ## 1e6 realisations takes ~8 seconds
    nrealisations <- 1e5
    lorenz63data <- matrix(nrow = nrealisations,
                           ncol = nparams + nx*2)

    set.seed(42)
    for(i in 1:nrealisations){
        xin <- rnorm(nx,0,2)
        xin <- runif(nx,0.5,4)
        ## xin <- c(1,1,1)
        ## params <- rgamma(nparams,2,2)
        params <- rgamma(nparams,2,2)
        xout <- .Fortran('model_lorenz63',
                         x = xin,
                         dt = dt,
                         nstep = nstep,
                         params = params)$x
        lorenz63data[i,] <- c(params,xin,xout)
    }
    if(any(!is.finite(lorenz63data))) stop('Ouch!')

    ftrain <- 0.6
    ntrain <- round(ftrain * nrealisations)
    ix <- 2
    set.seed(42)
    ## for(ix in 1:nx)
    {
        itrain <- sample(nrealisations, ntrain)
        dtrain <- xgb.DMatrix(lorenz63data[itrain,1:(nx + nparams)], label = lorenz63data[itrain,(nx+nparams) + ix])
        dtest <- xgb.DMatrix(lorenz63data[-itrain,1:(nx+nparams)], label = lorenz63data[-itrain,(nx+nparams) + ix])
        watchlist <- list(train = dtrain, eval = dtest)
        ## eta = 0.5:  eval-rmse:3.669825 
        ## eta = 0.05: eval-rmse:3.650622 
        param <- list(max_depth = 6, eta = 0.1, verbose = 0, nthread = 2,
                      objective = "reg:squarederror")
        bst <- xgb.train(params = param, data = dtrain, nrounds = 100, watchlist = watchlist, print_every_n = 10)
    }
    sd(lorenz63data[-itrain,(nx+nparams)+ix])
    pred <- predict(bst, newdata = dtest)
    cor(pred,lorenz63data[-itrain,(nx+nparams) + ix])
    smoothScatter(x = pred,y = lorenz63data[-itrain,(nx+nparams) + ix])
    abline(0,1)
    sd(pred) / sd(lorenz63data[-itrain,(nx+nparams) + ix]) ## ratio of SDs
    mean(pred) - mean(lorenz63data[-itrain,(nx+nparams) + ix]) ## bias
    
    ##########################
    ##

    library(EmulatorProject)
    library(xgboost)
    params <- c(8.0)
    nparams <- length(params)
    nx <- as.integer(40)
    nstep <- as.integer(20)
    dt = 0.1
    ## this takes ~8 seconds
    nrealisations <- 100000
    lorenz96data <- matrix(nrow = nrealisations,
                           ncol = nparams + nx*2)
    ##
    set.seed(42)
    for(i in 1:nrealisations){
        xin <- rnorm(nx)
        F <- rgamma(nparams,32,4)
        ##
        xout <- .Fortran('model_lorenz96',
                   x = xin,
                   n = nx,
                   dt = dt,
                   nstep = nstep,
                   F = F)$x
        ##
        if(any(!is.finite(xout))) stop('NaNs...')
        lorenz96data[i,] <- c(F,xin,xout)
    }

    ftrain <- 0.6
    ntrain <- ftrain * nrealisations
    ix <- 20
    set.seed(42)
    ## for(ix in 1:nx)
    {
        itrain <- sample(nrealisations, ntrain)
        dtrain <- xgb.DMatrix(lorenz96data[itrain,1:(nx+1)], label = lorenz96data[itrain,(nx+1) + ix])
        dtest <- xgb.DMatrix(lorenz96data[-itrain,1:(nx+1)], label = lorenz96data[-itrain,(nx+1) + ix])
        watchlist <- list(train = dtrain, eval = dtest)
        ## eta = 0.5:  eval-rmse:3.669825 
        ## eta = 0.05: eval-rmse:3.650622 
        param <- list(max_depth = 6, eta = 0.1, verbose = 0, nthread = 2,
                      objective = "reg:squarederror")
        bst <- xgb.train(params = param, data = dtrain, nrounds = 100, watchlist = watchlist, print_every_n = 10)
    }
    sd(lorenz96data[-itrain,1:(nx+1)])
    pred <- predict(bst, newdata = dtest)
    cor(pred,lorenz96data[-itrain,(nx+1) + ix])
    smoothScatter(x = pred,y = lorenz96data[-itrain,(nx+1) + ix])
    abline(0,1)
    sd(pred) / sd(lorenz96data[-itrain,(nx+1) + ix])
    mean(pred) - mean(lorenz96data[-itrain,(nx+1) + ix])
    
}

