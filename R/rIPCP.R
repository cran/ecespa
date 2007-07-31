`rIPCP` <-
function (x, lambda=NULL, type=1,lmax = NULL, win = owin(c(0, 1), c(0, 1)), ...) {
    x.ppp <- x$data
    xw <- x.ppp$w
    rho <- x$rho 
    sigma <- sqrt(x$sigma2)
    mu <- x.ppp$n/(rho*area.owin(x.ppp$w))

    if(is.null(lambda)) lambda <- x$lambda
    
    win <- if (is.im(lambda)) rescue.rectangle(as.owin(lambda)) else as.owin(win)
    
    x.ppp$window <- win
    
    if (is.null(lmax)) {
        imag <- as.im(lambda, win, ...)
        summ <- summary(imag)
        #lmax <- summ$max + 0.05 * diff(summ$range)
        lmax <- summ$max 
    }
    
    if (is.im(lambda)) {

        # mean retention probability of the thinning process (se podr<ed>a meter alg<fa>n corrector "+..." para afinar mejor)

          probm <- mean(eval.im(lambda/lmax))


        # "prethinnig" required n to obtain the original n after thinning
          n.preth <- x.ppp$n / probm

        # "prethinning" required mu and rho to obtain the original mu and rho after thinning
          mu.preth <- n.preth*mu/x.ppp$n #ojo, directamente mu.preth = mu/probm
          rho.preth <- rho/probm

        #pre-thining point pattern:
 
        if(type==1)  X <- rThomas (win=x.ppp$w, kappa=rho, sigma=sigma, mu.preth) #simula resultados exactos al modelo ajustado
 
        if(type==2)  X <- rThomas (win=x.ppp$w, kappa=rho.preth, sigma=sigma, mu=mu ) #simula resultados m<e1>s dispersos"

        if (X$n == 0) 
            return(X)

        # post-thining point pattern:

        prob <- lambda[X]/lmax
        u <- runif(X$n)
        retain <- (u <= prob)
        X <- X[retain, ]
        ##X$window <- xw        
        return(X)
    }
    stop(paste(sQuote("lambda"), "must be a constant, a function or an image"))
}

