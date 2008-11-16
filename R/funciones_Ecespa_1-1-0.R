

###############################################################################
###############################################################################
################ FUNCIÓN dixon 2002 (original ecespa) #####################################
########## y funciones accesorias: check, ginv, mNNinfo2, mNNinfo, mNNtest y NNid ###################
###############################################################################


dixon2002 <- function (datos, nsim = 99) 
{
    info = mNNinfo(xy = datos[, 1:2], label = datos[, 3])
    datos.test = mNNtest(info)
    Ni = rowSums(info$ON)
    S = log10((info$ON/(Ni - info$ON))/(info$EN/(Ni - info$EN)))
    ON = info$ON
    EN = info$EN
    Z = datos.test$Z
    Z.obs = Z
    pZas = 2 * (ifelse(Z >= 0, 1 - pnorm(Z), pnorm(Z)))
    C = datos.test$C[1]
    C.obs = C
    Ci = datos.test$Ci[, 1]
    Ci.obs = Ci
    pCas = datos.test$C[2]
    pCias = datos.test$Ci[, 2]
    pZr = NULL
    pCr = NULL
    pCir = NULL
    if (nsim > 0) {
        for (i in 1:nsim) {
            print(i)
            datos[, 3] = sample(datos[, 3])
            info = mNNinfo(xy = datos[, 1:2], label = datos[, 
                3])
            datos.test = mNNtest(info)
            Z = cbind(Z, datos.test$Z)
            C = c(C, datos.test$C[1])
            Ci = cbind(Ci, datos.test$Ci[, 1])
        }
        pZr = apply(Z, 1, p2colasr)
        pCr = 1 - rank(C)[1]/(length(C))
        pCir = apply(Ci, 1, function(x) 1 - rank(x)[1]/(length(x)))
    }
    St = as.data.frame(as.table(round(S, 2)))
    ONt = as.data.frame(as.table(ON))
    ENt = as.data.frame(as.table(round(EN, 2)))
    Zt = round(datos.test$Z, 2)
    round(pZas, 4)
    tableZ = cbind(ONt[order(ONt[, 1]), ], ENt[order(ENt[, 1]), 
        3], St[order(St[, 1]), 3], round(Z.obs, 2), round(pZas, 
        4))
    names(tableZ) = c("From", "To", "    Obs.Count", "    Exp. Count", 
        "S ", "Z ", "  p-val.as")
    if (length(pZr) != 0) {
        tableZ = cbind(tableZ, round(pZr, 4))
        names(tableZ) = c(names(tableZ)[-8], "  p-val.rnd")
    }
    rownames(tableZ) = NULL
    k = length(unique(datos[, 3]))
    df = c(k * (k - 1), rep(k - 1, k))
    nombres.test = c("Overall segregation", paste("From ", dimnames(EN)[[1]], 
        "          "))
    tablaC = data.frame(cbind(df, round(c(C.obs, Ci.obs), 2), 
        round(c(pCas, pCias), 4)))
    row.names(tablaC) = nombres.test
    names(tablaC) = c("  df ", "Chi-sq", "P.asymp")
    if (length(pCir) != 0) {
        tablaC = cbind(tablaC, round(c(pCr, pCir), 4))
        names(tablaC) = c(names(tablaC)[-4], "  P.rand")
    }
    return(list(ON = ON, EN = EN, Z = Z.obs, S = S, pZas = pZas, 
        pZr = pZr, C = C.obs, Ci = Ci.obs, pCas = pCas, pCias = pCias, 
        pCr = pCr, pCir = pCir, tablaZ = tableZ, tablaC = tablaC))
}





check <- function (x, v, l1, l2) 
{
    if (v[l1, l2] > 0) {
        print("WARNING from routine ", x, ": element", l1, l2, 
            "already non-zero")
    }
    v[l1, l2] <- x
    v
}

ginv <- function (m) 
{
    temp <- eigen(m, symmetric = T)
    va <- temp$values
    ve <- temp$vectors
    va <- ifelse((abs(va) < 1e-09), 0, 1/va)
    va2 <- 0 * m
    diag(va2) <- va
    ve %*% va2 %*% t(ve)
}


mNNinfo2 <-function (n, R, Q) 
{
    N <- sum(n)
    k <- length(n)
    l <- names(n)
    EN <- matrix(0, nrow = k, ncol = k)
    VN <- VarN <- matrix(0, nrow = k * k, ncol = k * k)
    for (i in 1:k) {
        for (j in 1:k) {
            EN[i, j] <- n[i] * (n[j] - (i == j))/(N - 1)
        }
    }
    for (l1 in 1:(k * k)) {
        i <- 1 + (l1 - 1)%/%k
        j <- 1 + (l1 - 1)%%k
        for (l2 in l1:(k * k)) {
            i2 <- 1 + (l2 - 1)%/%k
            j2 <- 1 + (l2 - 1)%%k
            if ((i == i2) & (j == j2)) {
                if (i == j) {
                  p2 <- n[i] * (n[i] - 1)/(N * (N - 1))
                  p3 <- p2 * (n[i] - 2)/(N - 2)
                  p4 <- p3 * (n[i] - 3)/(N - 3)
                  VN <- check(1, VN, l1, l2)
                  VarN[l1, l2] <- (N + R) * p2 + (2 * N - 2 * 
                    R + Q) * p3 + (N * (N - 3) - Q + R) * p4 - 
                    EN[i, j] * EN[i, j]
                }
                else {
                  p2 <- n[i] * n[j]/(N * (N - 1))
                  p3 <- p2 * (n[i] - 1)/(N - 2)
                  p4 <- p3 * (n[j] - 1)/(N - 3)
                  VN <- check(2, VN, l1, l2)
                  VarN[l1, l2] <- N * p2 + Q * p3 + (N * (N - 
                    3) - Q + R) * p4 - EN[i, j] * EN[i, j]
                }
            }
            else if ((i == j) & (i == i2) & (j != j2)) {
                p3 <- n[i] * (n[i] - 1) * n[j2]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[i] - 2)/(N - 3)
                VN <- check(3, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i2 == j2) & (i == i2) & (j != j2)) {
                p3 <- n[i2] * (n[i2] - 1) * n[j]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[i] - 2)/(N - 3)
                VN <- check(3, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i2 == j2) & (j == j2) & (i != i2)) {
                p3 <- n[j] * (n[j] - 1) * n[i]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[j] - 2)/(N - 3)
                VN <- check(4, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R + Q) * 
                  p3 + (N * (N - 3) - Q + R) * p4 - EN[i, j] * 
                  EN[i2, j2]
            }
            else if ((i == j) & (i == j2) & (i != i2)) {
                p3 <- n[i] * (n[i] - 1) * n[i2]/(N * (N - 1) * 
                  (N - 2))
                p4 <- p3 * (n[i] - 2)/(N - 3)
                VN <- check(14, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R + Q) * 
                  p3 + (N * (N - 3) - Q + R) * p4 - EN[i, j] * 
                  EN[i2, j2]
            }
            else if ((i == j) & (i2 == j2) & (i != i2)) {
                p4 <- n[i] * (n[i] - 1) * n[i2] * (n[i2] - 1)/(N * 
                  (N - 1) * (N - 2) * (N - 3))
                VN <- check(5, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i == j) & (i2 != i) & (j2 != j) & (i2 != 
                j2)) {
                p4 <- n[i] * (n[i] - 1) * n[i2] * n[j2]/(N * 
                  (N - 1) * (N - 2) * (N - 3))
                VN <- check(6, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i2 == j2) & (i2 != i) & (j2 != j) & (i != 
                j)) {
                p4 <- n[i2] * (n[i2] - 1) * n[i] * n[j]/(N * 
                  (N - 1) * (N - 2) * (N - 3))
                VN <- check(6, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i == i2) & (i != j) & (i2 != j2) & (j != 
                j2)) {
                p4 <- n[i] * (n[i] - 1) * n[j] * n[j2]/(N * (N - 
                  1) * (N - 2) * (N - 3))
                VN <- check(7, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i == j2) & (i2 == j) & (i != j)) {
                p2 <- n[i] * n[j]/(N * (N - 1))
                p3 <- p2 * (n[i] - 1 + n[j] - 1)/(N - 2)
                p4 <- p2 * (n[i] - 1) * (n[j] - 1)/((N - 2) * 
                  (N - 3))
                VN <- check(8, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- R * p2 + (N - 
                  R) * p3 + (N * (N - 3) - Q + R) * p4 - EN[i, 
                  j] * EN[i2, j2]
            }
            else if ((i != j) & (j == i2) & (i2 != j2) & (i != 
                j2)) {
                p3 <- n[i] * n[j] * n[j2]/(N * (N - 1) * (N - 
                  2))
                p4 <- p3 * (n[j] - 1)/(N - 3)
                VN <- check(9, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i != j) & (j == j2) & (i2 != j2) & (i != 
                i2)) {
                p3 <- n[i] * n[j] * n[i2]/(N * (N - 1) * (N - 
                  2))
                p4 <- p3 * (n[j] - 1)/(N - 3)
                VN <- check(10, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- Q * p3 + (N * 
                  (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
            else if ((i != j) & (i == j2) & (i2 != j2) & (j != 
                i2)) {
                p3 <- n[i] * n[j] * n[i2]/(N * (N - 1) * (N - 
                  2))
                p4 <- p3 * (n[i] - 1)/(N - 3)
                VN <- check(11, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N - R) * p3 + 
                  (N * (N - 3) - Q + R) * p4 - EN[i, j] * EN[i2, 
                  j2]
            }
            else if ((i != j) & (i != i2) & (i != j2) & (j != 
                i2) & (j != j2) & (i2 != j2)) {
                p4 <- n[i] * n[j] * n[i2] * n[j2]/(N * (N - 1) * 
                  (N - 2) * (N - 3))
                VN <- check(12, VN, l1, l2)
                VarN[l1, l2] <- VarN[l2, l1] <- (N * (N - 3) - 
                  Q + R) * p4 - EN[i, j] * EN[i2, j2]
            }
        }
    }
    v <- as.vector(t(outer(l, l, paste, sep = "")))
    dimnames(VN) <- dimnames(VarN) <- list(v, v)
    list(EN = EN, VN = VN, VarN = VarN)
}


mNNinfo <- function (xy, label, nnid = NULL, splancs = TRUE) 
{
    if (is.null(nnid)) {
        nnid <- NNid(xy, splancs)
    }
    n <- table(label)
    N <- sum(n)
    l <- names(n)
    k <- length(l)
    R <- sum((1:N) == nnid[nnid])
    Q1 <- table(c(1:6, table(nnid)))
    Q <- 2 * Q1[2] + 6 * Q1[3] + 12 * Q1[4] + 20 * Q1[5] + 30 * 
        Q1[6] - 70
    ON <- matrix(0, nrow = k, ncol = k)
    for (i in 1:k) {
        for (j in 1:k) {
            ON[i, j] <- sum((label == l[i]) & (label[nnid] == 
                l[j]))
        }
    }
    temp <- mNNinfo2(n, R, Q)
    rownames(ON) <- colnames(ON) <- rownames(temp$EN) <- colnames(temp$EN) <- l
    list(ON = ON, EN = temp$EN, VarN = temp$VarN, R = R, Q = Q)
}


mNNtest <- function (info, obsN = NULL) 
{
    if (is.null(obsN)) 
        obsN <- as.vector(t(info$ON))
    expN <- as.vector(t(info$EN))
    varN <- diag(info$VarN)
    Z <- (obsN - expN)/sqrt(varN)
    names(Z) <- dimnames(info$VarN)[[1]]
    k <- nrow(info$EN)
    delta <- as.matrix(obsN - expN)
    C <- t(delta) %*% ginv(info$VarN) %*% delta
    pC <- 1 - pchisq(C, k * (k - 1))
    Ci <- rep(0, k)
    for (i in 1:k) {
        i1 <- 1 + (i - 1) * k
        i2 <- i1 + (k - 2)
        Ci[i] <- t(delta[i1:i2, ]) %*% solve(info$VarN[i1:i2, 
            i1:i2]) %*% delta[i1:i2]
    }
    pCi <- 1 - pchisq(Ci, k - 1)
    list(Z = Z, C = c(C, pC), Ci = cbind(Ci, pCi))
}


NNid <- function (xy, splancs = TRUE) 
{
    if (splancs) {
        nndistG(xy)$neighs
    }
    else {
        temp <- find.neighbor(xy, k = 2)
        as.vector(temp[temp[, 1] != temp[, 2], 2])
    }
}
















###############################################################################
###############################################################################
################ FUNCIÓN getis (original ecespa) ########################################
# CHANGES nueva versión de getis que devuelve objetos de clase "ecespa.getis"
#CHANGES: S· method.
# rd modificado
################################################################################

getis <- function (mippp, nx = 30, ny = 30, R = 10) 
{
    require(spatstat)
    dataname <- deparse (substitute(mippp))
    cosagrid <- gridcenters(mippp$window, nx = nx, ny = ny)
    cosa1 <- ppp(x = mippp$x[1], y = mippp$y[1], window = mippp$window, marks = "1")
    cosa2 <- setmarks(mippp, "2")
    cosa12 <- superimpose(cosa1, cosa2)
    coso12 <- cosa12
    cosagrid <- as.data.frame(cosagrid)
    cosagrid <- rbind(cbind(x=mippp$x,y=mippp$y),cosagrid)

     klocalgrid <- apply(cosagrid, 1, function(x, cosa12=coso12) {
        cosa12$x[1] <- x[1]
        cosa12$y[1] <- x[2]
        Kcross(unique(cosa12), i = "1", j = "2", r = c(0, R), correction = "isotropic")$iso[2]
    })


    result <- list(x = cosagrid$x, y = cosagrid$y, klocal = klocalgrid[(1:mippp$n)], 
                   klocalgrid=klocalgrid, R = R, nx=nx, ny=ny, dataname=dataname, ppp = mippp)
    class(result) <- c("ecespa.getis", class(result))
   return(result)


}


# nueva S3 method función que sustituye a la antigua "getis.plot"
#rd modificado

plot.ecespa.getis <- function(getis.obj, type="k", interp=100, color=tim.colors(64), 
                    contour=TRUE , points=TRUE,...){
   require(fields)
   require (akima)
   require (spatstat)
   lambda = getis.obj$ppp$n/area.owin(getis.obj$ppp$window)
   if (type=="k") zg = getis.obj$klocalgrid
   if (type== "l") zg = sqrt(getis.obj$klocalgrid/pi)
   if (type== "n") zg = getis.obj$klocalgrid*lambda
   if (type== "d") zg = sqrt(getis.obj$klocalgrid/pi)-getis.obj$R
   seqx= seq(round(getis.obj$ppp$window$xrange)[1],
             round(getis.obj$ppp$window$xrange)[2],
             length=interp)
   seqy= seq(round(getis.obj$ppp$window$yrange)[1],
             round(getis.obj$ppp$window$yrange)[2],
             length=interp)
   image.plot(interp(x=getis.obj$x, y=getis.obj$y, z= zg,
                     xo=seqx, yo=seqy), col=color,...)
   if(contour==TRUE) contour(interp(x=getis.obj$x, y=getis.obj$y, z= zg,
                                     xo=seqx, yo=seqy), add=TRUE)
   if(points==TRUE) plot(getis.obj$ppp, pch=16,, cex=0.7, add=T)
}

print.ecespa.getis <- function(getis.obj,...)
{

cat("Getis local density function of the dataset", getis.obj$dataname,  "\n computed for a radius R=",getis.obj$R, 
", with a grid of", getis.obj$nx, "x",getis.obj$ny,".\n")
     
cat("Plot it to see the result.\n")
}




###############################################################################
###############################################################################
######## FUNCIÓN haz.ppp (original ecespa) ###############################################
######## Sigue siendo útil ###########################################################

haz.ppp <- function (W) 
{
    require(spatstat)
    if (dim(W)[2] == 2) 
        pepe = ppp(x = W[, 1], y = W[, 2], xrange = range(W[, 
            1]), yrange = range(W[, 2]))
    if (dim(W)[2] == 3) 
        pepe = ppp(x = W[, 1], y = W[, 2], xrange = range(W[, 
            1]), yrange = range(W[, 2]), marks = W[, 3])
    return(pepe)
}


###############################################################################
###############################################################################
######## FUNCIÓN K012 (original ecespa) ###############################################
######## Mantener siempre. Es la que se usa en Ecography 2008################################
#  Modificado K012.rd: Título: test against independent labelling
#  Modificado K012.rd: References: Referencia correcta del artículo (y del de Ecosistemas (páginas))
# !!!!!!!!!!!!! Modificar función: que detecte que X es ppp con marcas discretas y >=3 levels



K012 <- function (X, fijo, i, j, nsim = 99, nrank = 1, r = NULL, correction = "isotropic") 
{
    marx <- marks(X)
    fijo <- (marx == fijo)
    I <- (marx == i)
    J <- (marx == j)
    cosa <- Kmulti.ls(X, fijo, I, r, corre = correction)
    r = cosa$r
    k01.obs <- cosa[[3]]
    k02.obs <- Kmulti.ls(X, fijo, J, r, corre = correction)[[3]]
    k01.sim <- NULL
    k02.sim <- NULL
    cat("Generating simulations...")
    for (n in 1:nsim) {
        progressreport(n, nsim)
        X$marks[I | J] <- sample(X$marks[I | J])
        marx <- marks(X)
        I <- (marx == i)
        J <- (marx == j)
        k01.sim <- cbind(k01.sim, Kmulti.ls(X, fijo, I, r, corre = correction)[[3]])
        k02.sim <- cbind(k02.sim, Kmulti.ls(X, fijo, J, r, corre = correction)[[3]])
    }
    orderstat <- function(x, n) sort(x)[n]
    k01.lo <- apply(k01.sim, 1, orderstat, n = nrank)
    k01.hi <- apply(k01.sim, 1, orderstat, n = nsim - nrank + 
        1)
    k02.lo <- apply(k02.sim, 1, orderstat, n = nrank)
    k02.hi <- apply(k02.sim, 1, orderstat, n = nsim - nrank + 
        1)
    k01 <- cosa
    k01[[2]] <- k01.hi
    k01[[3]] <- k01.obs
    k01[[4]] <- k01.lo
    attributes(k01)$labl <- c(attributes(cosa)$labl[1], "hi(r)", 
        attributes(cosa)$labl[3], "lo(r)")
    attributes(k01)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k01)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k01)$desc <- c(attributes(cosa)$desc[1], "upper pointwise envelope of simulations", 
        attributes(cosa)$desc[3], "lower pointwise envelope of simulations")
    k02 <- cosa
    k02[[2]] <- k02.hi
    k02[[3]] <- k02.obs
    k02[[4]] <- k02.lo
    attributes(k02)$labl <- c(attributes(cosa)$labl[1], "hi(r)", 
        attributes(cosa)$labl[3], "lo(r)")
    attributes(k02)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k02)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k02)$desc <- c(attributes(cosa)$desc[1], "upper pointwise envelope of simulations", 
        attributes(cosa)$desc[3], "lower pointwise envelope of simulations")
    return(list(k01 = k01, k02 = k02))
}



###############################################################################
###############################################################################
######## FUNCIÓN K1K2 (original ecespa) ###############################################
######## Mantener siempre. Es la que se usa en Ecography 2008################################
# Modificado K1K2.rd: References: Referencia correcta del artículo (y del de Ecosistemas (páginas))
# !!!!!!!!!!!!! Modificar función: que detecte que X es ppp con marcas discretas 


K1K2 <- function (X, i, j, nsim = 99, nrank = 1, r = NULL, correction = "isotropic") 
{
    marx <- marks(X)
    I <- (marx == i)
    J <- (marx == j)
    k12 = Kmulti.ls(X, I, J, r, corre = correction)
    r = k12$r
    k1 = Kest(split(X)[names(split(X)) == i][[1]], r = r, correction = correction)
    k2 = Kest(split(X)[names(split(X)) == j][[1]], r = r, correction = correction)
    k1k2.o = k1[[3]] - k2[[3]]
    k112 = k1[[3]] - k12[[3]]
    k212 = k2[[3]] - k12[[3]]
    k1k2.s = NULL
    k112.s = NULL
    k212.s = NULL
    for (n in 1:nsim) {
        progressreport(n, nsim)
        X$marks[I | J] = sample(X$marks[I | J])
        marx <- marks(X)
        I <- (marx == i)
        J <- (marx == j)
        k12.s = Kmulti.ls(X, I, J, r = r, corre = correction)
        k1.s = Kest(split(X)[names(split(X)) == i][[1]], r = r, 
            correction = correction)
        k2.s = Kest(split(X)[names(split(X)) == j][[1]], r = r, 
            correction = correction)
        k1k2.s = cbind(k1k2.s, k1.s[[3]] - k2.s[[3]])
        k112.s = cbind(k112.s, k1.s[[3]] - k12.s[[3]])
        k212.s = cbind(k212.s, k2.s[[3]] - k12.s[[3]])
    }
    orderstat <- function(x, n) sort(x)[n]
    k1k2.lo <- apply(k1k2.s, 1, orderstat, n = nrank)
    k1k2.hi <- apply(k1k2.s, 1, orderstat, n = nsim - nrank + 
        1)
    k112.lo <- apply(k112.s, 1, orderstat, n = nrank)
    k112.hi <- apply(k112.s, 1, orderstat, n = nsim - nrank + 
        1)
    k212.lo <- apply(k212.s, 1, orderstat, n = nrank)
    k212.hi <- apply(k212.s, 1, orderstat, n = nsim - nrank + 
        1)
    k1k2 = k12
    k1k2[[2]] = k1k2.hi
    k1k2[[3]] = k1k2.o
    k1k2[[4]] = k1k2.lo
    attributes(k1k2)$labl = c(attributes(k12)$labl[1], "hi(r)", 
        "K1(r)-K2(r)", "lo(r)")
    attributes(k1k2)$ylab = expression(K[1] - K[2])
    attributes(k1k2)$names = c(attributes(k12)$names[1], "hi", 
        "K1-K2", "lo")
    attributes(k1k2)$desc = c(attributes(k12)$desc[1], "upper pointwise envelope of simulations", 
        paste("differences of ", attributes(k12)$desc[3]), "lower pointwise envelope of simulations")
    k1k12 = k12
    k1k12[[2]] = k112.hi
    k1k12[[3]] = k112
    k1k12[[4]] = k112.lo
    attributes(k1k12)$labl = c(attributes(k12)$labl[1], "hi(r)", 
        "K1(r)-K12(r)", "lo(r)")
    attributes(k1k12)$ylab = expression(K[1] - K[12]^"*")
    attributes(k1k12)$names = c(attributes(k12)$names[1], "hi", 
        "K1-K12", "lo")
    attributes(k1k12)$desc = c(attributes(k12)$desc[1], "upper pointwise envelope of simulations", 
        "difference of K1 and K12 functions", "lower pointwise envelope of simulations")
    k2k12 = k12
    k2k12[[2]] = k212.hi
    k2k12[[3]] = k212
    k2k12[[4]] = k212.lo
    attributes(k2k12)$labl = c(attributes(k12)$labl[1], "hi(r)", 
        "K2(r)-K12(r)", "lo(r)")
    attributes(k2k12)$ylab = expression(K[2] - K[12]^"*")
    attributes(k2k12)$names = c(attributes(k12)$names[1], "hi", 
        "K2-K12", "lo")
    attributes(k2k12)$desc = c(attributes(k12)$desc[1], "upper pointwise envelope of simulations", 
        "difference of K2 and K12 functions", "lower pointwise envelope of simulations")
    return(list(k1k2 = k1k2, k1k12 = k1k12, k2k12 = k2k12))
}





###############################################################################
###############################################################################
######## FUNCIÓN Kmm (original ecespa) ###############################################
######## Mantener siempre. Se usa en libro AEET 2008################################
#13/06/08: intento de añadir envueltas por permutación de las marcas
#29/06/08: print.ecespa.kmm
#28/07/08: Modificado Kmm.rd: S3 print and plot methods y  class; nueva estructura del resultado)
# 28/07/08: Modificado Kmm.rd: References: Referencia del libro AEET)
# 28/07/08: Modificados ejemplos: plot
# 28/07/08: Modificada función: que permita calcular con y sin envueltas



Kmm <- function (mippp, r = 1:10, nsim=NULL) 
{
    require(spatstat)
    dataname <- deparse(substitute(mippp))
    mippp1 = setmarks(mippp, 1)
    lambda = mippp$n/area.owin(mippp$window)
    mu = mean(mippp$marks)
    Kmm = NULL
    Kmm1 = NULL
    Kmmsim <- NULL
 

  if(!is.null(nsim)){
         pepesim <- list(mippp)
         pepesim <- rep(pepesim, nsim) 
         pepesim <- lapply(pepesim, function(x)  x = x %mark% sample(x$marks))
         # class(pepesim) <- c(class(pepesim), "splitppp")
        
         m0sim <- sapply(pepesim, markstat, sum, R=0)
	 }
       m0 = markstat(mippp, sum, R = 0)
       m01 = markstat(mippp1, sum, R = 0)
       
    


    cat("computing Kmm for distance class ")
    for (i in 1:length(r)) {
    	progressreport(i, length(r))
        mr <- markstat(mippp, sum, R = r[i])
        mr1 <- markstat(mippp1, sum, R = r[i])
			
        sumatorio <- (mr - m0) * m0
        sumatorio1 <- (mr1 - m01) * m01
        	    
	E0 <- mean(sumatorio)
        E01 <- mean(sumatorio1)
		
        Kmm = c(Kmm, E0/(lambda * mu^2))
        Kmm1 = c(Kmm1, E01/(lambda))
	
	if(!is.null(nsim)){
	   mrsim <- sapply(pepesim,markstat,sum, R=r[i])
	   sumatoriosim <- (mrsim - m0sim) * m0sim
	   E0sim <- apply(sumatoriosim,2, mean)
	   Kmmsim <- rbind(Kmmsim,E0sim/(lambda*mu^2))
	}
		
   }
    
   if(!is.null(nsim)) Kmmsim.n <- apply(Kmmsim,2, function(x) x/Kmm1) else Kmmsim.n <-NULL
	   
	   result=list(dataname=dataname, r = r, nsim=nsim, kmm = Kmm, 
                   kmm.n = Kmm/Kmm1, kmmsim=Kmmsim, 
                   kmmsim.n=Kmmsim.n)

   class(result) =c("ecespa.kmm",class(result))
	return(result)   
}




plot.ecespa.kmm <- 
function(kmm.ob,type="Kmm.n", q=0.025, 
            xlime=NULL, ylime=NULL,  maine=NULL, add=F, kmean=TRUE,
            ylabe=NULL, xlabe=NULL, lty=c(1,2,3), col=c(1,2,3), lwd=c(1,1,1),
             ...)
{
   if(!is.null(kmm.ob$nsim)){	
	if (type == "Kmm.n") {        
	    cosa <- cbind(kmm.ob$kmm.n, 
                  t(apply(kmm.ob$kmmsim.n, 1, quantile, c(q,1-q), na.rm=T)))
            cosamean <- t(apply(kmm.ob$kmmsim.n, 1, mean, na.rm=T))
	    if(is.null(ylabe)) ylabe <- expression(normalized.K[mm](r)) 
            
	}
        if (type == "Kmm") {        
	    cosa <- cbind(kmm.ob$kmm, 
                  t(apply(kmm.ob$kmmsim, 1, quantile, c(q,1-q), na.rm=T)))
            cosamean <- t(apply(kmm.ob$kmmsim, 1, mean, na.rm=T))
            if(is.null(ylabe)) ylabe <- expression(K[mm](r))
            
	
	}	
       if(is.null(maine)) maine <- kmm.ob$dataname	
       if(is.null(xlabe)) xlabe <- "distance"
       matplot(kmm.ob$r, cosa, type="l", lty=c(lty[1],lty[2],lty[2]), 
           col=c(col[1],col[2],col[2]), lwd=c(lwd[1],lwd[2],lwd[2]),
           main= maine,xlab= xlabe, ylab= ylabe, ylim=ylime, xlim=xlime,
           mgp=c(2.5,1,0), add=add)
        if (kmean==TRUE) lines(kmm.ob$r, cosamean, lty=lty[3], col=col[3],lwd=lwd[3])
   }
   
   if(is.null(kmm.ob$nsim)){
     if (type == "Kmm.n") {        
	    cosa <- kmm.ob$kmm.n 
	    if(is.null(ylabe)) ylabe <- expression(normalized.K[mm](r)) 
            
	}
        if (type == "Kmm") {        
	    cosa <- kmm.ob$kmm
	     if(is.null(ylabe)) ylabe <- expression(K[mm](r))
	}	
       if(is.null(maine)) maine <- kmm.ob$dataname	
       if(is.null(xlabe)) xlabe <- "distance"
       matplot(kmm.ob$r, cosa, type="l", lty=lty[1], 
           col=col[1], lwd=lwd[1],
           main= maine,xlab= xlabe, ylab= ylabe, ylim=ylime, xlim=xlime,
           mgp=c(2.5,1,0), add=add)
        
   }
   
}

print.ecespa.kmm <- function(kmm.ob, ...)
{
    if(!is.null(kmm.ob$nsim)){
	cat("mark-weighted K function computed on",kmm.ob$dataname,"\nand tested with",
               kmm.ob$nsim,"random permutations of marks.\n")
    }
    if(is.null(kmm.ob$nsim)){
	cat("mark-weighted K function computed on",kmm.ob$dataname,"\n")
   }	
cat("plot it to see the result.\n")
}



###############################################################################
###############################################################################
######## FUNCIÓN Kmulti.ls (original ecespa) ###############################################
######## Mantener siempre. Es llamada por K012 y K1K2 ################################
# !!!!!!!!!!!!! Modificar función: que detecte que X es ppp con marcas discretas 


Kmulti.ls <- function (X, I, J, r = NULL, corre = "isotropic") 
{
    n1 = sum(I)
    n2 = sum(J)
    cosa12 <- Kmulti(X, I, J, r, correction = corre)
    cosa21 <- Kmulti(X, J, I, r, correction = corre)
    K12ls = ((n2 * cosa12[[3]]) + (n1 * cosa21[[3]]))/(n1 + n2)
    cosa12[[3]] = K12ls
    return(cosa12)
}



###############################################################################
######### FUNCIÓN marksum (original ecespa) ##############################################################
## nueva versión de marksum que devuelve objetos de clase "ecespa.marksum"
# modificado el rd
#13/06/08
 #modificación a marksum para que acepte ventanas poligonales y 
 # para que dibuje apoyado en smooth.ppp. 
 #CHANGES: Ha cambiado bastante plot.marksum. AHora no es interpolación sino smooth (se consigue
 # poder aplicarlo a ventanas irregulares)
 ###############################################################################
 
marksum <- function (mippp, R = 10, nx = 30, ny = 30) 
{
    require(spatstat)
    dataname <- deparse (substitute(mippp))
    verifyclass(mippp, "ppp")
    if (is.marked(mippp) != T) 
        stop("marksum only implemented for **marked** patterns")
   #generate a grid of dimensions nx ny inside the window of the pattern
    grid <- gridcenters(mippp$window, nx, ny)
          options(warn= -1) # do not warning on trimming points
    grid.ppp <- ppp(x = grid$x, y = grid$y, marks = rep(0, length(grid$x)), 
                         window = mippp$window)[mippp$window]
         options(warn= 0)
    prueba <- superimpose(grid.ppp, mippp)

   ##now count all the points within a distance R of the grid points
    marksum <- markstat(prueba, sum, R = R)
    marksum <- marksum[1:grid.ppp$n] #we want only the sums around the grid points!
    pointsum <- markstat(prueba, length, R = R)
    pointsum <- pointsum[1:grid.ppp$n] #we want only the counts around the grid points!
## correction of excess points in pointsum. We must substract the grid points summed
## in the previous markstat
    minus <- markstat(grid.ppp, length, R = R)
    pointsum <- pointsum - minus
    normalized.marksum <- marksum/pointsum
    normalized.marksum[marksum == 0] <- 0
    result <- list(normalized = normalized.marksum, marksum = marksum,
                     pointsum = pointsum,  minus = minus, grid = grid, grid.ppp=grid.ppp,nx = nx, 
                     ny = ny, R = R, window = mippp$window, dataname=dataname)
   class(result) <- c("ecespa.marksum", class(result))
   return(result)
}


plot.ecespa.marksum <- 
function(mkobject, what="normalized", contour=F,...)
{
 require (spatstat)
 if(what=="normalized") {cosa <- mkobject$normalized; what="normalized mark-sum"}
 if(what=="pointsum") {cosa <-  mkobject$pointsum; what="point-sum"}
 if(what=="marksum") {cosa <-  mkobject$marksum; what="mark-sum"}
 
 plot(smooth.ppp(setmarks(mkobject$grid.ppp, cosa),...), main="")
 title(main=paste(mkobject$dataname,"\n",noquote(what), "measure; R=", mkobject$R))
 if (contour==TRUE) contour(smooth.ppp(setmarks(mkobject$grid.ppp, cosa),...), add=T)
}

print.ecespa.marksum <- function(obj,...)
{
cat("Mark-sum measure of the dataset", obj$dataname,  "\n computed for a radius R=",obj$R, 
", with a grid of", obj$nx, "x",obj$ny,".\n")
cat("Plot it to see the result.\n")
}


###############################################################################
######### FUNCIÓN p2colasr (original ecespa) ##############################################
######### Mantener siempre: es llamada por dixon 2002 #######################################

p2colasr <- function (Z, nsim = length(Z)) 
{
    if (Z[1] < 0) {
        p1 = rank(Z)[1]/(nsim + 1)
        ptodos = 1 - (rank(Z)/(nsim + 1))
        if (sum(ptodos[ptodos < p1]) == 0) 
            p2 = p1
        else p2 = p1 + max(ptodos[ptodos < p1])
    }
    else {
        p1 = 1 - rank(Z)[1]/(nsim + 1)
        ptodos = rank(Z)/(nsim + 1)
        if (sum(ptodos[ptodos < p1]) == 0) 
            p2 = p1
        else p2 = p1 + max(ptodos[ptodos < p1])
    }
    return(p2)
}









###############################################################################
###############################################################################
######## FUNCIONES pc.estK  y Kclust (original ecespa) ###############################################
######## Mantener por  usarse  en el libro original         #####################################
# ######  Modificado pc.estK.rd  para hacer referencia a ipc.estK ####################################

pc.estK <-function (Kobs, r, sigma2 = NULL, rho = NULL) 
{
    theta = c(sigma2, rho)
    if (is.null(theta)) 
        theta = c(1, 1)
    D.theta = function(theta, Kobs, r) {
        sum((Kobs^(0.25) - Kclust(r, theta[1], theta[2])^(0.25))^2)
    }
    nsfit = optim(par = theta, fn = D.theta, Kobs = Kobs, r = r)
    return(list(sigma2 = nsfit$par[1], rho = nsfit$par[2]))
}



Kclust <- function (r, sigma2, rho) 
{
    (pi * r^2) + ((1 - exp(-(r^2)/(4 * sigma2)))/rho)
}


###############################################################################
###############################################################################
######## FUNCION sim.poissonc (original ecespa) ###############################################
######## Mantener por  usarse  en el libro original         #####################################
### Hecha referencia a que sólo simula en ventanas cuadradas ##################################
# Modificado sim.poissonc.rd  para hacer referencia a rIPCP ####################################

sim.poissonc <- function (x.ppp, rho, sigma) 
{
    mirpoispp = function(lambda, win) {
        n = lambda * area.owin(win)
        x = runif(n, win$xrange[1], win$xrange[2])
        y = runif(n, win$yrange[1], win$yrange[2])
        return(ppp(x, y, window = win))
    }
    ventana = x.ppp$window
    npadres = rho * area.owin(ventana)
    media = x.ppp$n/npadres
    padres = mirpoispp(lambda = rho, win = ventana)
    padres.pois = rpois(npadres, media)
    hijosx = rep(padres$x, padres.pois)
    hijosy = rep(padres$y, padres.pois)
    desviacionesx = rnorm(sum(padres.pois), mean = 0, sigma)
    desviacionesy = rnorm(sum(padres.pois), mean = 0, sigma)
    hijosx = hijosx + desviacionesx
    hijosy = hijosy + desviacionesy
    for (i in 1:length(hijosx)) {
        if (hijosx[i] < ventana$x[1]) 
            hijosx[i] = ventana$x[2] - (abs(hijosx[i]) - abs(ventana$x[1]))
        if (hijosx[i] > ventana$x[2]) 
            hijosx[i] = ventana$x[1] + (abs(hijosx[i]) - abs(ventana$x[2]))
        if (hijosy[i] < ventana$y[1]) 
            hijosy[i] = ventana$y[2] - (abs(hijosy[i]) - abs(ventana$y[1]))
        if (hijosy[i] > ventana$y[2]) 
            hijosy[i] = ventana$y[1] + (abs(hijosy[i]) - abs(ventana$y[2]))
    }
    return(ppp(x = hijosx, y = hijosy, window = ventana))
}





####################################################################################################################
###################################################################################################################
#######          ipc.estK                                                             #####################################################################
#######                                                                                  ####################################################################
### ###   CHANGES: esta función es nueva para ecespa 1.1-0.          ###############################################################3
####### CHANGES: esta función debería reemplazar a pc.estK de ecespa 1.0-3 (que debería eliminarse, además de Kclust) ###########################
##################################################################################################################

ipc.estK <- function (mippp, lambda=NULL, correction="iso", r=NULL, sigma2 = NULL, rho = NULL, q = 1/4, p = 2) 
{
    Kclust <- function (r, sigma2, rho) {
       (pi * r^2) + ((1 - exp(-(r^2)/(4 * sigma2)))/rho)
    }
   
     D.theta <- function(theta, Kobs, r) {
        sum((Kobs^q - Kclust(r, theta[1], theta[2])^q)^p)
     }

     lambdaname <- deparse(substitute(lambda))
     if(is.null(lambda)){
           lambda <- predict(ppm(mippp), type="trend")
           lambdaname <- NULL
     }
     Kobs <- Kinhom(mippp, r=r, correction=correction, lambda=lambda)

     if(is.null(r)) r <- Kobs$r 
     Kobs <- Kobs[[3]]
     theta <- c(sigma2, rho)
     if (is.null(theta)) theta <- c(1, 1)

     nsfit <- optim(par = theta, fn = D.theta, Kobs = Kobs, r = r)

     Kfit <- Kclust(r, nsfit$par[1], nsfit$par[2])
     dataname <- deparse(substitute(mippp))
    
     dtheta <- sum((Kobs^q - Kfit^q)^p)

     result <- list(sigma2 = nsfit$par[1], rho = nsfit$par[2], d.theta=dtheta, Kobs=Kobs,
                    Kfit=Kfit, r=r, data=mippp, lambda=lambda, dataname=dataname,
                    lambdaname=lambdaname, p=p, q=q)
     class(result) <- c("ecespa.minconfit", class(result))
     return(result)
}



#S3 method para dibujar resultados de ipcestK (y otros de clase "miconfit")

plot.ecespa.minconfit <- function(x, type="L", add=F, xlim=NULL, ylim=NULL,
                            lwd=c(1,1),lty=c(1,2), col=c(1,2), main=NULL, ...){
   if(type!="L"){
       if(add==F){
           plot(x$r, x$Kobs, xlab="r", ylab="K(r)", type="l", xlim=xlim, ylim=ylim,
                     lwd=lwd[1], lty =lty[1], col=col[1], ...)
           if(!is.null(main)) title(main=main) else title(main=x$dataname)
       }
       if(add==T) lines(x$r, x$Kobs, lwd=lwd[1], lty =lty[1], col=col[1])
       lines(x$r, x$Kfit, lwd=lwd[2], lty =lty[2], col=col[2])
       print("dashed line is K fited")
    }
    else if(type =="L"){
            if(add==F) {
            if(is.null(ylim)) ylim= range(c(range(sqrt(x$Kobs/pi)- x$r), range(sqrt(x$Kfit/pi)- x$r)))
            plot(x$r, sqrt(x$Kobs/pi)- x$r, xlab="r", ylab="L(r)" , 
                 type="l", xlim=xlim, ylim=ylim,
                 lwd=lwd[1], lty =lty[1], col=col[1], ...)
            if(!is.null(main)) title(main=main) else title(main=x$dataname)
         }
         if(add==T) lines(x$r, sqrt(x$Kobs/pi)- x$r, lwd=lwd)
         lines(x$r, sqrt(x$Kfit/pi)- x$r, lwd=lwd[2], lty =lty[2], col=col[2])
         print("dashed line is L fited")
   }
    invisible(NULL)
}




# S3 method para imprimir objetos de tipo "ecespa.minconfit"

print.ecespa.minconfit <-function (x, ...) {
    cat(paste("Minimum contrast fit ", "(", "object of class ", 
        dQuote("ecespa.minconfit"), ")", "\n", sep = ""))
    da <- x$dataname
    cat(paste("Model:PCP\n"))
    if(is.null(x$lambdaname)){
    cat(paste("Fitted by matching theoretical Kest function to", 
            x$dataname))
    }
    else cat(paste("Fitted by matching theoretical Kinhom function to", 
            x$dataname, "with lambda estimated by ", x$lambdaname))
    cat("\n")
    cat("Parameters fitted by minimum contrast ($sigma2 and $rho):\n")
    print(c(sigma2=x$sigma2, rho=x$rho,...))
    cat(paste("Domain of integration:", "[", signif(min(x$r), 
        4), ",", signif(max(x$r), 4), "]\n"))
    cat(paste("Exponents:", "p=", paste(signif(x$p, 3), ",", 
        sep = ""), "q=", signif(x$q, 3), "\n"))
    invisible(NULL)
}




####################################################################################################################
###################################################################################################################
###### FUNCIÖN       rIPCP #########################################################################
###########################################################################
## function to simulate  (I)PCP from a fitted object of class mi.minconfit
## nombre original : rpoisCipp2
# 12/06/08: introducido salvar del original  polygonal window (xw) y reponer al resultado final
# para que pueda ser calculada la iso correction. Comparar con trans y ventana mask
#!!!!!!!!!!! #Rd; faltan ejemplos




rIPCP <- function (x, lambda=NULL, type=1,lmax = NULL, win = owin(c(0, 1), c(0, 1)), ...) {
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

        # mean retention probability of the thinning process (se podría meter algún corrector "+..." para afinar mejor)

          probm <- mean(eval.im(lambda/lmax))


        # "prethinnig" required n to obtain the original n after thinning
          n.preth <- x.ppp$n / probm

        # "prethinning" required mu and rho to obtain the original mu and rho after thinning
          mu.preth <- n.preth*mu/x.ppp$n #ojo, directamente mu.preth = mu/probm
          rho.preth <- rho/probm

        #pre-thining point pattern:
 
        if(type==1)  X <- rThomas (win=x.ppp$w, kappa=rho, sigma=sigma, mu.preth) #simula resultados exactos al modelo ajustado
 
        if(type==2)  X <- rThomas (win=x.ppp$w, kappa=rho.preth, sigma=sigma, mu=mu ) #simula resultados más dispersos"

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
















## Función para calcular envueltas simulando partir de modelos ajustados, tanto de tipo "ppm" como de tipo
## "ecespa.minconfit"
## si simu="both", simula los dos modelos;
## las envueltas de K12 se construyen dejando fijo el patrón 1 y simulando el patrón 2
## las envueltas de K21 se construyen dejando fijo el patrón 2 y simulando el patrón 1
## Nombre original = Kci.envYYYYY. Se ha eliminado el arg "nlarge" de Kmulti.inhom
## podrían meterse más argumentos para controlar la simulación del rM-H algorithm
## 12/06/08: introducido argumento "ngrid" para permitir seleccionar el grid de los trends generados
##            a partir de modelos ppm. Efectos importantes en datos de urkiola ngrid=50 vs. ngrid>=200.
## 12/06/08: introducido arg "verbose" para controlar los mensajes de rmh (default=F).
## 29/07/08: sacado arg "verbose" de los argumentos y puesto en la primera línea de la función.
##                puede modificarse editando la función.
## 29/07/08: Cambiado argumento "Kmean" por "kmean" en plot.ecespa.kci, por coherencia con plot.ecespa.kmm
## 30/07/08: Hecho RD. Incluye Ki y plot.ecespa.kci



Kci <- function(mod1,mod2,correction="trans", nsim=99,
                         ngrid=200, nrep=1e5, r=NULL, simu="both", 
                         spctype=1)
{
   verbose <- FALSE # 
   ## datos básicos
   modnamea <- deparse(substitute(mod1))
   modnameb <- deparse(substitute(mod2))
   
   if(inherits(mod1,"ppm")){  
      I.ppp <- mod1$Q$data
      lambdaI <- predict(ppm(mod1$Q$data, mod1$trend), type="trend", ngrid=ngrid)
      Isim <- "ppm"
      dataname.a <- mod1$call[[2]]
   }
   if(inherits(mod1,"ecespa.minconfit")){  
      I.ppp <- mod1$data
      lambdaI <- mod1$lambda
      Isim <- "spc"
      dataname.a <- mod1$dataname
   }

   if(inherits(mod2,"ppm")){  
      J.ppp <- mod2$Q$data
      lambdaJ <- predict(ppm(mod2$Q$data, mod2$trend), type="trend", ngrid=ngrid)
      Jsim <- "ppm"
      dataname.b <- mod2$call[[2]]
   }
   if(inherits(mod2,"ecespa.minconfit")){  
      J.ppp <- mod2$data
      lambdaJ <- mod2$lambda
      Jsim <- "spc"
      dataname.b <- mod2$dataname
   }

   ## Cálculo de los Ki de cada patrón
   Kia <- Kinhom(I.ppp, lambdaI, correction=correction, r=r)
   mi.r <- Kia$r
   Kia <- Kia[[3]]
   Kib <- Kinhom(J.ppp, lambdaJ, correction=correction, r = mi.r)
   Kib <- Kib[[3]]
   
   ## generación del patrón multivariado y cálculo de Kci observada
   ## ojo: por el borde y la inhomogeneidad no es lo mismo Kab que Kba
   
   IJ.ppp <- superimpose(a=I.ppp, b=J.ppp)   
   Kci.ab.o <- Kcross.inhom(IJ.ppp, i="a", j="b", lambdaI, lambdaJ,
                            correction= correction, r=mi.r)[[3]]
   Kci.ba.o <- Kcross.inhom(IJ.ppp, i="b", j="a", lambdaJ, lambdaI,
                            correction= correction, r=mi.r)[[3]]

   Kia.s <- NULL ##para la K univariada del patrón que se simula
   Kib.s <- NULL
   Kci.ab.s <- NULL
   Kci.ba.s <- NULL
   # Isim.ppp <- I.ppp
   ## start simulations
   for (i in 1: nsim){ 
      progressreport(i,nsim)
      ## simulate from second model:
      if (Jsim=="ppm"){     Jsim.ppp <- rmh(mod2, start=list(x.start=J.ppp),
                                            control=list(p=1, nrep=nrep),
                                             verbose=verbose)}
      else if (Jsim=="spc") Jsim.ppp <- rIPCP (mod2, type=spctype)
      
      ## aseguramos que no haya NAs en el vector simulado de lambdas
      dentro <- !is.na(lambdaJ[Jsim.ppp, drop=F]) 
      Jsim.ppp <- Jsim.ppp[dentro] 

      ## simulated multivariate PP
      IJs.ppp <- superimpose(a=I.ppp, b=Jsim.ppp, W=I.ppp$w)
      IsJ.ppp <- IJs.ppp # si solo se simula J, creamos dos multipatrones iguales

      if(simu=="both"){
          if(Isim=="ppm"){    Isim.ppp <- rmh(mod1, start=list(x.start=I.ppp),
                                               control=list(p=1, nrep=nrep),
                                               verbose=verbose)}
          else if (Isim=="spc") Isim.ppp <- rIPCP (mod1, type=spctype)
          
          ## aseguramos que no haya NAs en el vector simulado de lambdas
          dentro <- !is.na(lambdaI[Isim.ppp, drop=F]) 
          Isim.ppp <- Isim.ppp[dentro]

          ## simulated multivariate PP:
          IsJ.ppp <- superimpose(a=Isim.ppp, b=J.ppp, W=I.ppp$w)
      } 

      ## K simuladas
      Kib.s <- cbind(Kib.s, Kinhom(Jsim.ppp, lambdaJ, correction=correction,
                                  r = mi.r, nlarge=Inf)[[3]])
      Kci.ab.s <- cbind(Kci.ab.s, Kcross.inhom(IJs.ppp, i="a", j="b",
                                     lambdaI, lambdaJ, r=mi.r,
                                     correction=correction)[[3]])
       Kci.ba.s <- cbind(Kci.ba.s, Kcross.inhom(IsJ.ppp, i="b", j="a",
                                     lambdaJ, lambdaI, r=mi.r,
                                     correction=correction)[[3]])
       if(simu=="both"){Kia.s <- cbind(Kia.s, Kinhom(Isim.ppp, lambdaI, 
                                                 correction=correction, r = mi.r, 
                                                 nlarge=Inf)[[3]])}
   }
   

   result <- list(r=mi.r, kia = Kia, kib=Kib, kci.ab.o=Kci.ab.o,
                kci.ba.o=Kci.ba.o, kci.ab.s=Kci.ab.s, kci.ba.s=Kci.ba.s,
                kib.s=Kib.s, kia.s=Kia.s, datanamea=dataname.a, datanameb=dataname.b,
		modnamea=modnamea, modnameb=modnameb, type="Kci")
   class(result) <- c("ecespa.kci", class(result))
   return(result)
}



## Función para calcular envueltas simulando partir de modelos ajustados, tanto de tipo "ppm" como de tipo
## "ecespa.minconfit"
## Nombre original = Ki.envYYYYY
## podrían meterse más argumentos para controlar la simulación del rM-H algorithm

#!!!!!!!!!!! #FALTA RD;


Ki <-function(mod1,correction="trans", nsim=99, ngrid=200,
                         nrep=1e5, r=NULL, spctype=1)
{
## datos básicos
 modnamea <- deparse(substitute(mod1))
 
if(inherits(mod1,"ppm")){  
   I.ppp= mod1$Q$data
   lambdaI = predict(ppm(mod1$Q$data, mod1$trend), type="trend", ngrid=ngrid)
   Isim="ppm"
   dataname.a <- mod1$call[[2]]
}
if(inherits(mod1,"ecespa.minconfit")){  
   I.ppp= mod1$data
   lambdaI = mod1$lambda
   Isim="spc"
   dataname.a <- mod1$dataname
}


## corrección de los problemas de descoordinación entre inside.owin y density en
## la versión 1.11.1 de spatstat( hay que asegurarse de que no hay puntos fuera de 
## las ventanas de intensidad
 #      dentro = !is.na(lambdaI[I.ppp, drop=F])
 #      I.ppp=I.ppp[dentro]
 #      dentro = !is.na(lambdaJ[J.ppp, drop=F])
 #      J.ppp=J.ppp[dentro]
##### FIN DE LA CORRECCIÓN

  ## Cálculo de los Ki de cada patrón
   Kia=Kinhom(I.ppp, lambdaI, correction=correction, r=r)
   mi.r = Kia$r
   Kia = Kia[[3]]
   Kia.s= NULL
   Isim.ppp = I.ppp
      for (i in 1: nsim){ 
           progressreport(i,nsim)
        ## aseguramos que no haya NAs en el vector simulado de lambdas
            if(Isim=="ppm"){      Isim.ppp = rmh(mod1, start=list(x.start=I.ppp), 
                                                control=list(p=1, nrep=nrep), verbose=FALSE) }
            else if (Isim=="spc") Isim.ppp = rIPCP (mod1, type=spctype)

            dentro = !is.na(lambdaI[Isim.ppp, drop=F]) 
            Isim.ppp = Isim.ppp[dentro]
            ## K simuladas
            Kia.s= cbind(Kia.s, Kinhom(Isim.ppp, lambdaI, 
                           correction=correction, r = mi.r, nlarge=Inf)[[3]])
    }
    
   result=(list(r=mi.r, kia = Kia, kia.s=Kia.s, lambda=lambdaI, datanamea=dataname.a, modnamea=modnamea, type="Ki"))
   class(result)<-c("ecespa.kci", class(result))
   return(result)
}





# S3 method para objetos resultantes de  KCI

plot.ecespa.kci <- function(kci.ob, type=1, q=0.025, kmean=TRUE, add=F,
	                             maine=NULL, xlabe=NULL, ylabe=NULL, xlime=NULL, ylime=NULL,    
				     lty=c(1,2,3), col=c(1,2,3), lwd=c(1,1,1), ...)
{
  
   if (type == 1) {                      ## univariado típico L1
       cosa <- cbind(kci.ob$kia, 
                  t(apply(kci.ob$kia.s,1,quantile, c(q,1-q),na.rm=T)))
       cosamean <- sqrt(apply(kci.ob$kia.s,1,mean,na.rm=T)/pi)- kci.ob$r
       cosa <- sqrt(cosa/pi)-kci.ob$r
       if(is.null(ylabe)) ylabe <- expression(hat(L)[1])
       if(is.null(maine)) maine <- kci.ob$datanamea
   }
   if (type == 2) {                      ## univariado típico L2
       cosa <- cbind(kci.ob$kib, 
                  t(apply(kci.ob$kib.s,1,quantile, c(q,1-q),na.rm=T)))
       cosamean <- sqrt(apply(kci.ob$kib.s,1,mean,na.rm=T)/pi)- kci.ob$r
       cosa <- sqrt(cosa/pi)-kci.ob$r
       if(is.null(ylabe)) ylabe <- expression(hat(L)[2])
       if(is.null(maine)) maine <- kci.ob$datanameb
   }
   if (type == 12) {                      ## bivariado típico L12
       cosa <- cbind(kci.ob$kci.ab.o, 
                  t(apply(kci.ob$kci.ab.s,1,quantile, c(q,1-q),na.rm=T)))
       cosa <- sqrt(cosa/pi)-kci.ob$r
       cosamean <- apply(kci.ob$kci.ab.s,1,mean,na.rm=T)
       cosamean <- sqrt(cosamean/pi)-kci.ob$r
       if(is.null(ylabe)) ylabe <- expression(hat(L)[12])
       if(is.null(maine)) maine <- paste(kci.ob$datanamea, " vs. ", kci.ob$datanameb) 
   }
   
   if (type == 21) {                       ## bivariado típico L21
       cosa <- cbind(kci.ob$kci.ba.o, 
                  t(apply(kci.ob$kci.ba.s,1,quantile, c(q,1-q),na.rm=T)))
       cosa=sqrt(cosa/pi)-kci.ob$r
       cosamean <- apply(kci.ob$kci.ba.s,1,mean,na.rm=T)
       cosamean <- sqrt(cosamean/pi)-kci.ob$r

       if(is.null(ylabe)) ylabe <- expression(hat(L)[21])
       if(is.null(maine)) maine <- paste(kci.ob$datanamea, " vs. ", kci.ob$datanameb) 
  }

   if (type == 112) {                       ## segregación primer ppp (K1-K12)
       d1.12.o <- kci.ob$kia - kci.ob$kci.ab.o
       d1.12.s <- kci.ob$kia - kci.ob$kci.ab.s
       cosa <- cbind(d1.12.o, t(apply(d1.12.s,1,quantile, c(q,1-q),na.rm=T)))
       cosamean <- apply(d1.12.s,1,mean,na.rm=T)
       if(is.null(ylabe)) ylabe <- expression(hat(K)[1]- hat(K)[12])
       if(is.null(maine)) maine <- paste(kci.ob$datanamea, " vs. ", kci.ob$datanameb) 
   }

   if (type == 221) {                       ## segregación segundo ppp (K2-K21)
       d2.21.o = kci.ob$kib - kci.ob$kci.ba.o
       d2.21.s = kci.ob$kib.s - kci.ob$kci.ba.s
       cosa <- cbind(d2.21.o, t(apply(d2.21.s,1,quantile, c(q,1-q),na.rm=T)))
       cosamean <- apply(d2.21.s,1,mean,na.rm=T)
       if(is.null(ylabe)) ylabe <- expression(hat(K)[2]- hat(K)[21])
       if(is.null(maine)) maine <- paste(kci.ob$datanamea, " vs. ", kci.ob$datanameb) 

   }
   if(is.null(xlabe)) xlabe <-"distance"
   matplot(kci.ob$r,cosa, type="l", lty=c(lty[1],lty[2],lty[2]), 
           col=c(col[1],col[2],col[2]), lwd=c(lwd[1],lwd[2],lwd[2]),
           main= maine,xlab= xlabe, ylab= ylabe, ylim=ylime, xlim=xlime,
           mgp=c(2.5,1,0), add=add)
  if (kmean==TRUE) lines(kci.ob$r, cosamean, lty=lty[3], col=col[3],lwd=lwd[3])

}





print.ecespa.kci <- function(kci.ob, ...)
{
   if(kci.ob$type=="Kinhom.log") { 
        if (kci.ob$probname=="null" & is.null(kci.ob$modtrend)){
            cat ("K function and envelopes of", kci.ob$datanamea,
                   "\nsimulated by", kci.ob$nsim, "random thinnings of the original pattern,
		    \n using ecespa function Kinhom.log\n\n")
        }
   
        if( kci.ob$probname=="null" & !is.null(kci.ob$modtrend)){
            cat ("inhomogeneous K function and envelopes of", kci.ob$datanamea, 
		   "\nsimulated by", kci.ob$nsim, "random thinnings of the original pattern, 
	            \nusing function ecespa Kinhom.log with spatial trend as", deparse(kci.ob$modtrend),"\n\n")
        }

        if(kci.ob$probname!="null" & !is.null(kci.ob$modtrend)){
            cat ("inhomogeneous K function and envelopes of", kci.ob$datanamea, 
                   "\nsimulated by", kci.ob$nsim, "thinnings of the original pattern,\n", 
	           "according to the probability vector",kci.ob$probname, 
	           ",\n using ecespa function Kinhom.log with spatial trend as", deparse(kci.ob$modtrend),"\n\n")
        }

        if( kci.ob$probname!="null" & is.null(kci.ob$modtrend)){
            cat ("K function and envelopes of", kci.ob$datanamea, 
                   "\nsimulated by", kci.ob$nsim, "thinnings of the original pattern, 
	            \naccording to the probability vector",kci.ob$probname, 
	           ",\n using ecespa function Kinhom.log\n\n")
        }
   }

   if(kci.ob$type=="Ki"){
      cat("K function and envelopes of", kci.ob$datanamea,
            "\nbuilt by", kci.ob$nsim,"simulations from", kci.ob$modnamea,
             "model,\n using ecespa function Ki\n")
  }

  if(kci.ob$type=="Kci"){
   cat("Cross K functions and envelopes of", kci.ob$datanamea, "and", kci.ob$datanameb,
   "\nbuilt by", kci.ob$nsim,"simulations from", kci.ob$modnamea,"and", kci.ob$modnameb,
    "models,\n using ecespa function Ki\n")
  }
}










#######################################################################################
######################################################################################
#####Kinhom.log. Derivado de Kinhom3.log (c:\2007\trabajos\URKIOLA) ######################################
##### Función y RD elaborados el 30/08/2008



Kinhom.log <- function(A, lambda=NULL, mod=NULL, lifemark="0", prob=NULL,
                                r=NULL, nsim=99, correction="trans", ngrid=200){
  
   require(spatstat)
    dataname <- deparse(substitute(A))
    if (is.null(prob)) probname <- "null" else probname <- deparse(substitute(prob))
   Avivos <- unmark(A[A$marks == lifemark])
   Atotal <- unmark(A)
   
   # if there is not probability vector, give everyone the same probability
   # (get it from the frequency of "live" cases):
   
   if (is.null(prob)) prob <- Avivos$n / Atotal$n
   
   # if there is no lambda nor mod, generate a constant lambda
   # this is useful for homogeneous patterns:
   if (is.null(lambda) & is.null (mod)){
         mod <- ppm(Avivos, ~1, Poisson())
         lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
   }
   
   # if there is not lambda and there is mod, predict lambda from mod
   if (is.null(lambda) & !is.null(mod)){
         mod <- update(mod, Q=Avivos)
         lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
    }
   
   # if there are both lambda and mod, use mod and forget lambda
    if (!is.null(lambda) & !is.null(mod)){
         mod <- update(mod, Q=Avivos)
         lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
    }

# if there is lambda but not model, approximate the model with lambda as a covariate
    if (!is.null(lambda) & is.null(mod)){
         mod <- ppm(Avivos,~1 + lam, covariates=list(lam=lambda))
	 lambda <- predict.ppm(mod, type="trend", ngrid=ngrid)
    }

    ## recálculo de los que están dentro de la ventana que genera density.ppp
    ## son los problemas típicos de decordinación con inside.owin en la versión
    ## 1.11.1 de spatstat

    dentro <- !is.na(lambda[Avivos,drop=F])
    Avivos <- Avivos[dentro]

    ## FIN DEL RECÁLCULO/CORRECCIÓN

    Kao <- Kinhom(Avivos , lambda, correction=correction, r=r)
    mi.r <- Kao$r
    Kao <- Kao[[3]] # 3 es el lugar de la lista resultante de Kinhom que ocupa la K observada
    Kas <- NULL
    
    ##simulations
    for(i in 1:nsim){
      progressreport(i,nsim)
      Avivos <- rthin(Atotal, prob)  
      mod <- update(mod, Q=Avivos)
      lambda <- predict.ppm(mod, type="trend")
      
      dentro <- !is.na(lambda[Avivos,drop=F])
      Avivos <- Avivos[dentro]
      
      Kas=cbind(Kas, Kinhom(Avivos , lambda,  correction=correction, 
                             r=mi.r)[[3]]) 

   }
   result <- (list(r=mi.r, kia = Kao, kia.s=Kas,  datanamea=dataname, type="Kinhom.log", probname=probname, modtrend=mod$trend, nsim=nsim))
   class(result)<-c("ecespa.kci", class(result))
   return(result)
   
}





















########################################################################################################
########################################################################################################
# COnjunto de funciones para calcular el test de Syrjala. 
########################################################################################################
# El cálculo se hace invocando "syrjala(ppp1,ppp2)
# pp1 y pp2 son patrones marcados, -idénticos en las coordenadas-
# ya que representan diferentes observaciones (sexos, especies, abundancias en dif. años)
# en los mismos _puntos_

#28/06/08 intento de remplazar el for loop de intercambio de valores por un lapply
#28/06/08 todas las funciones dentro de la principal ("syrjala.test")
#28/06/08 calculo de gammaf siguiendo al pie de la letra la deficnición de Srykjala 1996:77
#28/06/08 desarrollo de un print.ecespa.syrjala y plot.ecespa.syrjala
# OJO: en la implementación del trabajo de aves nocturnas, lo que se permutaba era la abundancia
# bruta, en vez de la abundancia normalizada, que es lo que se permuta ahora (de acuerdo con
# el trabajo de Syrjala)

syrjala.test <- function(ppp1, ppp2, nsim=999){
   
   require(spatstat)
   datanames <- c(deparse(substitute(ppp1)), deparse(substitute(ppp2)))   
   
   # function gamma: computes cummulative distribution function at the location(xk, yk)
   gammaf <- function(cosa.ppp){
              # First, normalize the observed density data
		#cosa.ppp$marks <- cosa.ppp$marks/sum(cosa.ppp$marks)
   
	      # Second, compute cumulative distribution at each location "k"
		gamma <- NULL
		for(k in 1: cosa.ppp$n){
    		    gamma <- c(gamma,
                        		sum(cosa.ppp$marks[cosa.ppp$x<=cosa.ppp$x[k] & cosa.ppp$y<=cosa.ppp$y[k]]))
                }
                return(gamma)
        }

    # function psi: computes the squared difference between the two cumulative distribution functions 
    psi <- function(ppp1, ppp2){
		gamma1 <- gammaf (ppp1)
		gamma2 <- gammaf (ppp2)
		psi <- sum( (gamma1-gamma2)^2)
		return(psi)
	}


    # function psimean: computes the mean psi (i.e among different psi's computed on cummulative distributions
    # relative to the different four corners of the study region). 
    psimean <- function(ppp1, ppp2){

			psi1 <- psi(ppp1,ppp2)
			psi2 <- psi(affine(ppp1, mat=diag(c(-1,-1))), affine(ppp2, mat=diag(c(-1,-1))))
			psi3 <- psi(affine(ppp1, mat=diag(c(1,-1))), affine(ppp2, mat=diag(c(1,-1))))
			psi4 <- psi(affine(ppp1, mat=diag(c(-1,1))), affine(ppp2, mat=diag(c(-1,1))))
			psimean <- (psi1 + psi2 + psi3 + psi4)/4
			return(psimean)
	}

      #function rpsi: pairwise permutation of the normalized density observation
      ## OJO: en la implementación original para las nocturnas, lo que se permutó fué la observación bruta!!      
      # pepesim es una lista de listas (cada elemento lleva los dos point patterns)
      rpsi <- function(pepesim){     
                         marcas <- cbind(ppp1$marks,ppp2$marks)
                         marcas <- t(apply(marcas,1,sample))
                         ppp1$marks <- marcas[,1]
                         ppp2$marks <- marcas[,2]
                       #cálculo del psi simulado
                       psi.sim <-  psimean (ppp1, ppp2)
        }
	
	
# Compute syrjala test and simulations:	
	
# First, normalize the observed density data
		ppp1$marks <- ppp1$marks/sum(ppp1$marks)
   	        ppp2$marks <- ppp2$marks/sum(ppp2$marks)
	

#Second, compute observed psi
               psi.obs <- psimean(ppp1,ppp2)
   

#Third, instead of the classical loop, generate a list of pairs of populations, 
# permute  pairwisely the density and compute psi sim

     pepesim <- list(pepesim=list(ppp1=ppp1, ppp2=ppp2))
     pepesim <- rep(pepesim, nsim) 
     psi.sim <- sapply(pepesim, rpsi)
     
       

   result <- list(psi.obs=psi.obs, psi.sim=psi.sim, datanames=datanames, nsim=nsim)
   class(result) <- c("ecespa.syrjala", class(result))
   return(result)


}


plot.ecespa.syrjala <-function(syrjala.ob, ...){
    hist(syrjala.ob$psi.sim, col="grey", border="grey", main="", xlab=expression(psi),...)
    title(main=paste("Syrjala test\n for the difference between the spatial distributions of \n",
                            syrjala.ob$datanames[1], " and ", syrjala.ob$datanames[2]))
    abline(v=syrjala.ob$psi.obs, lwd=3)
    print(syrjala.ob)
}

print.ecespa.syrjala <-function(syrjala.ob, ...){
    cat("Syrjala test for the difference between the spatial distributions of \n",
                            syrjala.ob$datanames[1], "and ", syrjala.ob$datanames[2], ", based on", syrjala.ob$nsim,"simulations\n\n")
    cat("   psi:     ",syrjala.ob$psi.obs,"\n")
    if(mean(syrjala.ob$psi.sim) <= syrjala.ob$psi.obs)
       cat("   p-value: ", sum(syrjala.ob$psi.sim >= syrjala.ob$psi.obs)/(syrjala.ob$nsim+1), "\n\n")
    else cat("   p-value: ", sum(syrjala.ob$psi.sim <= syrjala.ob$psi.obs)/(syrjala.ob$nsim+1), "\n\n")
}

