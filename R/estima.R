#'Estimation of parameters of autologistic regression model for data on a grid
#'@param data dataset with the coordinates in the two first columns.
#'@param covariate1 spatio-temporal covariate. The covariate dataframe must have \code{dim(data)[1] = dim(covariate)[1]} (same numbers of individuals) and \code{dim(data)[1] = dim(covariate)[1] + 3} as the covaiate dataset must not contain coordinates, but must match the coodinates of the dataset; and \code{T-1} years (\code{T} is the number of years in the dataset \code{"data"}) as the model needs the first year to initialize. See \code{"User guides, package vignettes and other documentation"} the \code{"estima"} vignette.
#'@param covariate2 spatio-temporal covariate. The covariate dataframe must have \code{dim(data)[1] = dim(covariate)[1]} (same numbers of individuals) and \code{dim(data)[1] = dim(covariate)[1] + 3} as the covaiate dataset must not contain coordinates, but must match the coodinates of the dataset; and \code{T-1} years (\code{T} is the number of years in the dataset \code{"data"}) as the model needs the first year to initialize. See \code{"User guides, package vignettes and other documentation"} the \code{"estima"} vignette.
#'@param covariate3 spatio-temporal covariate. The covariate dataframe must have \code{dim(data)[1] = dim(covariate)[1]} (same numbers of individuals) and \code{dim(data)[1] = dim(covariate)[1] + 3} as the covaiate dataset must not contain coordinates, but must match the coodinates of the dataset; and \code{T-1} years (\code{T} is the number of years in the dataset \code{"data"}) as the model needs the first year to initialize. See \code{"User guides, package vignettes and other documentation"} the \code{"estima"} vignette.
#'@param norm  \code{"euclidean"}, \code{"inf"}, \code{"abs"}, \code{"lin"}. \code{norm = "euclidean"} by default. See vignette \code{Build}.
#'@param dx positive real : distance between sites on x-axis. \code{dx = 1} by default.
#'@param dy positive real : distance between sites on y-axis. \code{dy = 1} by default.
#'@param swpresent if \code{TRUE} the programm will test all possible neighborhood for the spatial autocorrelation (coefficient \code{rho1}) with parameters \code{vxmaxpresent}, \code{vymaxpresent}, \code{dx} and \code{dy}, otherwise the programm will test the neighborhood with the parameters \code{vxpresent} and \code{vypresent}. \code{swpresent = TRUE} by default.
#'@param swpast if \code{TRUE} the programm will test all possible neighborhood for the autoregression on on the sum of the \code{Zi,t-1} (coefficient Betapast) with parameters \code{vxmaxpast}, \code{vymaxpast}, \code{dx} and \code{dy}, otherwise the programm will test the neighborhood with the parameters \code{vxpast} and \code{vypast}. \code{swpast = TRUE} by default.
#'@param vxpresent positive real. Parameter of the ellipse for the tested neighborhood on x-axes in norm \code{"norm"} if \code{swpresent = FALSE}. If \code{swpresent = TRUE}, \code{vxpresent} will be the upper bound of the tested neighborhoods on x-axes in norm \code{norm}. See \code{swpresent}.
#'@param vypresent positive real. Parameter of the ellipse for the tested neighborhood on y-axes in norm \code{"norm"} if \code{swpresent = FALSE}. If \code{swpresent = TRUE}, \code{vypresent} will be the upper bound of the tested neighborhoods on y-axes in norm \code{norm}. See \code{swpresent}.
#'@param vxpast positive real. Parameter of the ellipse for the tested neighborhood on x-axes in norm \code{"norm"} if \code{swpast = FALSE}. If \code{swpast = TRUE}, \code{vxpast} will be the upper bound of the tested neighborhoods on x-axes in norm \code{norm}. See \code{swpast}. Only use if \code{pastcov = TRUE}.
#'@param vypast positive real. Parameter of the ellipse for the tested neighborhood on y-axes in norm \code{"norm"} if \code{swpast = FALSE}. If \code{swpast = TRUE}, \code{vypast} will be the upper bound of the tested neighborhoods on y-axes in norm \code{norm}. See \code{swpast}. Only use if \code{pastcov = TRUE}.
#'@param graph if \code{graph = TRUE}, the program will also return the plot of the dataset for the last time (and the year before if \code{estima = 3}). \code{graph = FALSE} by default.
#'@param pastcov boolen. If \code{pastcov = TRUE}, the function will use the past neighborhood as a covariate. See \code{"User guides, package vignettes and other documentation"} the \code{"estima"} vignette. \code{pastcov = FALSE} by default.
#'@param buildpres boolean which allow the use of a custom neighborhood matrix. \code{buildpres = NULL} by default.
#'@param buildpast boolean which allow the use of a custom neighborhood matrix. \code{buildpast = NULL} by default.
#'
#'
#'@return list : estimate parameters using the pseudo-likelihood.
#'
#'@details See \code{"User guides, package vignettes and other documentation"} the \code{"estima"} vignette.
#'
#'@export
#'
#'@import stats
#'@import utils
#'
#'@examples
#'
#'data <- plantillness
#'v <- which(data$NRang <= 10)
#'data <- data[v,]
#'v <- which(data$NCep <= 10)
#'data<-data[v,]
#'result <- estima(data = data)
#'
#'
#'
#'\donttest{
#'#Example in "lin" norm, with a fixed neighborhood :
#'
#'result <- estima(data = plantillness, norm = "lin",swpresent = FALSE,vxpresent = 3, vypresent = 4)
#'
#'
#'
#'#Example with a spatial covariate (adapted to the dimension of the dataset) :
#'
#'cov <- covplant[,1]
#'for (i in (1:(dim(plantillness)[2] - 4))){
#'  cov <- cbind(cov,covplant[,1])
#'}
#'result <- estima(data = plantillness,covariate1 = cov)
#'
#'
#'
#'#Example with the past neighborhood as covariate:
#'
#'result <- estima(data = plantillness,pastcov = TRUE)
#'
#'
#'
#'#Exemple with a custom neighborhood matrix
#'
#'custompres <- build(data = plantillness)
#'custompast <- build(data = plantillness, vx = 5,vy = 6)
#'result <- estima(data = plantillness,pastcov = TRUE,buildpres = custompres,buildpast = custompast)
#'
#'}
#'

estima <- function(data = 0,covariate1 = NULL, covariate2 = NULL,covariate3 = NULL,norm = "euclidean",vxpresent = 3, vypresent = 3,vxpast = 3, vypast = 3,  dx = 1,dy = 1,swpresent = TRUE,swpast = TRUE,graph = FALSE,pastcov = FALSE,buildpres = NULL,buildpast = NULL){


norm1 = norm
norm2 = norm

vxmaxpresent = vxpresent
vymaxpresent = vypresent

vxmaxpast = vxpast
vymaxpast = vypast



if (!pastcov){

  cov1 <- covariate1
  cov2 <- covariate2
  cov3 <- covariate3
  norme <- norm1
  vlmax <- vxmaxpresent
  vcmax <- vymaxpresent
  dl <- dx
  dc <- dy
  sweep <- swpresent
  vlp <- vxpresent
  vcp <- vypresent
  buil <- buildpres


  if (!(graph %in% c(TRUE,FALSE))){stop("graph must be TRUE or FALSE")}
  if (!(sweep %in% c(TRUE,FALSE))){stop("sweep must be TRUE or FALSE")}

  if (!("ggplot2" %in% rownames(installed.packages())) && graph){stop("ggplot2 package must be installed before.")}


  if (!(norme %in% c("euclidean","inf","abs","lin"))){stop("Norm must be \"euclidean\" \"inf\" \"abs\" or \"lin\" ")}

  if ((dl < 0 )){stop("dl must be positive")}
  if ((dc < 0)){stop("dc must be positive")}

  if (sweep && (vlmax < 0)){stop("vlmax must be positive")}
  if (sweep &&  (vcmax < 0 )){stop("vcmax must be positive")}
  if (!sweep && (vlp < 0)){stop("vlp must be positive")}
  if (!sweep && (vcp < 0)){stop("vcp must be positive")}




  data5 <- as.matrix(data)
  data2<-data5[,-(1:2)]

  if (is.null(cov1) && is.null(cov2) && is.null(cov3)){  #0 covariable
    cov = 0
  }

  if (!is.null(cov1) && is.null(cov2) && is.null(cov3)){  #1 covariable
    cov = 1
    x1 <- as.matrix(cov1)
    if (!((dim(x1)[1] == dim(data2)[1]) && (dim(x1)[2] == (dim(data2)[2]-1)))){stop("Dimensions of covariate1 does not match the dimensions of the dataset.")}
  }

  if (!is.null(cov1) && !is.null(cov2) && is.null(cov3)){  #2 covariable
    cov = 2
    x1 <- as.matrix(cov1)
    x2 <- as.matrix(cov2)
    if (!((dim(x1)[1] == dim(data2)[1]) && (dim(x1)[2] == (dim(data2)[2]-1)))){stop("Dimensions of covariate1 does not match the dimensions of the dataset.")}
    if (!((dim(x2)[1] == dim(data2)[1]) && (dim(x2)[2] == (dim(data2)[2]-1)))){stop("Dimensions of covariate2 does not match the dimensions of the dataset.")}
  }

  if (!is.null(cov1) && !is.null(cov2) && !is.null(cov3)){  #3 covariable
    cov = 3
    x1 <- as.matrix(cov1)
    x2 <- as.matrix(cov2)
    x3 <- as.matrix(cov3)
    if (!((dim(x1)[1] == dim(data2)[1]) && (dim(x1)[2] == (dim(data2)[2]-1)))) {stop("Dimensions of covariate1 does not match the dimensions of the dataset.")}
    if (!((dim(x2)[1] == dim(data2)[1]) && (dim(x2)[2] == (dim(data2)[2]-1)))) {stop("Dimensions of covariate2 does not match the dimensions of the dataset.")}
    if (!((dim(x3)[1] == dim(data2)[1]) && (dim(x3)[2] == (dim(data2)[2]-1)))) {stop("Dimensions of covariate3 does not match the dimensions of the dataset.")}
  }


  if (is.null(buil)){

    if (sweep){
      tours <- vlmax*vcmax
      modvoisB<-vector("list", vlmax*vcmax)
      k<-1
      for (vl1 in (1:vlmax)){
        for (vc1 in (1:vcmax)){
          modvoisB[[k]]<-build(data = data5,vx = vl1,vy = vc1,norm = norme)
          k<-k+1
        }
      }
    }
    else{
      tours <- 1
      modvoisB<-vector("list", 1)
      modvoisB[[1]]<-build(data = data5,vx = vlp,vy = vcp,norm = norme)
    }
  }
  else{
    tours <- 1
    modvoisB<-vector("list", 1)
    modvoisB[[1]]<-buil
    sweep = FALSE
  }


  ################################################################
  ################ Cas des balayges ##############################



  n <- dim(data5)[1]
  tp<-ncol(data2)-1
  zav<-as.matrix(data2[,-(tp+1)])
  z<-as.matrix(data2[,-1])
  mattotB<-NULL

  PLEstim<-vector("list",tours)
  rw <- 9 + (2*cov)
  w <- rw*tours
  a = rep(0,w)
  matpasB<-matrix(a,ncol = rw)

  if (sweep){    #cas balayage
    matpasB[,1]<-rep(1:vlmax,each=vcmax)
    matpasB[,2]<-rep(1:vcmax,vlmax)
  }

  else{    # cas non balayage
    matpasB[,1]<-vlp
    matpasB[,2]<-vcp
  }


  #Calcul glm

  if (cov == 0){
    for (k in 1:tours){

      data <- data2
      voisin <-  modvoisB[[k]]
      present <- (voisin)%*%z

      p <- glm(c(z)~1+c(present)+c(zav), family = binomial(link = "logit"))$coef

      ppast<-c(1,1,1)
      delta<-sum((p-ppast)*(p-ppast))
      i<-1

      while(delta>0.00001){
        ppast <-p

        x_beta <- p[1]
        present <-voisin%*%(z - (exp(x_beta+(p[3])*zav )/(1+exp(x_beta+(p[3])*zav ))))

        p<- glm(c(z)~1+c(present)+c(zav), family = binomial(link = "logit"))$coef
        delta<-sum((p-ppast)*(p-ppast))
        i<-i+1
      }
      s<-summary(glm(c(z)~1+c(present)+c(zav), family = binomial(link = "logit")))
      p<-s$coefficients[,1]
      x_beta <- p[1]
      present <-voisin%*%(z - (exp(x_beta+(p[3])*zav )/(1+exp(x_beta+(p[3])*zav ))))
      LPL<-0
      for (i in 1:length(z)){
        LPL<-LPL+z[i]*(p[1] + p[2]*present[i] + p[3]*zav[i]) + log(1 - (exp(p[1] + p[2]*present[i] + p[3]*zav[i])/(1+exp(p[1] + p[2]*present[i] + p[3]*zav[i]))))
      }


      PLEstim[[k]]<-list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
      matpasB[k,3:5]<-PLEstim[[k]][[1]]
      matpasB[k,6:8]<-PLEstim[[k]][[2]]
      matpasB[k,9]<-PLEstim[[k]][[4]]
    }
    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vx','vy','beta0','rho1','rho2','sdbeta0','sdrho1','sdrho2','LPL')
  }


  if (cov == 1){
    for (k in 1:tours){

      data <- data2
      voisin <- modvoisB[[k]]
      present <- (voisin)%*%z

      p <- glm(c(z)~1+c(x1)+c(present)+c(zav), family = binomial(link = "logit"))$coef

      ppast<-c(1,1,1,1)
      delta<-sum((p-ppast)*(p-ppast))
      i<-1

      while(delta>0.00001){
        ppast<-p
        beta0 <- p[1]
        beta1 <- p[2]
        rho2 <- p[4]
        x_beta<- (beta0) + (beta1)*x1

        present<-voisin%*%(z - (exp(x_beta+(rho2)*zav )/(1+exp(x_beta+(rho2)*zav ))))
        p<- glm(c(z)~1+c(x1)+c(present)+c(zav), family = binomial(link = "logit"))$coef
        delta<-sum((p-ppast)*(p-ppast))
        i<-i+1
      }
      s<-summary(glm(c(z)~1+c(x1)+c(present)+c(zav), family = binomial(link = "logit")))
      p<-s$coefficients[,1]
      beta0 <- p[1]
      beta1 <- p[2]
      rho2 <- p[4]
      x_beta<- (beta0) + (beta1)*x1

      present<-voisin%*%(z - (exp(x_beta+(rho2)*zav)/(1+exp(x_beta+(rho2)*zav))))
      LPL<-0
      for (i in 1:length(z)){
        LPL<-LPL+z[i]*(p[1] + p[2]*x1[i] + p[3]*present[i] + p[4]*zav[i]) + log(1 - (exp(p[1] + p[2]*x1[i] + p[3]*present[i] + p[4]*zav[i])/(1+exp(p[1] + p[2]*x1[i] + p[3]*present[i] + p[4]*zav[i]))))
      }

      PLEstim[[k]]<-list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
      matpasB[k,3:6]<-PLEstim[[k]][[1]]
      matpasB[k,7:10]<-PLEstim[[k]][[2]]
      matpasB[k,11]<-PLEstim[[k]][[4]]
    }
    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vx','vy','beta0','beta1','rho1','rho2','sdbeta0','sdbeta1','sdrho1','sdrho2','LPL')
  }

  if (cov == 2){
    for (k in 1:tours){
      data <- data2
      voisin <- modvoisB[[k]]
      present <- (voisin)%*%z

      p <- glm(c(z)~1+c(x1)+c(x2)+c(present)+c(zav), family = binomial(link = "logit"))$coef

      ppast<-c(1,1,1,1,1)
      delta<-sum((p-ppast)*(p-ppast))
      i<-1

      while(delta>0.00001){
        ppast<-p
        beta0 = p[1]
        beta1 = p[2]
        beta2 = p[3]
        rho2 = p[5]
        x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2

        present<-voisin%*%(z - (exp(x_beta+(rho2)*zav )/(1+exp(x_beta+(rho2)*zav ))))
        p<- glm(c(z)~1+c(x1)+c(x2)+c(present)+c(zav), family = binomial(link = "logit"))$coef
        delta<-sum((p-ppast)*(p-ppast))
        i<-i+1
      }
      s<-summary(glm(c(z)~1+c(x1)+c(x2)+c(present)+c(zav), family = binomial(link = "logit")))
      p<-s$coefficients[,1]
      beta0 = p[1]
      beta1 = p[2]
      beta2 = p[3]
      rho2 = p[5]
      x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2

      present<-voisin%*%(z - (exp(x_beta+(rho2)*zav )/(1+exp(x_beta+(rho2)*zav ))))
      LPL<-0
      for (i in 1:length(z)){
        LPL<-LPL+z[i]*(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*present[i] + p[5]*zav[i]) + log(1 - (exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*present[i] + p[5]*zav[i])/(1+exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*present[i] + p[5]*zav[i]))))
      }

      PLEstim[[k]]<-list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
      matpasB[k,3:7]<-PLEstim[[k]][[1]]
      matpasB[k,8:12]<-PLEstim[[k]][[2]]
      matpasB[k,13]<-PLEstim[[k]][[4]]
    }
    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vx','vy','beta0','beta1','beta2','rho1','rho2','sdbeta0','sdbeta1','sdbeta2','sdrho1','sdrho2','LPL')
  }

  if (cov == 3){
    for (k in 1:tours){
      data <- data2
      voisin <- modvoisB[[k]]
      present <- (voisin)%*%z

      p <- glm(c(z)~1+c(x1)+c(x2)+c(x3)+c(present)+c(zav), family = binomial(link = "logit"))$coef

      ppast<-c(1,1,1,1,1,1)
      delta<-sum((p-ppast)*(p-ppast))
      i<-1

      while(delta>0.00001){
        ppast<-p
        beta0 <- p[1]
        beta1 <- p[2]
        beta2 <- p[3]
        beta3 <- p[4]
        rho2 <- p[6]
        x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2 + (beta3)*x3

        present<-voisin%*%(z - (exp(x_beta+(rho2)*zav )/(1+exp(x_beta+(rho2)*zav ))))
        p<- glm(c(z)~1+c(x1)+c(x2)+c(x3)+c(present)+c(zav), family = binomial(link = "logit"))$coef
        delta<-sum((p-ppast)*(p-ppast))
        i<-i+1
      }
      s<-summary(glm(c(z)~1+c(x1)+c(x2)+c(x3)+c(present)+c(zav), family = binomial(link = "logit")))
      p<-s$coefficients[,1]
      beta0 <- p[1]
      beta1 <- p[2]
      beta2 <- p[3]
      beta3 <- p[4]
      rho2 <- p[6]
      x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2 + (beta3)*x3

      present<-voisin%*%(z - (exp(x_beta+(rho2)*zav )/(1+exp(x_beta+(rho2)*zav ))))
      LPL<-0
      for (i in 1:length(z)){
        LPL<-LPL+z[i]*(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*x3[i] + p[5]*present[i] + p[6]*zav[i]) + log(1 - (exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*x3[i] + p[5]*present[i] + p[6]*zav[i])/(1+exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*x3[i] + p[5]*present[i] + p[6]*zav[i]))))
      }

      PLEstim[[k]]<- list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
      matpasB[k,3:8]<-PLEstim[[k]][[1]]
      matpasB[k,9:14]<-PLEstim[[k]][[2]]
      matpasB[k,15]<-PLEstim[[k]][[4]]
    }
    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vx','vy','beta0','beta1','beta2','beta3','rho1','rho2','sdbeta0','sdbeta1','sdbeta2','sdbeta3','sdrho1','sdrho2','LPL')
  }


  if (graph){
    tfin <- dim(data5)[2]
    plot <- ggplot2::ggplot(data = as.data.frame(data5), ggplot2::aes(y = data5[,2], x = data5[,1],col = data5[,tfin ])) + ggplot2::geom_point() +  ggplot2::scale_color_gradient(low="white", high="red") + ggplot2::xlab("NRang") + ggplot2::ylab("NCep") + ggplot2::guides(color= ggplot2::guide_legend(title = "Levels"))
    if (!(is.null(buil))){
      r <- mattotB[,-(1:2)]
      r <- r[rev(order(r[,"LPL"])),]
    }
    else{r <- mattotB[rev(order(mattotB[,"LPL"])),]}
    result <- list('matrix' = r,'plot' = plot)
  }
  else{
    if (!(is.null(buil))){
      result <- mattotB[,-(1:2)]
    }
    else{result <- mattotB[rev(order(mattotB[,"LPL"])),]}
  }

  remove("mattotB","matpasB","data2","PLEstim","zav","z","data5")



  }


if (pastcov){


  cov1 = covariate1
  cov2 = covariate2
  cov3 = covariate3
  vlmaxpresent = vxmaxpresent
  vcmaxpresent = vymaxpresent
  vlpresent = vxpresent
  vcpresent = vypresent
  vlmaxpast = vxmaxpast
  vcmaxpast = vymaxpast
  vlpast = vxpast
  vcpast = vypast
  dl = dx
  dc = dy


  if (!(norm1 %in% c("euclidean","inf","abs","lin"))){stop("Norm must be \"euclidean\" \"inf\" \"abs\" or \"lin\" ")}
  if (!(norm2 %in% c("euclidean","inf","abs","lin"))){stop("Norm must be \"euclidean\" \"inf\" \"abs\" or \"lin\" ")}


  if (!(graph %in% c(TRUE,FALSE))){stop("graph must be TRUE or FALSE")}
  if (!(swpresent %in% c(TRUE,FALSE))){stop("swpresent must be TRUE or FALSE")}
  if (!(swpast %in% c(TRUE,FALSE))){stop("swpast must be TRUE or FALSE")}

  if (!("ggplot2" %in% rownames(installed.packages())) && graph){stop("ggplot2 package must be installed before.")}

  if ((dl < 0 )){stop("dl must be positive")}
  if ((dc < 0)){stop("dc must be positive")}

  if (swpresent && (vlmaxpresent < 0)){stop("vlmaxpresent must be positive")}
  if (swpresent && (vcmaxpresent < 0)){stop("vcmaxpresent must be positive")}
  if (!swpresent && (vlpresent < 0)){stop("vlpresent must be positive")}
  if (!swpresent && (vcpresent < 0)){stop("vcpresent must be positive")}

  if (swpast && (vlmaxpast < 0)){stop("vlmaxpast must be positive")}
  if (swpast && (vcmaxpast < 0)){stop("vcmaxpast must be positive")}
  if (!swpast && (vlpast < 0)){stop("vlpast must be positive")}
  if (!swpast && (vcpast < 0)){stop("vcpast must be positive")}



  data5 <- as.matrix(data)
  data2<-data5[,-(1:2)]

  if (is.null(cov1) && is.null(cov2) && is.null(cov3)){  #0 covariable
    cov = 0
  }

  if (!is.null(cov1) && is.null(cov2) && is.null(cov3)){  #1 covariable
    cov = 1
    x1 <- as.matrix(cov1)
    if (!((dim(x1)[1] == dim(data2)[1]) && (dim(x1)[2] == (dim(data2)[2]-1)))){stop("Dimensions of covariate1 does not match the dimensions of the dataset.")}
  }

  if (!is.null(cov1) && !is.null(cov2) && is.null(cov3)){  #2 covariable
    cov = 2
    x1 <- as.matrix(cov1)
    x2 <- as.matrix(cov2)
    if (!((dim(x1)[1] == dim(data2)[1]) && (dim(x1)[2] == (dim(data2)[2]-1)))){stop("Dimensions of covariate1 does not match the dimensions of the dataset.")}
    if (!((dim(x2)[1] == dim(data2)[1]) && (dim(x2)[2] == (dim(data2)[2]-1)))){stop("Dimensions of covariate2 does not match the dimensions of the dataset.")}
  }

  if (!is.null(cov1) && !is.null(cov2) && !is.null(cov3)){  #3 covariable
    cov = 3
    x1 <- as.matrix(cov1)
    x2 <- as.matrix(cov2)
    x3 <- as.matrix(cov3)
    if (!((dim(x1)[1] == dim(data2)[1]) && (dim(x1)[2] == (dim(data2)[2]-1)))) {stop("Dimensions of covariate1 does not match the dimensions of the dataset.")}
    if (!((dim(x2)[1] == dim(data2)[1]) && (dim(x2)[2] == (dim(data2)[2]-1)))) {stop("Dimensions of covariate2 does not match the dimensions of the dataset.")}
    if (!((dim(x3)[1] == dim(data2)[1]) && (dim(x3)[2] == (dim(data2)[2]-1)))) {stop("Dimensions of covariate3 does not match the dimensions of the dataset.")}
  }



  if (is.null(buildpres)){

    if (swpresent){
      tours1 <- vlmaxpresent*vcmaxpresent
      modvois1<-vector("list", vlmaxpresent*vcmaxpresent)
      k<-1
      for (vl1 in (1:vlmaxpresent)){
        for (vc1 in (1:vcmaxpresent)){
          modvois1[[k]]<-build(data = data5,vx = vl1,vy = vc1,norm = norm1)
          k<-k+1
        }
      }
    }
    else{
      tours1 <- 1
      modvois1<-vector("list", 1)
      modvois1[[1]]<-build(data = data5,vx = vlpresent,vy = vcpresent,norm = norm1)
    }
  }
  else {
    tours1 <- 1
    modvois1<-vector("list", 1)
    modvois1[[1]]<-buildpres
    swpresent = FALSE
  }




  ##############################################
  ############ Voisin 2 ########################
  ##############################################


  if (is.null(buildpast)){

    if (swpast){
      tours2 <- vlmaxpast*vcmaxpast
      modvois2<-vector("list", vlmaxpast*vcmaxpast)
      k<-1
      for (vl1 in (1:vlmaxpast)){
        for (vc1 in (1:vcmaxpast)){
          modvois2[[k]]<-build(data = data5,vx = vl1,vy = vc1,norm = norm2)
          k<-k+1
        }
      }
    }
    else{
      tours2 <- 1
      modvois2<-vector("list", 1)
      modvois2[[1]]<-build(data = data5,vx = vlpast,vy = vcpast,norm = norm2)
    }
  }
  else {
    tours2 <- 1
    modvois2<-vector("list", 1)
    modvois2[[1]]<-buildpast
    swpast = FALSE
  }

  ##########################################
  ########## Estimation ####################
  ##########################################



  tp<-ncol(data2)-1
  mattotB<-NULL

  PLEstim<-vector("list",tours1*tours2)
  rw <- 13 + (2*cov)
  tor <- tours1*tours2*rw
  matpasB<-matrix(rep(0,tor),ncol=rw)

  #balayages

  if (swpresent){

    if (swpast){
      k = 1
      for (v1 in (1:vlmaxpresent)){
        for (v2 in (1:vcmaxpresent)){
          for (v3 in (1:vlmaxpast)){
            for (v4 in (1:vcmaxpast)){

              matpasB[k,1]<-v1
              matpasB[k,2]<-v2
              matpasB[k,3]<-v3
              matpasB[k,4]<-v4

              k = k+1
            }
          }
        }
      }
    }

    if (!swpast){

      k = 1
      for (v1 in (1:vlmaxpresent)){
        for (v2 in (1:vcmaxpresent)){

          matpasB[k,1]<-v1
          matpasB[k,2]<-v2

          if (is.null(buildpast)){
            matpasB[k,3]<-vlpast
            matpasB[k,4]<-vcpast
          }
          else{
            matpasB[k,3]<-"custom"
            matpasB[k,4]<-"custom"
          }

          k = k+1

        }
      }
    }
  }

  if (!swpresent){
    if(swpast){

      k = 1
      for (v1 in (1:vlmaxpast)){
        for (v2 in (1:vcmaxpast)){

          if (is.null(buildpres)){
            matpasB[k,1]<-vlpresent
            matpasB[k,2]<-vcpresent
          }
          else{
            matpasB[k,1]<-"custom"
            matpasB[k,2]<-"custom"
          }
          matpasB[k,3]<-v1
          matpasB[k,4]<-v2

          k = k+1

        }
      }

    }
    if(!swpast){

      if (is.null(buildpres)){
        matpasB[,1]<-vlpresent
        matpasB[,2]<-vcpresent
      }
      else{
        matpasB[,1]<-"custom"
        matpasB[,2]<-"custom"
      }
      if (is.null(buildpast)){
        matpasB[,3]<-vlpast
        matpasB[,4]<-vcpast
      }
      else{
        matpasB[,3]<-"custom"
        matpasB[,4]<-"custom"
      }

    }
  }


  #Calcul final


  #covariables


  zav<-as.matrix(data2[,-(tp+1)])
  z<-as.matrix(data2[,-1])



  if (cov == 0){
    k = 0
    for (k1 in 1:tours1){
      for (k2 in 1:tours2){
        k = k+1
        data = data2
        voisin = modvois1[[k1]]
        voisin2 = modvois2[[k2]]

        present <- voisin%*%z
        past <- voisin2%*%zav
        p <- glm(c(z)~1 + c(past) + c(present)+c(zav), family = binomial(link = "logit"))$coef

        ppast<-c(1,1,1,1)
        delta<-sum((p-ppast)*(p-ppast))
        i<-1

        while(delta>0.00001){
          ppast<-p
          beta0 = p[1]
          betapast = p[2]
          rho2 = p[4]
          x_beta<- (beta0) + (betapast)*past

          present<-voisin%*%(z - (exp(x_beta + (rho2)*zav )/(1+exp(x_beta + (rho2)*zav ))))
          p<- glm(c(z) ~ 1 + c(past) + c(present) + c(zav), family = binomial(link = "logit"))$coef
          delta<-sum((p-ppast)*(p-ppast))
          i<-i+1
        }
        s<-summary(glm(c(z) ~ 1 + c(past) + c(present) + c(zav), family = binomial(link = "logit")))
        p<-s$coefficients[,1]
        beta0 = p[1]
        betapast = p[2]
        rho2 = p[4]
        x_beta<- (beta0) + (betapast)*past

        present<-voisin%*%(z - (exp(x_beta + (rho2)*zav )/(1+exp(x_beta + (rho2)*zav ))))
        LPL<-0
        for (i in 1:length(z)){
          LPL<-LPL+z[i]*(p[1] + p[2]*past[i] + p[3]*present[i] + p[4]*zav[i]) + log(1 - (exp(p[1] + p[2]*past[i] + p[3]*present[i] + p[4]*zav[i])/(1+exp(p[1] + p[2]*past[i] + p[3]*present[i] + p[4]*zav[i]))))
        }

        PLEstim[[k]] <- list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
        matpasB[k,5:8] <- PLEstim[[k]][[1]]
        matpasB[k,9:12] <- PLEstim[[k]][[2]]
        matpasB[k,13] <- PLEstim[[k]][[4]]
      }
    }

    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vlpast','vcpast','vlpresent','vcpresent','beta0','betapast','rho1','rho2','sdbeta0','sdbetapast','sdrho1','sdrho2','LPL')
  }

  if (cov == 1){
    k = 0
    for (k1 in 1:tours1){
      for (k2 in 1:tours2){
        k = k+1
        data <- data2
        voisin <- modvois1[[k1]]
        voisin2 <- modvois2[[k2]]

        present<-voisin%*%z
        past <- voisin2%*%zav
        p <- glm(c(z)~1 + c(x1) + c(past) + c(present)+c(zav), family = binomial(link = "logit"))$coef

        ppast <- c(1,1,1,1,1)
        delta <- sum((p-ppast)*(p-ppast))
        i <- 1

        while(delta>0.00001){
          ppast <- p
          beta0 <- p[1]
          beta1 <- p[2]
          betapast <- p[3]
          rho2 <- p[5]
          x_beta <- (beta0) + (beta1)*x1 + (betapast)*past

          present <- voisin%*%(z - (exp(x_beta + (rho2)*zav)/(1+exp(x_beta + (rho2)*zav))))
          p <- glm(c(z) ~ 1 + c(x1) + c(past) + c(present) + c(zav), family = binomial(link = "logit"))$coef
          delta <- sum((p-ppast)*(p-ppast))
          i <- i+1
        }
        s <- summary(glm(c(z) ~ 1 + c(x1) + c(past) + c(present) + c(zav), family = binomial(link = "logit")))
        p <- s$coefficients[,1]
        beta0 <- p[1]
        beta1 <- p[2]
        betapast <- p[3]
        rho2 <- p[5]
        x_beta <- (beta0) + (beta1)*x1 + (betapast)*past

        present <- voisin%*%(z - (exp(x_beta + (rho2)*zav)/(1+exp(x_beta + (rho2)*zav))))
        LPL<-0
        for (i in 1:length(z)){
          LPL <- LPL+z[i]*(p[1] + p[2]*x1[i] + p[3]*past[i] + p[4]*present[i] + p[5]*zav[i]) + log(1 - (exp(p[1] + + p[2]*x1[i] + p[3]*past[i] + p[4]*present[i] + p[5]*zav[i])/(1+exp(p[1] + + p[2]*x1[i] + p[3]*past[i] + p[4]*present[i] + p[5]*zav[i]))))
        }

        PLEstim[[k]]<-list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
        matpasB[k,5:9]<-PLEstim[[k]][[1]]
        matpasB[k,10:14]<-PLEstim[[k]][[2]]
        matpasB[k,15]<-PLEstim[[k]][[4]]
      }
    }

    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vlpast','vcpast','vlpresent','vcpresent','beta0','beta1','betapast','rho1','rho2','sdbeta0','sdbeta1','sdbetapast','sdrho1','sdrho2','LPL')
  }


  if (cov == 2){
    k = 0
    for (k1 in 1:tours1){
      for (k2 in 1:tours2){
        k = k+1
        data <- data2
        voisin <- modvois1[[k1]]
        voisin2 <- modvois2[[k2]]

        present<-voisin%*%z
        past <- voisin2%*%zav

        p <- param<- glm(c(z)~1 + c(x1) + c(x2) + c(past) + c(present)+c(zav), family = binomial(link = "logit"))$coef
        ppast <- c(1,1,1,1,1,1)
        delta <- sum((p-ppast)*(p-ppast))
        i <- 1
        while(delta>0.00001){
          ppast<-p
          beta0 <- p[1]
          beta1 <- p[2]
          beta2 <- p[3]
          betapast <- p[4]
          rho2 <- p[6]
          x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2 + (betapast)*past

          present<-voisin%*%(z - (exp(x_beta + (rho2)*zav )/(1+exp(x_beta + (rho2)*zav ))))
          p<- glm(c(z) ~ 1 + c(x1) + c(x2) + c(past) + c(present) + c(zav), family = binomial(link = "logit"))$coef
          delta<-sum((p-ppast)*(p-ppast))
          i<-i+1
        }
        s<-summary(glm(c(z) ~ 1 + c(x1) + c(x2) + c(past) + c(present) + c(zav), family = binomial(link = "logit")))
        p<-s$coefficients[,1]
        beta0 <- p[1]
        beta1 <- p[2]
        beta2 <- p[3]
        betapast <- p[4]
        rho2 <- p[6]
        x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2 + (betapast)*past

        present<-voisin%*%(z - (exp(x_beta + (rho2)*zav )/(1+exp(x_beta + (rho2)*zav ))))
        LPL<-0
        for (i in 1:length(z)){
          LPL<-LPL+z[i]*(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*past[i] + p[5]*present[i] + p[6]*zav[i]) + log(1-(exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*past[i] + p[5]*present[i] + p[6]*zav[i])/(1+exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*past[i] + p[5]*present[i] + p[6]*zav[i]))))
        }
        PLEstim[[k]] <- list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
        matpasB[k,5:10]<-PLEstim[[k]][[1]]
        matpasB[k,11:16]<-PLEstim[[k]][[2]]
        matpasB[k,17]<-PLEstim[[k]][[4]]
      }
    }

    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vlpast','vcpast','vlpresent','vcpresent','beta0','beta1','beta2','betapast','rho1','rho2','sdbeta0','sdbeta1','sdbeta2','sdbetapast','sdrho1','sdrho2','LPL')
  }



  if (cov == 3){
    k = 0
    for (k1 in 1:tours1){
      for (k2 in 1:tours2){
        k = k+1
        data <- data2
        voisin <- modvois1[[k1]]
        voisin2 <- modvois2[[k2]]

        present<-voisin%*%z
        past <- voisin2%*%zav
        p <- glm(c(z)~1 + c(x1) + c(x2) + c(x3) + c(past) + c(present)+c(zav), family = binomial(link = "logit"))$coef
        ppast <- c(1,1,1,1,1,1,1)
        delta <- sum((p-ppast)*(p-ppast))
        i<-1
        while(delta>0.00001){
          ppast<-p
          beta0 <- p[1]
          beta1 <- p[2]
          beta2 <- p[3]
          beta3 <- p[4]
          betapast <- p[5]
          rho2 <- p[7]
          x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2 + (beta3)*x3 + (betapast)*past

          present<-voisin%*%(z - (exp(x_beta + (rho2)*zav )/(1+exp(x_beta + (rho2)*zav ))))
          p<- glm(c(z) ~ 1 + c(x1) + c(x2) + c(x3) + c(past) + c(present) + c(zav), family = binomial(link = "logit"))$coef
          delta<-sum((p-ppast)*(p-ppast))
          i<-i+1
        }
        s <- summary(glm(c(z) ~ 1 + c(x1) + c(x2) + c(x3) + c(past) + c(present) + c(zav), family = binomial(link = "logit")))
        p<-s$coefficients[,1]
        beta0 <- p[1]
        beta1 <- p[2]
        beta2 <- p[3]
        beta3 <- p[4]
        betapast <- p[5]
        rho2 <- p[7]
        x_beta<- (beta0) + (beta1)*x1 + (beta2)*x2 + (beta3)*x3 + (betapast)*past

        present<-voisin%*%(z - (exp(x_beta + (rho2)*zav )/(1+exp(x_beta + (rho2)*zav ))))
        LPL<-0
        for (i in 1:length(z)){
          LPL<-LPL+z[i]*(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*x3[i] + p[5]*past[i] + p[6]*present[i] + p[7]*zav[i]) + log(1 - (exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*x3[i] + p[5]*past[i] + p[6]*present[i] + p[7]*zav[i])/(1+exp(p[1] + p[2]*x1[i] + p[3]*x2[i] + p[4]*x3[i] + p[5]*past[i] + p[6]*present[i] + p[7]*zav[i]))))
        }


        PLEstim[[k]] <- list(s$coefficients[,1],s$coefficients[,2],s$aic,LPL)
        matpasB[k,5:11]<-PLEstim[[k]][[1]]
        matpasB[k,12:18]<-PLEstim[[k]][[2]]
        matpasB[k,19]<-PLEstim[[k]][[4]]
      }
    }

    mattotB<-rbind(mattotB,matpasB)
    colnames(mattotB) <- c('vlpast','vcpast','vlpresent','vcpresent','beta0','beta1','beta2','beta3','betapast','rho1','rho2','sdbeta0','sdbeta1','sdbeta2','sdbeta3','sdbetapast','sdrho1','sdrho2','LPL')
  }



  if (graph){
    tfin <- dim(data5)[2]
    x_present = data5[,1]
    y_present = data5[,2]
    x_past = data5[,1]
    y_past = data5[,2]
    plotpresent <- ggplot2::ggplot(data = as.data.frame(data5), ggplot2::aes(y = y_present, x = x_present,col = data5[,tfin ])) + ggplot2::geom_point() +  ggplot2::scale_color_gradient(low="white", high="red") + ggplot2::xlab("NRang") + ggplot2::ylab("NCep") + ggplot2::guides(color= ggplot2::guide_legend(title = "Levels"))
    plotpast <- ggplot2::ggplot(data = as.data.frame(data5), ggplot2::aes(y = y_past, x = x_past,col = data5[,(tfin-1)])) + ggplot2::geom_point() +  ggplot2::scale_color_gradient(low="white", high="red") + ggplot2::xlab("NRang") + ggplot2::ylab("NCep") + ggplot2::guides(color= ggplot2::guide_legend(title = "Levels"))
    if ((!(is.null(buildpast)))&&(is.null(buildpres))){
      r <- mattotB[,-(1:2)]
      r <- r[rev(order(r[,"LPL"])),]
    }
    else if ((!(is.null(buildpres)))&&(is.null(buildpast))){
      r <- mattotB[,-(3:4)]
      r <- r[rev(order(r[,"LPL"])),]
    }
    else if ((!(is.null(buildpres)))&&(!(is.null(buildpast)))){
      r <- mattotB[,-(1:4)]
      r <- r[rev(order(r[,"LPL"])),]
    }
    else{
      r <- mattotB
      r <- r[rev(order(r[,"LPL"])),]
    }
    result <- list('matrix' = r,'plotpresent ' = plotpresent,'plotpast' = plotpast)
  }


  else{
    if ((!(is.null(buildpast)))&&(is.null(buildpres))){
      result <- mattotB[,-(1:2)]
      result <- result[rev(order(result[,"LPL"])),]
    }
    else if ((!(is.null(buildpres)))&&(is.null(buildpast))){
      result <- mattotB[,-(3:4)]
      result <- result[rev(order(result[,"LPL"])),]
    }
    else if ((!(is.null(buildpres)))&&(!(is.null(buildpast)))){
      result <- mattotB[,-(1:4)]
    }
    else{
      result <- mattotB
      result <- result[rev(order(result[,"LPL"])),]
    }
  }


  }

result
}
