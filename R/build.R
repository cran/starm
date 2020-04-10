#'Creation of the neighbourhood's matrix.
#'@param data dataset with first column the X-coordinates of the sites and the second the Y-coodinates of the sites.
#'@param vx integer, first parameter of the neighbourhood  ( i.e. first parameter of ellipse if \code{norm = "euclidean"}  for instance). \code{vx = 3} by default.
#'@param vy integer, second parameter of the neighbourhood  ( i.e. second parameter of ellipse if \code{norm = "euclidean"}  for instance). \code{vy = 3} by default.
#'@param dx positive real, distance between sites on a row. \code{dx = 1} by default.
#'@param dy positive real, distance between sites on a column. \code{dy = 1} by default.
#'@param t double. If \code{selec = TRUE}, each neighborhood will contain only elements which the type is in \code{t} (see \code{examples}).
#'@param norm Response type : "euclidean" "inf" "abs" "lin".  \code{norm = "euclidean"} by default.
#'@param returnplot If \code{TRUE}, will return the plot of the most recent neighborhhod in addition to the neighborhood matrix.
#'@param selec see \code{t}.
#'
#'@return The neighborhood matrix
#'
#'
#'
#'@export
#'
#'@details The function will return the neighborhood matrix of a dataset which must contain coodinates in the two first columns and a third column at least with the "type" of each site (it can be only "0" or "1" for example). The parameter \code{norm} let you choose between 4 sorts of neighborhood : 3 ellipses in norm 1, 2 or infinite (resp "abs","euclidean" and "inf") with the parameters \code{vx} and \code{vy} which are the width and the height of the ellipse, and the norm \code{lin} will condider only sites on the same row and column with the same parameters \code{vx} and \code{vy}.
#'
#'@value The neighborhood matrix. It is a matrix of  dimension $n$ the number of lines in the dataset. It is an adjacency matrix that contains 1 at position (i,j) if i and j are neighbors and 0 if not.
#'
#'@note If \code{returnplot = TRUE}, \code{variable$plot} will return an exemple of the choosen neighborhood on a center point of the dataset.
#'
#'
#'@examples
#'
#'data <- plantillness
#'v <- which((data$NRang <= 20))
#'data <- data[v,]
#'v <- which(data$NCep <= 20)
#'data<-data[v,]
#'res <- build(data = data)
#'
#'
#'\donttest{
#'#Example with the plantillness dataset and the plot available :
#'
#'res <- build(data = plantillness,returnplot = TRUE,vx = 5,vy = 5)
#'
#'
#'
#'#Example with the plantillness dataset, only considering the sites of the type "0" :
#'
#'res <- build(data = plantillness, selec = TRUE, t = c(0),vx = 5,vy = 7,norm = "inf")
#'
#'}


build <- function(data = 0,vx = 3,vy = 3,dx = 1,dy = 1,selec = FALSE,t = 0,norm = "euclidean",returnplot = FALSE){

  if (!(returnplot %in% c(TRUE,FALSE))){stop("graph must be TRUE or FALSE")}
  if (!("ggplot2" %in% rownames(installed.packages())) && returnplot){stop("ggplot2 package must be installed before.")}
  if (!(norm %in% c("euclidean","inf","abs","lin"))){stop("Norm must be \"euclidean\" \"inf\" \"abs\" or \"lin\".")}


  if (selec == FALSE){
    temp <- unique(data[,dim(data)[2]])
  }
  else{temp = t}

  if (norm == "lin"){

    n <- max(data[,1])
    m <- max(data[,2])

    A1 = matrix(rep(-1,n*m),nrow = m,ncol = n)

    v1 = vx/dx
    v2 = vy/dy
    A = matrix(0, m * n, m * n)
    for (i in 1:(m * n))
    {
      if (v1 >= 1){
        left = i - 1
        right = i + 1
        leftbol = FALSE
        rightbol = FALSE
        if (left %% n != 0)
          A[i, left] = 1
        else leftbol = TRUE
        if (i %% n != 0)
          A[i, right] = 1
        else rightbol = TRUE
        if (v1>1)
        {
          for   (j in 2:v1)
          {
            left = i - j
            right = i + j
            if (left%%n == 0 | leftbol)
              leftbol = TRUE
            else
              A[i, left] = 1
            if (right%%n == 1 | rightbol)
              rightbol = TRUE
            else
              A[i, right] = 1
          }
        }
      }
      if (v2 >= 1){
        for (j in 1:v2)
        {
          up = i - j*n
          down = i + j*n
          if (up > 0)
            A[i, up] = 1
          if (down <= (m * n))
            A[i, down] = 1
        }
      }
    }

    for (i in (1:(dim(data)[1]))){
      A1[data[i,2],data[i,1]] = (data[i,1]-1)*m + data[i,2]
    }

    x<- A1[,1]

    for (i in (2:n)){
      x <- c(x,A1[,i])
    }

    y <- which(x >0)

    A <- A[y,y]

    for (i in 1:dim(A)[1]){
      A[i,] = (A[i,]*(data[,dim(data)[2]]%in%temp))}

    B = A

  }

  else{
    l = data[,1]
    c = data[,2]
    v1 = vx/dx
    v2 = vy/dy
    n<-length(l)
    B = matrix(0,n,n)
    for (i in (1:n)){
      voisl<-l%in%((l[i]-v1):(l[i]+v1))
      voisr<-c%in%((c[i]-v2):(c[i]+v2))
      A<-(which((voisl*voisr)<=1))
      A<-A[which(A>i)]
      for (j in A){
        if (norm == "euclidean"){
          d<-((l[i]-l[j])/v1)^2 + ((c[i]-c[j])/v2)^2
        }

        if (norm == "abs"){
          d<-(abs(l[i]-l[j])/v1) + (abs(c[i]-c[j])/v2)
        }

        if (norm == "inf"){
          d<-max((abs(l[i]-l[j])/v1),(abs(c[i]-c[j])/v2))
        }

        if (d <= 1 && data[j,dim(data)[2]]%in%temp){
          B[i,j]<-1
        }
        if (d <= 1 && data[i,dim(data)[2]]%in%temp){
          B[j,i]<-1
        }
      }
    }

    if (returnplot){
      n <- max(data[,1])
      m <- max(data[,2])
      n1 <- which(data[,1] == floor(n/2))
      m1 <- which(data[,2] == floor(m/2))
      f <- intersect(n1,m1)
      a <- B[f,]
      plot <- ggplot2::ggplot(data = data, ggplot2::aes(y = data[,2], x = data[,1],col = a)) + ggplot2::geom_point() + ggplot2::xlab("NRang") + ggplot2::ylab("NCep")
      result <- list('matrix' = B,'plot' = plot)
    }
    else{result <- B}

    return(result)
  }
}
