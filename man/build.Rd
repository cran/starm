% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build.R
\name{build}
\alias{build}
\title{Creation of the neighbourhood's matrix.}
\usage{
build(data = 0, vx = 3, vy = 3, dx = 1, dy = 1, selec = FALSE,
  t = 0, norm = "euclidean", returnplot = FALSE)
}
\arguments{
\item{data}{dataset with first column the X-coordinates of the sites and the second the Y-coodinates of the sites.}

\item{vx}{integer, first parameter of the neighbourhood  ( i.e. first parameter of ellipse if \code{norm = "euclidean"}  for instance). \code{vx = 3} by default.}

\item{vy}{integer, second parameter of the neighbourhood  ( i.e. second parameter of ellipse if \code{norm = "euclidean"}  for instance). \code{vy = 3} by default.}

\item{dx}{positive real, distance between sites on a row. \code{dx = 1} by default.}

\item{dy}{positive real, distance between sites on a column. \code{dy = 1} by default.}

\item{selec}{see \code{t}.}

\item{t}{double. If \code{selec = TRUE}, each neighborhood will contain only elements which the type is in \code{t} (see \code{examples}).}

\item{norm}{Response type : "euclidean" "inf" "abs" "lin".  \code{norm = "euclidean"} by default.}

\item{returnplot}{If \code{TRUE}, will return the plot of the most recent neighborhhod in addition to the neighborhood matrix.}
}
\value{
The neighborhood matrix
}
\description{
Creation of the neighbourhood's matrix.
}
\details{
The function will return the neighborhood matrix of a dataset which must contain coodinates in the two first columns and a third column at least with the "type" of each site (it can be only "0" or "1" for example). The parameter \code{norm} let you choose between 4 sorts of neighborhood : 3 ellipses in norm 1, 2 or infinite (resp "abs","euclidean" and "inf") with the parameters \code{vx} and \code{vy} which are the width and the height of the ellipse, and the norm \code{lin} will condider only sites on the same row and column with the same parameters \code{vx} and \code{vy}.
}
\note{
If \code{returnplot = TRUE}, \code{variable$plot} will return an exemple of the choosen neighborhood on a center point of the dataset.
}
\examples{

data <- plantillness
v <- which((data$NRang <= 20))
data <- data[v,]
v <- which(data$NCep <= 20)
data<-data[v,]
res <- build(data = data)


\donttest{
#Example with the plantillness dataset and the plot available :

res <- build(data = plantillness,returnplot = TRUE,vx = 5,vy = 5)



#Example with the plantillness dataset, only considering the sites of the type "0" :

res <- build(data = plantillness, selec = TRUE, t = c(0),vx = 5,vy = 7,norm = "inf")

}
}
