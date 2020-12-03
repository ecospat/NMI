mve <- function (x, m.cov = covMcd(x), thresh=0.95, robust=FALSE) {
  
  ellips <- function(loc, cov) {
    dist <- sqrt(qchisq(thresh, 2)) 
    A <- solve(cov)
    eA <- eigen(A)
    ev <- eA$values
    lambda1 <- max(ev)
    lambda2 <- min(ev)
    eigvect <- eA$vectors[, order(ev)[2]]
    z <- seq(0, 2 * pi, by = 0.01)
    z1 <- dist/sqrt(lambda1) * cos(z)
    z2 <- dist/sqrt(lambda2) * sin(z)
    alfa <- atan(eigvect[2]/eigvect[1])
    r <- matrix(c(cos(alfa), -sin(alfa), sin(alfa), cos(alfa)), 
                ncol = 2)
    t(loc + t(cbind(z1, z2) %*% r))
  }
  
  if (is.data.frame(x)) 
    x <- data.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) 
    stop("x is not a numeric dataframe or matrix.")
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (p != 2) 
    stop("Dimension {= ncol(x)} must be 2!")
  if (!is.numeric(m.cov$center) || !is.numeric(m.cov$cov)) 
    stop("argument 'm.cov' must have numeric components 'center' and 'cov'")
  x.loc <- m.cov$center
  x.cov <- n/(n - 1) * m.cov$cov
  xM <- colMeans(x)
  z1 <- ellips(loc = xM, cov = n/(n - 1) * cov.wt(x)$cov)
  z2 <- ellips(loc = x.loc, cov = x.cov)
  if(robust) out <- z2 else out <- z1
  out.poly <- SpatialPolygons(list(Polygons(list(Polygon(out)), ID=1))) 
  return(out.poly)
}