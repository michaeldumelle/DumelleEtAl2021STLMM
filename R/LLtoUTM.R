#' Transverse Mercator Projection
#'
#' Project latitude and longitude spatial coordinates using a Transverse Mercator
#' projection. Units are in kilometers, the minimum projected coordinates in
#' the x and y direction are 0, and projected coordinates are scaled to 1 kilometer.
#'
#' @param cm iCentral meridian. A common choice is the
#' mean of the longitude values in your data set
#'
#' @param lat Vector of latitudes.
#'
#' @param lon Vector of longitudes.
#'
#' @param xcol Name of the transformed x coordinates output by the function.
#'
#' @param ycol Name of the transformed y coordinates output by the function.
#'
#' @param minx A minimum x coordinate to center by. Defaults to setting the minimum
#'   transformed coordinate to 0.
#'
#' @param miny A minimum y coordinate to center by. Defaults to setting the minimum
#'   transformed coordinate to 0.
#'
#' @return A list returning the transformed coordinates, central meridian argument,
#'   minimum x coordinate, and minimum y coordinate.
#'
#' @export
LLtoUTM <- function(cm, lat, lon, xcol = "x", ycol = "y", minx = NULL, miny = NULL)
{
  # check if any longitude values straddle the -180, +180 longitude line
  # if so, convert minus longitude values to longitude values > 180
  if(any(lon > 90 & lon < 180) & any(lon > -180 & lon < -90))
    lon[lon < 0] <- 360 + lon[lon < 0]

  # initialize some variables
  e2 <- 0.00676865799729
  a <- 6378206.4
  ep2 <- e2 / (1-e2)
  drc <- pi / 180
  sc <- 0.9996
  fe <- 500000
  ftm <- 0.30480371
  #calculate some frequently used values
  lar <- lat * drc
  ls <- sin(lar)
  ls2 <- ls^2
  els2 <- ep2 * ls2
  lc <- cos(lar)
  lc2 <- lc^2
  lc3 <- lc^3
  lc5 <- lc^5
  elc2 <- ep2 * lc2
  lt2 <- tan(lar)^2
  lt4 <- lt2^2
  # do the transformation
  v <- a/sqrt(1 - e2*ls2)
  p <- drc*(cm - lon)
  temp <- 5104.57388 - (lc2*(21.73607 - 0.11422*lc2))
  r1 <- 6367399.689*(lar - ls*lc*0.000001*temp)
  r2 <- (v*ls*lc*p^2)/2
  temp <- 5 - lt2 + 9*elc2 + (2*elc2)^2
  r3 <- (v*ls*lc3*p^4*temp)/24
  r4 <- v*lc*p
  temp <- 1 - lt2 + elc2
  r5 <- (v*lc3*p^3*temp)/6
  temp <- 61 - 58*lt2 + lt4 + 270*elc2 - 330*els2
  ra6 <- (v*ls*lc5*p^6*temp)/720
  temp <- 5 - 18*lt2 + lt4 + 14*elc2 - 58*els2
  rb5 <- (v*lc5*p^5*temp)/120
  northing <- sc*(r1 + r2 + r3 + ra6)
  easting <- -sc*(r4 + r5 + rb5)
  if(is.null(miny)) miny <- min(northing)
  y <- (northing - miny)/1000
  if(is.null(minx)) minx <- min(easting)
  x <- (easting - minx)/1000

  out <- cbind(x,y)
  colnames(out) <- c(xcol, ycol)
  list(xy = out, cm = cm, minx = minx, miny = miny)
}
