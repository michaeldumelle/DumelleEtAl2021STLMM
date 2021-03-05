#' Oregon Maximum Tempertaure Data
#'
#' A dataset containing maxiumum temperature data in Oregon duing July 2019. A
#' link to NOAA's Climate Data Online (CDO) tool used to obtain the data is
#' also provided
#'
#' @format A data frame with 2972 rows and 13 variables:
#' \describe{
#'   \item{STATION}{A unique station identifier.}
#'   \item{NAME}{Name of the station.}
#'   \item{LATITUDE}{Transformed latitude coordinates (using a Transverse)
#'     Mercator Projection with a central meridian of -121.5158.}
#'   \item{LONGITUDE}{Transformed longitude coordinates (using a Transverse)
#'     Mercator Projection with a central meridian of -121.5158.}
#'   \item{ELEVATION}{Elevation in meters above mean sea level.}
#'   \item{DATE}{The date.}
#'   \item{PRCP}{Precipitation in millimeters.}
#'   \item{TMAX}{Daily maximum temperature in Farenheit.}
#'   \item{MONTH}{Month of the year.}
#'   \item{DAY}{Day of the month.}
#'   \item{TIMES}{Temporal index for unique days.}
#'   \item{SPINDEX}{Spatial index for unique sites.}
#'   \item{TYPE}{Type identifier. Equal to \code{"TRAIN"} for the training data
#'     or \code{"TEST"} for the test data.}
#'   ...
#' }
#' @source \url{https://www.ncdc.noaa.gov/cdo-web/#t=secondTabLink}
"or_data"
