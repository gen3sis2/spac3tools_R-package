## METADATA ===============================================================
## Description: just a list of cost functions
##
## R version: 4.2.2 for Windows
## Date: 2022-12-20 13:13:40
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##
# list of cost functions gcfs = gen3sis cost funtions list
gcf <- list()

# Available Cost Functions:
gcf$only_dist_m <- function(source, habitable_src, dest, habitable_dest){
  return(1)
}
gcf$only_dist_Km <- function(source, habitable_src, dest, habitable_dest){
  return(1/1000)
}
gcf$XXHarderNA_dist_m <- function(source, habitable_src, dest, habitable_dest){
  if (!all(habitable_src, habitable_dest)) {
    return(2)
  } else {
    return(1)
  }
}
gcf$XXHarderNA_dist_Km <- function(source, habitable_src, dest, habitable_dest){
  if (!all(habitable_src, habitable_dest)) {
    return(2/1000)
  } else {
    return(1/1000)
  }
}
gcf$PhysDiv_Water_Km <- function(source, habitable_src, dest, habitable_dest) {
  if (!all(habitable_src, habitable_dest)) {
    return(2/1000)
  } else {
    return((max(source["physdiv"], dest["physdiv"])+1/1000))
  }
}


## Example from AS on equal area projection. Attention: This seams to no work with the
## spherical corrections of the borders.
## equal area projection:
# # "I just use the projectRaster function to change the projection.
# # The geometry will always be distorted when you plot a globe onto a grid,
# # the different projections will just distort it in different ways.
# # To me, they are all equally valid."
# projectRaster(ras, res=100*1000, crs = crs('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +no_defs'))
