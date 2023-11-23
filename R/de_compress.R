## METADATA ===============================================================
## Description: functions to recreate (aka decompress) gen3sis input only from
## gen3sis former inputs
##
## R version: 4.2.2 for Windows
## Date: 2022-12-20 13:00:20
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##


#' decompress a compressed gen3sis_space
#'
#' @param dir_input the input directory containing the gen3sis space.rds
#' and ideally the METADATA.txt
#' @param cost_function a constant function to be used Normaly
#' Default is the declared cost_functions list at cost_lists
#' Note that different const_functions can be used and it's computation
#' depends only on gen3sis::create_input_landscape function
#' @param dir_output output directory to save output. Default: dir_output=NULL
#'  is to creates a temporary directory and print it.
#' @param remove_temp_rasters boolean. If FALSE, delete temporary folder and files,
#' i.e. temporary raster brinks inside directory temp_rasters
#'
#' @return passing confirmation statement
#' @export
#'
#' @examples inst/examples/run_recreate_space.R
#' # example

decompress_space <- function(dir_input=NULL,
                             cost_function=gcfl$only_dist_Km,
                             dir_output=NULL,
                             remove_temp_rasters=TRUE
                             #, timestoMa=6 # TODO
                             # SOLVE TEMPORAL ISSUE WITH FILE STRUCTUR
                             # PROB WITH TEMPORAL AND spatial resolution
){

  if (is.null(dir_output)){
    # error("Please provide an output directory for the recreated landscape")
    dir_output <- tempdir()
  }

  cat(paste0("Using output directory: [", dir_output, "]"))

  space_file_loc <- file.path(dir_input,"landscapes.rds")
  if (file.exists(space_file_loc)){
    print(paste0("[OK] landscapes.rds found: [", space_file_loc, "]"))
  } else{
    stop(paste0("[MISSING] landscapes.rds was not found: [", space_file_loc, "]"))
  }

  gen3sis_space=readRDS(space_file_loc)
  ls <- gen3sis_space

  if (is.null(gen3sis_space)){
    error("Please provide a loaded gen3sis landscape.rds as gen3sis_space")
  }

  # create temp dir
  dir_temp_raster=file.path(dir_output, "temp_rasters")

  # create raster bricks
  # browser() FINISH THIS PIPELINE LATER

  # od <- dir_output
  dir.create(dir_output, showWarnings = FALSE)
  # create temp dir for temp rasters
  dir.create(dir_temp_raster, showWarnings = FALSE)

  # create temp rasters
  for (var_i in names(ls)){
    # var_i <- names(ls)[1]
    rstack <- NULL
    for (c_i in 3:ncol(ls[[var_i]])){
      rstack[[c_i-2]] <- rasterFromXYZ(ls[[var_i]][,c(1,2,c_i)])
    }
    rstack <- brick(rstack)
    writeRaster(rstack, filename = file.path(dir_temp_raster, paste0(var_i, ".grd")), overwrite=TRUE)
  }
  print(paste("Temporary Raster Bricks Saved to: ", dir_temp_raster))

  # load raster bricks
  #list all temp raster bricks
  bf <- list.files(dir_temp_raster, pattern=".grd")
  print(paste("Raster Bricks are: ", paste(bf, collapse = "; ")))
  #load all temp raster bricks
  b <- NULL
  for (i in 1:length(bf)){
    b[[i]] <- brick(file.path(dir_temp_raster, bf[i]))
  }
  #prepare list
  lsn <- lapply(ls, function(x){x <- NULL})
  #attribute to list
  for (i in 1:length(lsn)){
    for (il in 1:nlayers(b[[i]])) {
      lsn[[i]] <- c(lsn[[i]], b[[i]][[il]])
    }
  }
  # TODO WHY ARE THE VALUES CHANGING BEFORE AND AFTER SAVING? Check this../. between lines 35 and 38!
  # myt <- (as.numeric(colnames(ls[[1]][-c(1,2)])))/timestoMa
  # string_time_step <- paste0(formatC(round(myt,2), width=5, flag="0", digits=2, format="f"),"Ma" )
  create_input_landscape(landscapes = lsn, cost_function = cost_function, directions=8,
                         output_directory = file.path(dir_output, basename(dir_input)), #timesteps = string_time_step,
                         calculate_full_distance_matrices = T, crs=s_wgs84, verbose=T)

  #updating the landscape.rds to avoid raster FUCK-UPS
  saveRDS(ls, file.path(dir_output, basename(dir_input), "landscapes.rds"))
  print(paste0("landscapes.rds moved to [", file.path(dir_output, "landscapes.rds"), "]"))

  # remove temp raster in case remove_temp_rasters is TRUE
  if (remove_temp_rasters){
    unlink(dir_temp_raster, recursive = TRUE)
    print(paste(dir_temp_raster, "removed sucessfully"))
  }

  metadata_file_loc <- file.path(dir_input,"METADATA.txt")
  if (file.exists(metadata_file_loc)){
    print(paste0("[OK] METADATA.txt found! [", metadata_file_loc, "]"))
    # update metadata.txt NOTE: This might demand manual changes
    file.copy(metadata_file_loc,
              file.path(dir_output, basename(dir_input), "METADATA.txt"),
              overwrite=TRUE)
    print("Metadata.tx moved. Please update METADATA.txt properly!")

  } else {
    warning(paste0("[MISSING] METADATA.txt is missing! \n Create METADATA.txt manually \n at [", metadata_file_loc, "]"))
  }
  return(paste0("Landscape decompressed sucessfully to [", file.path(dir_output, basename(dir_input)),"]" ))
}









#' Title
#'
#' @param dir_input
#' @param dir_output
#' @param cost_function
#' @param remove_input
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#' # get gcf containing list of functions
source("./R/cost_functions.R")

compress_space <- function(dir_input=NULL,
                           dir_output=NULL){

  prepare_dirs(dir_input, dir_output)


  if (is.null(dir_output)){
    # error("Please provide an output directory for the recreated landscape")
    dir_output <- tempdir()
  }

  cat(paste0("Using output directory: [", dir_output, "]"))

  space_file_loc <- file.path(dir_input,"landscapes.rds")
  if (file.exists(space_file_loc)){
    print(paste0("[OK] landscapes.rds found: [", space_file_loc, "]"))
  } else{
    stop(paste0("[MISSING] landscapes.rds was not found: [", space_file_loc, "]"))
  }

  gen3sis_space=readRDS(space_file_loc)

  if (is.null(gen3sis_space)){
    error("Please provide a loaded gen3sis landscapes.rds as gen3sis_space")
  }
  saveRDS(gen3sis_space, file.path(dir_output, "landscapes.rds"), compress=T)

} #, timestoMa=6
