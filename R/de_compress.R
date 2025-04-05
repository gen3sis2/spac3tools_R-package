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
#' @param cost_function the index of the cost function to be used and as provided
#' at spaces.rds at meta$cost_function. Default is the first cost function.
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
#' @examples
#' # Refer to inst/examples/run_recreate_space.R
#' # TODO example

decompress_space <- function(dir_input=NULL,
                             cost_function_index=1,
                             dir_output=NULL,
                             remove_temp_rasters=FALSE
                             #, timestoMa=6 # TODO
                             # SOLVE TEMPORAL ISSUE WITH FILE STRUCTUR
                             # PROB WITH TEMPORAL AND spatial resolution
){

  if (is.null(dir_output)){
    # error("Please provide an output directory for the recreated landscape")
    dir_output <- tempdir()
  }

  cat(paste0("Using output directory: [", dir_output, "]"))

  space_file_loc <- file.path(dir_input,"spaces.rds")
  if (file.exists(space_file_loc)){
    print(paste0("[OK] spaces.rds found: [", space_file_loc, "]"))
  } else{
    stop(paste0("[MISSING] spaces.rds was not found: [", space_file_loc, "]"))
  }

  gen3sis_space=readRDS(space_file_loc)
  ls <- gen3sis_space$env

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
      # c_i <- 4
      # rstack[[c_i-2]] <- raster::rasterFromXYZ(ls[[var_i]][,c(1,2,c_i)])
      rstack[[c_i-2]] <- terra::rast(ls[[var_i]][,c(1,2,c_i)], type="xyz")
    }
    rstack <- terra::rast(rstack)
                #raster::brick(rstack)
    terra::writeRaster(rstack, filename = file.path(dir_temp_raster, paste0(var_i, ".grd")), overwrite=TRUE)
  }
  print(paste("Temporary Raster Bricks Saved to: ", dir_temp_raster))
  # load raster bricks
  #list all temp raster bricks ending with .grd
  #bf <- list.files(dir_temp_raster, pattern=".grd$")
  bf <- paste0(names(ls), ".grd")
  print(paste("Raster Bricks are: ", paste(bf, collapse = "; ")))
  #load all temp raster bricks

  b <- NULL
  for (i in 1:length(bf)){
    #b[[i]] <- raster::brick(file.path(dir_temp_raster, bf[i]))
    b[[i]] <- terra::rast(file.path(dir_temp_raster, bf[i]))
  }
  #prepare list
  lsn <- lapply(ls, function(x){x <- NULL})
  #attribute to list
  for (i in 1:length(lsn)){
    for (il in 1:dim(b[[i]])[3]) {
      lsn[[i]] <- c(lsn[[i]], b[[i]][[il]])
    }
  }
    # myt <- (as.numeric(colnames(ls[[1]][-c(1,2)])))/timestoMa
  # string_time_step <- paste0(formatC(round(myt,2), width=5, flag="0", digits=2, format="f"),"Ma" )
  gsd <- gen3sis_space$meta$duration
  gen3sis2::create_spaces_raster(raster_list = lsn,
                                 cost_function = gen3sis_space$meta$cost_function[[cost_function_index]],
                                 directions=8, output_directory = file.path(dir_output,"decompressed"),
                                 duration = gsd, gsd$unit, full_dists = T, geodynamic=gen3sis_space$meta$geodynamic,
                                crs=gen3sis_space$meta$crs, verbose=T, overwrite_output = TRUE)
  # TODO remove this line after create_input_landscape is fixed to create_space_raster
  unlink(file.path(dir_output, "landscapes.rds"))
  saveRDS(gen3sis_space, file.path(dir_output, "spaces.rds"))
  print(paste0("spaces.rds moved to [", file.path(dir_output, "spaces.rds"), "]"))

  # remove temp raster in case remove_temp_rasters is TRUE
  if (remove_temp_rasters){
    unlink(dir_temp_raster, recursive = TRUE)
    print(paste(dir_temp_raster, "removed sucessfully"))
  }

  # metadata_file_loc <- file.path(dir_input,"METADATA.txt")
  # if (file.exists(metadata_file_loc)){
  #   print(paste0("[OK] METADATA.txt found! [", metadata_file_loc, "]"))
  #   # update metadata.txt NOTE: This might demand manual changes
  #   file.copy(metadata_file_loc,
  #             file.path(dir_output, basename(dir_input), "METADATA.txt"),
  #             overwrite=TRUE)
  #   print("Metadata.tx moved. Please update METADATA.txt properly!")
  #
  # } else {
  #   warning(paste0("[MISSING] METADATA.txt is missing! \n Create METADATA.txt manually \n at [", metadata_file_loc, "]"))
  # }
  return(paste0("Space decompressed sucessfully to [", dir_output,"]" ))
}









#' compress a gen3sis_space, by removing the cost distances
#'
#' @param dir_input
#' @param dir_output
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

  gen3sis2:::prepare_dirs(dir_input, dir_output)


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
