## METADATA ===============================================================
## Description: functions to recreate (aka decompress) gen3sis input only from
## gen3sis former inputs
##
## R version: 4.2.2 for Windows
## Date: 2022-12-20 13:00:20
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##


#' Decompress a compressed gen3sis2 space object
#'
#' This function takes a `spaces.rds` file (compressed gen3sis2 format), reconstructs the
#' simulation-ready full costs distances, and saves them to the specified output directory.
#' This functions primary use is undoing the compact storage of distances in the `spaces.rds` object.
#'
#' @param dir_input The input directory containing the gen3sis2 spaces.rds
#' @param cost_function_index The index of the cost function to be used and as provided
#' at spaces.rds at meta$cost_function. Default is the first cost function.
#' Default is the declared cost_functions list at cost_lists
#' Note that different const_functions can be used and it's computation
#' depends only on gen3sis::create_spaces_raster function
#' @param dir_output Output directory.
#' @param remove_temp_rasters Boolean. If TRUE, delete temporary environment tif folder and files.
#'
#' @return passing confirmation statement
#' @export
#'
#' @example inst/examples/decompress_space_help.R

decompress_space <- function(dir_input=NULL,
                             cost_function_index=1,
                             dir_output=NULL,
                             remove_temp_rasters=FALSE
){

  if (is.null(dir_output)){
    stop("Please provide an output directory for the recreated landscape")
  }

  cat(paste0("Using output directory: [", dir_output, "]"))

  space_file_loc <- file.path(dir_input,"spaces.rds")
  if (file.exists(space_file_loc)){
    cat(paste0('[OK] spaces.rds found: "', space_file_loc, '"'))
  } else{
    stop(paste0('[MISSING] spaces.rds was not found: "', space_file_loc, '"'))
  }

  gen3sis_space=readRDS(space_file_loc)
  ls <- gen3sis_space$env

  # create temp dir
  dir_temp_raster=file.path(dir_output, "temp_rasters")
  # create raster bricks
  dir.create(dir_output, showWarnings = FALSE)
  # create temp dir for temp rasters
  dir.create(dir_temp_raster, showWarnings = FALSE)

  # create temp rasters
  for (var_i in names(ls)){
    rstack <- list()
    for (c_i in 3:ncol(ls[[var_i]])){
      rstack <- c(rstack, terra::rast(ls[[var_i]][,c(1,2,c_i)], type="xyz"))
    }
    rstack <- terra::rast(rstack)

    terra::writeRaster(
      rstack,
      filename = file.path(
        dir_temp_raster,
        paste0(var_i, ".tif")
      ),
      filetype = "GTiff",
      overwrite = TRUE
    )
  }
  print(paste("Temporary Raster Bricks Saved to: ", dir_temp_raster))
  # load raster bricks
  # list all raster bricks ending with .tif
  bf <- paste0(names(ls), ".tif")
  print(paste("Raster Bricks are: ", paste(bf, collapse = "; ")))

  # load all temp raster bricks
  b <- vector(mode = "list", length = length(bf))
  for (i in seq_along(bf)){
    b[[i]] <- terra::rast(file.path(dir_temp_raster, bf[i]))
  }
  # prepare list
  lsn <- vector(mode = "list", length = length(ls))
  lsn <- setNames(lsn, names(ls))
  # attribute to list
  for (i in 1:length(lsn)){
    for (j in seq(terra::nlyr(b[[i]]))) {
      lsn[[i]] <- c(lsn[[i]], b[[i]][[j]])
    }
  }

  gsd <- gen3sis_space$meta$duration
  gen3sis2::create_spaces_raster(raster_list = lsn,
                                 cost_function = ifelse(cost_function_index == 0, gen3sis_space$meta$cost_function, gen3sis_space$meta$cost_function[[cost_function_index]]),
                                 directions=8, output_directory = file.path(dir_output,"decompressed"),
                                 duration = gsd, gsd$unit, full_dists = T, geodynamic=gen3sis_space$meta$geodynamic,
                                crs=gen3sis_space$meta$crs, verbose=T, overwrite_output = TRUE)

  saveRDS(gen3sis_space, file.path(dir_output, "spaces.rds"))
  print(paste0("spaces.rds moved to [", file.path(dir_output, "spaces.rds"), "]"))

  # remove temp raster in case remove_temp_rasters is TRUE
  if (remove_temp_rasters){
    unlink(dir_temp_raster, recursive = TRUE)
    print(paste(dir_temp_raster, "removed sucessfully"))
  }

  return(paste0("Space decompressed sucessfully to [", dir_output,"]" ))
}

#' Compress a gen3sis2 Space object by removing cost distances
#'
#' This function removes large cost distance matrices from a `spaces.rds` file,
#' creating a lighter version suitable for storage and distribution. All other
#' metadata and environmental data are preserved.
#' @name compress_space
#' @param dir_input Directory containing the `spaces.rds` to compress.
#' @param dir_output Directory where the compressed version will be saved.
#'
#' @return Saves a compressed `spaces.rds` file to `dir_output`. Does not return an R object.
#' @export
#'
#' @example inst/examples/compress_space_help.R

compress_space <- function(dir_input=NULL,
                           dir_output=NULL){

  gen3sis2:::prepare_dirs(dir_input, dir_output)


  if (is.null(dir_output)){
    stop("Please provide an output directory for the recreated space")
  }

  cat(paste0("Using output directory: [", dir_output, "]"))

  space_file_loc <- file.path(dir_input,"spaces.rds")
  if (file.exists(space_file_loc)){
    print(paste0("[OK] spaces.rds found: [", space_file_loc, "]"))
  } else{
    stop(paste0("[MISSING] spaces.rds was not found: [", space_file_loc, "]"))
  }

  gen3sis_space=readRDS(space_file_loc)

  if (is.null(gen3sis_space)){
    stop("Please provide a loaded gen3sis2 spaces.rds as gen3sis_space")
  }
  saveRDS(gen3sis_space, file.path(dir_output, "spaces.rds"), compress=T)

}
