#' @title spac3tools: Tools for gen3sis input manipulation including conversion,
#' compression and decompression. This package is meat to provide support for
#' gen3sis previous and future generations, facilitating storage, modification
#' and reproducibility
#' General Engine for Eco-Evolutionary Simulations
#' @name spac3tools
#' @description
#' @references Development team
#' @details Gen3sis is implemented in a mix of R and C++ code, and wrapped into
#' an R-package. All high-level functions that the user may interact with are
#' written in R, and are documented via the standard R / Roxygen help files for
#' R-packages.
#' Runtime-critical functions are implemented in C++ and coupled to R via the
#' Rcpp framework.
#' Additionally, the package provides several convenience functions to generate
#' input data, configuration files and plots, as well as tutorials in the form
#' of vignettes that illustrate how to declare models and run simulations.
#' @seealso \code{\link{landscape_to_space}}   \code{\link{check_space}}  \code{\link{compress_space}}  \code{\link{decompress_space}}
#' @keywords programming IO iteration methods utilities
#' @concept spacial tools for inputs used by gen3sis modeling eco-evolutionary
#' macroevolution macroecology mechanisms or other flavours.
#' @examples
#' \dontrun{
#'
#' # 1. Load gen3sis and all necessary input data is set (landscape and config).
#'
#' library(gen3sis)
#'
#' # get path to example input inside package
#' datapath <- system.file(file.path("extdata", "WorldCenter"), package = "gen3sis")
#' path_config <- file.path(datapath, "config/config_worldcenter.R")
#' path_landscape <- file.path(datapath, "landscape")
#'
#' # 2. Run simulation
#'
#'sim <- run_simulation(config = path_config, landscape = path_landscape)
#'
#' # 3. Visualize the outputs
#'
#' # plot summary of entire simulation
#' plot_summary(sim)
#'
#' # plot richness at a given time-step
#' # this only works if species is saved for this time-step
#' landscape_t_150 <- readRDS(file.path(datapath,
#' "output", "config_worldcenter", "landscapes", "landscape_t_150.rds"))
#' species_t_150 <- readRDS(file.path(datapath,
#' "output", "config_worldcenter", "species", "species_t_150.rds"))
#' plot_richness(species_t_150, landscape_t_150)
#'
#' }
#' @docType package
#' @useDynLib gen3sis, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import Matrix
NULL





#' Create a space.rds with a gen3sis_space object
#' from a dir_input containing a landscape.rds and ideally a METADATA.txt
#' to a a space.rds into the dir_output.
#'
#' @param dir_input locatino of landscape.rds to be converted
#' @param dir_output location to store the converted space.rds
#' @param duration see \code{?create_space()}
#' @param area.unit see \code{?create_space()} only unit necessary, rest is calculated.
#' @param crs see \code{?create_space()}
#' @param cost_function
#' @param ... see \code{?create_space(...)}
#'
#' @return
#' @export
#'
#' @examples conv_landscape_to_space_help.R
landscape_to_space <- function(dir_input="C:/temp/4ds",
                               dir_output=dir_input,
                               duration=list(from=-65, to=0, by=1, unit="Ma"),
                               area.unit="km2",
                               crs="+proj=longlat +datum=WGS84 +no_defs",
                               cost_function=list(xx=gcf$XXHarderNA_dist_Km),...){

  prepare_dirs(dir_input, dir_output)
  landscape_file_loc <- file.path(dir_input,"landscapes.rds")
  if (file.exists(landscape_file_loc)){
    print(paste0("[OK] landscapes.rds found: [", landscape_file_loc, "]"))
  } else{
    stop(paste0("[MISSING] landscapes.rds was not found: [", landscape_file_loc, "]"))
  }
  lc <- readRDS(landscape_file_loc)
  all_timesteps <- rev(paste0(seq(duration$from, duration$to, duration$by), duration$unit))
  lc <- lapply(lc, function(x){
    colnames(x)[-c(1,2)] <- all_timesteps
    return(x)
  })
  ex_r <- rasterFromXYZ(lc[[1]][,1:3])
  crs(ex_r) <- crs
  total_area <- conv_unit(sum(area(ex_r)[]), "km2", area.unit)
  n_sites <- ex_r@ncols*ex_r@nrows
  gs <- create_space(type="raster",
                     duration=duration,
                     area=list(total_area=total_area,
                               n_sites=n_sites,
                               unit=area.unit),
                     env=lc,
                     cost_function = cost_function,
                     geo_dynamic=TRUE
                     )

  check_space(gs)
  saveRDS(gs$env, file.path(dir_output, "space3.rds"), compress=T)


}
