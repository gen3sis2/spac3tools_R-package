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
#' @seealso \code{\link{landscape_to_space}}   \code{\link{gen3sis2::check_spaces}}  \code{\link{compress_space}}  \code{\link{decompress_space}}
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





#' Create a space.rds object
#' from a dir_input containing a landscape.rds
#' to a a space.rds into the dir_output.
#'
#' @param dir_input Location of landscape.rds to be converted
#' @param dir_output location to store the converted space.rds
#' @param duration see \code{?gen3sis2::create_spaces}
#' @param crs see \code{?gen3sis2::create_spaces}
#' @param cost_function see \code{?create_spaces} Default it a cost of 2 for sites with NA.
#' @param ... see \code{?gen3sis2::create_spaces}
#'
#' @return
#' @export
#'
#' @examples
#' # Refer to inst/examples/conv_landscape_to_space_help.R
#' # TODO example
landscape_to_space <- function(dir_input=NA,
                               dir_output=dir_input,
                               duration=list(from=NA, to=NA, by=NA, unit="Ma"),
                               crs="+proj=longlat +datum=WGS84 +no_defs",
                               cost_function=list(xx=spac3tools::gcf$XXHarderNA_dist_Km),...){

  gen3sis2:::prepare_dirs(dir_input, dir_output)
  landscape_file_loc <- file.path(dir_input,"landscapes.rds")
  if (file.exists(landscape_file_loc)){
    print(paste0("[OK] landscapes.rds found: [", landscape_file_loc, "]"))
  } else{
    stop(paste0("[MISSING] landscapes.rds was not found: [", landscape_file_loc, "]"))
  }
  lc <- readRDS(landscape_file_loc)
  # create order flipped to match old landscape object order
  all_timesteps <- rev(paste0(seq(duration$from, duration$to, duration$by), duration$unit))
  if (length(colnames(lc[[1]])[-c(1,2)])!=length(all_timesteps)){
    stop("Mismatch on number of time steps, check your duration and env parameters")
  }
  lc <- lapply(lc, function(x){
    colnames(x)[-c(1,2)] <- all_timesteps
    # revert the order of the columns, so that space objects time, moves from left to right
    x <- x[,c(colnames(x)[1:2],rev(colnames(x)[-c(1,2)]))]
    return(x)
  })
  # get dummy raster for area calculations
  ex_r <- terra::rast(lc[[1]][,1:3], type="xyz")
  terra::crs(ex_r) <- crs
  total_area <- sum(terra::cellSize(ex_r)[]) # area is fix to Km2
  # get number of rows
  n_sites <- dim(ex_r)[1]
  gs <- gen3sis2::create_spaces(env=lc,
                     type="raster",
                     duration=duration,
                     area=list(extent=terra::ext(ex_r)[],
                               total_area=total_area,
                               n_sites=n_sites,
                               unit="km2"),
                     cost_function = cost_function,
                     geodynamic=NULL,
                     type_spec=list("res"=terra::res(ex_r)),
                     ...
                     )
  gen3sis2::check_spaces(gs)
  saveRDS(gs, file.path(dir_output, "spaces.rds"), compress=T)
  print(paste0("space.rds type=", gs$type , " saved to [", file.path(dir_output, "spaces.rds"), "]"))
}


#' Convert a space.rds type raster to type h3
#' It is recommended to use this on top of the smallest available raster resolution
#'@param dir_input Location of space.rds type raster
#'@param dir_output location to store the converted space.rds of type h3
#'@param res resolution of the h3 grid. This parameter is passed directly to the
#' `res` parameter of `hrjsr::get_children()`. For details refer to h3 documentation
#'@param verbose integer, 0 for no messages, 1 for some messages, 2 for all messages
#'@compute_distances logical, if TRUE, the cost distances are computed and stored (default TRUE)
#'@return
#'@export
#'@examples conv_space_raster_to_h3_help.R
space_raster_to_h3 <- function(dir_input="C:/temp/decompressed_spaces/world60by10at4d",
                               dir_output=file.path(dir_input, "h3_v2"),
                               res=0,
                               verbose=0,
                               compute_distances=TRUE){
  # res is the resolution of the h3 grid in m
  # level is the h3 level
  # space_raster is a raster object
  # returns a data.frame with h3 index and the corresponding value
  # the value is the mean
  gen3sis2:::prepare_dirs(dir_input, dir_output)

  # origin space
  o_gls_all <- readRDS(file.path(dir_input, "spaces.rds"))
  # origin full distances
  o_cd_fl  <- list.files(file.path(dir_input, "distances_full"))
  if (length(o_cd_fl)==0){
    stop("No full distances found")
  }
  # order o_cd_fl
  numb_ts <- as.numeric(gsub("[^0-9.-]+", "", o_cd_fl))
  o_cd_fl <- o_cd_fl[order(numb_ts)]

  # create dir_output distances_full folder
  dir.create(file.path(dir_output, "distances_full"), showWarnings = FALSE)
  time_steps <- colnames(o_gls_all$env[[1]][-c(1,2)])
  n_ts <- length(time_steps)

  if (o_gls_all$meta$geodynamic & n_ts!=length(o_cd_fl)){
    #stop function and show error message
    stop("Mismatch on number of time steps and full distances.rds")
  }
  all_pts <- o_gls_all$env[[1]][,c("x", "y")]
  # set as coordinates
  all_pts_sf <- sf::st_as_sf(all_pts, coords = c("x","y"), crs = o_gls_all$meta$crs)
  all_pts_4326 <- sf::st_transform(all_pts_sf, 4326) # transform to 4326

  # h3_fill: cells corresponding to a convex polygon around all points
  h3_fill <- all_pts_4326 |>
    sf::st_union() |>
    sf::st_convex_hull() |>
    h3jsr::polygon_to_cells(res=res)

  h3_fill <- unlist(h3_fill) |>
    h3jsr::get_disk(2) |>
    unlist() |>
    unique()

  # h3_cell_idx: cells corresponding only to the points itself
  h3_cell_idx <- h3jsr::point_to_cell(all_pts_4326, res = res)

  # if (sum(duplicated(h3_cell_idx)) > 0){
  #   warning("Consider increasing the resolution")
  # }

  # prepare variables
  p1 <- h3jsr::cell_to_point(h3_address = unlist(h3_fill), simple = FALSE)

  # make H3 simple
  h3rezz_pts <- as.data.frame(sf::st_coordinates(p1$geometry))
  table_sites <- unlist(lapply(list("all_pts"=all_pts, "h3rezz_pts"=h3rezz_pts), nrow))
  if (verbose>0){
    cat("# Number of sites in space type:")
    cat(names(table_sites))
    cat(table_sites)
  }
  #####  run loop over each h3rezz point to find closest raster point
  # number of points
  np <- table_sites["h3rezz_pts"] # same as length(p1$h3_address)
  h3index <- p1$h3_address
  # get all raster indexes
  all_r_index <- rownames(o_gls_all$env[[1]])
  # closest raster index, empty vector
  c_r_index <- rep(NA, np)
  #create numerical empty vector
  error_m <- c_r_index # store error distance in meters
  # create a dictionary for h3 index
  cell_dict <- list()
  ## short initial loop to finds closest points... only do once
  for (h3pi in 1:np){ # loop for each point
    # h3pi <- 1
    if (verbose>0){
      print(paste0("h3pi or h3 index of points is: ", h3pi, "/", np))
    }
    #get closest point
    distances <- geosphere::distGeo(h3rezz_pts[h3pi,], all_pts)
    m_dist <- min(distances)
    closest_p_index <- all_r_index[distances==m_dist]
    if (length(closest_p_index)!=1){
      closest_p_index <- closest_p_index[1]
      # giver warning message
      warning("Error: more than one point found")
    }
    c_r_index[h3pi] <- closest_p_index
    error_m[h3pi] <- m_dist

    # store h3 index in dictionary
    h3_code <- unlist(h3_fill)[h3pi]
    cell_dict[[h3_code]] <- closest_p_index
  }

  if(c_r_index |> duplicated() |> sum() != 0){
    cat("Resolving duplicated cells...")

    inverted_list <- list()
    for (name in names(cell_dict)) {
      h3_indexes <- cell_dict[[name]]
      if (is.null(inverted_list[[h3_indexes]])) {
        inverted_list[[h3_indexes]] <- c(name)
      } else {
        inverted_list[[h3_indexes]] <- c(inverted_list[[h3_indexes]], name)
      }
    }

    for (forgotten_point in which(!(1:nrow(all_pts) %in% as.numeric(names(inverted_list))))){
      inverted_list[[as.character(forgotten_point)]] <- NA
    }

    dist_resolved <- rep(NA, length(c_r_index))
    new_error <- c()
    #for (raster_idx in c_r_index) {
    for(raster_idx in names(inverted_list)){
      print(raster_idx)
      raster_point <- all_pts[as.numeric(raster_idx),]
      distances <- geosphere::distGeo(raster_point, h3rezz_pts)
      min_dist <- min(distances)
      new_error <- c(new_error, min_dist)
      closest_cell_idx <- which(distances == min_dist)

      dist_resolved[closest_cell_idx] <- raster_idx
    }

    if((dist_resolved[!is.na(dist_resolved)] |> duplicated() |> sum() > 0)
       || (length(dist_resolved[!is.na(dist_resolved)]) != nrow(all_pts))){
      stop("Could not resolve duplicates. Try increasing the resolution.")
    } else {
      c_r_index <- dist_resolved
      error_m <- new_error
    }
  }

  #closest_raster_index <<- c_r_index # debug
  if (verbose>1){
    hist(error_m/1000, main="Error in Km")
  }
  # prepare input variables...
  dummyM <- matrix(rep(NA, np*n_ts), ncol=n_ts)
  rownames(dummyM) <- c_r_index # latest change
  if(any(is.na(rownames(dummyM)))){
    point_starter <- sum(!is.na(rownames(dummyM)))+1
    for(i in which(is.na(row.names(dummyM)))){
      rownames(dummyM)[i] <- as.character(point_starter)
      point_starter <- point_starter + 1
    }
  }
  colnames(dummyM) <- time_steps
  envs <- vector(mode = "list", length = length(o_gls_all$env))
  envs_names <- c(names(o_gls_all$env))
  names(envs) <- envs_names
  envs <- lapply(envs, function(x){x <- dummyM})

  # loop over time steps and update envs
  for (ti in 1:n_ts){
    # ti <- 1
    num_ts <- (n_ts-1):0
    if(verbose>0){
      print(paste(ti,"of ", n_ts, "time_steps:", time_steps[ti]))
    }
    # store gen3sis value
    # get envs of closest index
    for (env_i in envs_names){
      # env_i <- envs_names[1]
      envs[[env_i]][,time_steps[ti]] <- o_gls_all$env[[env_i]][c_r_index,time_steps[ti]]
    }
    # set the presences for the time-steps
    lp1_m <- !is.na(envs[[1]][,time_steps[ti]])
    # t_p <- time_steps[ti]
    # lp1_m <- apply(envs[[1]][,c(t_p, t_f)], 1, function(x) any(!is.na(x)))
    # old lp1_m <- !is.na(envs[[1]][,time_steps[ti]]) # per time step


    # load distances and manipulate it...
    # if geodynamic is true or geo is false and ti is 1
    if (compute_distances & (o_gls_all$meta$geodynamic | (!o_gls_all$meta$geodynamic&ti==1))){
      if (o_gls_all$meta$geodynamic){
        tiis <- rev(o_cd_fl)[ti]
      }else{
        tiis <- "distances_full_0.rds"
      }
      # lfd_ti = landscapes full distances at ti
      lfd_ti <- readRDS(file.path(dir_input, "distances_full", tiis)) # TODO FIX THE
      print(paste("Loaded:", rev(o_cd_fl)[ti]))
      # sub select cost distances
      #ib <- colnames(lfd_ti)%in%c_r_index[lp1_m] # TODO USE h3c1 index later here too, do match#
      # ix <- which(c_r_index%in%colnames(cd))
      # length(ib)
      new_cost_dist_full <- lfd_ti[c_r_index[lp1_m],c_r_index[lp1_m]]
      #dimnames(new_cost_dist_full) <- list(h3index[lp1_m], h3index[lp1_m])

      # h3_ii <- which(c_r_index%in%colnames(new_cost_dist_full))
      # dimnames(new_cost_dist_full) <- list(h3index[h3_ii], h3index[h3_ii])

      # save cost function ti
      saveRDS(new_cost_dist_full, file.path(dir_output, "distances_full", tiis))
      if (verbose>0){
        print(paste("Saved:", file.path(dir_output, "distances_full", tiis)))
      }
    }
    plot_stuff <- FALSE
    if (plot_stuff){
      #possible points
      terra::plot(terra::rast(o_gls_all$env[[1]][,1:3], type="xyz"))
      # points(all_pts, pch=2, col=rgb(1,1,1,1,1))
      points(h3rezz_pts[lp1_m, ], pch=3)
      # plot points in original raster

      # plot all possible points
      plot_points(xy=h3rezz_pts,
                  vcol=1,
                  cols=terrain.colors(2),
                  vcex=1)

      # plot only points with temperature on h3
      plot_points(xy=h3rezz_pts[lp1_m, ],
                  vcol=envs[[1]][lp1_m],
                  cols=gen3sis::color_richness(5),
                  vcex=1.5)
      # plot points with temperature on original raster input
      plot_points(xy=all_pts_sf$geometry[!is.na(o_gls_all[[1]][,3]) ],
                  vcol=na.omit(o_gls_all[[1]][,3]),
                  cols=gen3sis::color_richness(30),
                  vcex=1)

      # plot cost distance values with the original location...
      # if (compute_distances) {
      # points(all_pts[colnames(new_cost_dist_full),], pch=3)
      # # plot cost distance location....
      # points(all_pts[colnames(lfd_ti),], pch=3)
      #
      # # plot random distances
      #
      # # plot only points with temperature and plot all points with distances on top...
      # df <- all_pts[colnames(new_cost_dist_full),]
      # plot_points(xy=h3rezz_pts,
      #             vcol=envs$gtemp, # use 1 for all points
      #             cols=gen3sis::color_richness(5),
      #             vcex=1.5)
      #
      # #  plot(df, col=rainbow(800)) # :) have a nice day
      #
      # # select random points...
      # points <- sample(colnames(new_cost_dist_full), 2)
      # points(df, pch=3)
      # segments(x0=df[points[1], "x"],y0=df[points[1], "y"],x1=df[points[2], "x"],y1=df[points[2], "y"])
      # round(new_cost_dist_full[points[1],points[2]], 2)
      # text(50, 50, )
      #
      #
      # df <- all_pts[colnames(new_cost_dist_full),]
      # plot_points(xy=df,
      #             vcol=1, # use 1 for all points
      #             cols=gen3sis::color_richness(5),
      #             vcex=1.5)
      # }
    }
  } # END LOOP OVER TIME STEPS
  # save final landscape
  h3pts <- h3rezz_pts
  #rownames(h3pts) <- h3c1
  colnames(h3pts) <- c("x", "y")
  final_envs <- lapply(envs, function(x){
    #row.names(x) <- NULL # set rownames to h3c1 # debug
    return(cbind(h3pts, x))
  })
  # copy old landscape as ref.
  final_space <- o_gls_all
  final_space$env <- final_envs
  final_space$meta$type <- "h3"
  final_space$meta$type_spec <- list("res"=res)
  final_space$meta$area$total_area <- sum(h3jsr::cell_area(h3_fill, final_space$area$unit, simple=TRUE))
  final_space$meta$area$n_sites <- as.numeric(np)
  saveRDS(final_space, file.path(dir_output, "spaces.rds"))
  return(paste0("Space converted to h3 and saved to [", dir_output,"]" ))
}

#' Do local conversion of a space.rds type h3 to type points Locally
#' @param dir_input Location of space.rds type h3 to be converted to points
space_h3_to_points <- function(dir_input){
  space <- readRDS(file.path(dir_input, "spaces.rds"))
  space$meta$type <- "points"
  saveRDS(space, file.path(dir_input, "spaces.rds"))
}

