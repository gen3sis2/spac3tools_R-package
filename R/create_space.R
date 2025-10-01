#' create an spaces input from a named list of H3 datasets and user defined cost function
#'
#' @param h3_list named list  of H3 datasets. Starting from the past towards the present.
#' NOTE: the list names are important since these are the environmental names
#' @param cost_function function that returns a cost value between a pair of sites (neighbors) that should have the following signature:
#'  \code{cost_function <-
#'  function(src, src_habitable, dest, dest_habitable){
#'  rules for environmental factors to be considered (e.g. elevation)
#'  return(cost value)
#' }}
#' where: **src** is a vector of environmental conditions for the origin sites,
#' **src_habitable** (TRUE or FALSE) for habitable condition of the origin sites,
#' **dest** is a vector of environmental conditions for the destination site, dest_habitable  (TRUE or FALSE) for habitable condition of the destination cell
#' @param output_directory path for storing the gen3sis ready space (i.e. space.rds, metadata.txt and full- and/or local_distance folders)
#' @param timesteps vector of names for every time-step to represent the time-step at gen3sis ready space.
#' If timesteps=NULL (default), time-steps are sequentially numbered from 0 to the latest time-step.
#' @param full_dists should a full distance matrix be calculated? TRUE or FALSE? Default is FALSE.
#' If TRUE calculates the entire distance matrix for every time-step and between all habitable cells
#' (faster CPU time, higher storage required).
#' If FALSE (default), only local distances are calculated (slower CPU time when simulating but smaller gen3sis space size)
#' @param crs the coordinate reference system in crs format (see raster::crs). Default is defined by \code{\link{gen3sis2::create_spaces}}
#' @param overwrite_output TRUE or FALSE
#' @param verbose print distance calculation progress (default: FALSE)
#' @param duration list with from, to, by and unit. Default is from -latest time to zero by 1 Ma
#' @param geodynamic True or False, if the space is dynamic (e.g. sea-level change) or static. Default is NULL,
#' i.e. deciding final value based on the input data using \code{?is_geodynamic}.
#' @param author author of the space, see \code{?create_spaces}
#' @param source source of the space, see \code{?create_spaces}
#' @param description list with env and methods, see \code{?create_spaces}
#' @returns no return object. This function saves the space input files for gen3sis at the output_directory
#'
#' @importFrom gen3sis2 create_spaces check_spaces
#' @importFrom h3jsr cell_area get_res
#'
#' @export
#' @examples # TODO
create_spaces_h3 <- function(
    h3_list,
    cost_function,
    output_directory,
    # timesteps = NULL,
    full_dists = FALSE,
    overwrite_output = FALSE,
    verbose = FALSE,
    duration=list(from=NA, to=NA, by=NA, unit="Ma"),
    geodynamic=NULL,
    ...
    ) {
  # habitability masks removed for now.. assumin NA's from input data as non-habitable (applies if only one env. variable has NA in the cell)!
  habitability_masks = NULL

  # prepare directories
  gen3sis2:::create_directories(output_directory, overwrite_output, full_dists)

  # compute time-steps
  if (any(is.na(duration)) || all(!names(duration) %in% c("from", "to", "by", "unit"))) {
    warning("Duration is ideally informed as a list with from, to, by and unit.")
    timesteps <- -(length(setdiff(names(h3_list[[1]]), c("x","y","h3_address"))) - 1):0

    if(is.list(duration) && "unit" %in% names(duration)){
      warning("Assuming default duration from -latest time to zero by 1 in given unit.")
      duration <- list(from=timesteps[1], to=0, by=1, unit=duration$unit)
    } else {
      warning("Assuming default duration from -latest time to zero by 1 Ma.")
      duration <- list(from=timesteps[1], to=0, by=1, unit="Ma")
    }
  }

  timesteps <- paste0(seq(duration$from, duration$to, by = duration$by), duration$unit)

  # check h3 addresses consistency
  h3_cells <- sapply(names(h3_list), function(v){h3_list[[v]][["h3_address"]]})

  if(all(apply(h3_cells,1,function(cell){length(unique(cell)) == 1}))) {
    h3_cells <- h3_cells[,1]
  } else {
    stop("Inconsistencies found in H3 addresses between variables.")
  }

  # prepare and save spaces
  compiled_env <- lapply(h3_list, function(v){
    v <- v[,-which(colnames(v) == "h3_address")]
    v <- v[,c("x","y",setdiff(colnames(v),c("x","y")))]
    colnames(v) <- c("x","y",timesteps)
    v
  })

  ts_habitabilty <- list()
  for (ts in timesteps) {
    hab_location <- sapply(compiled_env, function(v){
      !is.na(v[,ts])
    })

    hab_mask <- apply(hab_location, 1, all)
    l <- list(hab_mask)
    names(l) <- ts
    ts_habitabilty <- append(ts_habitabilty, l)

    for (v in 1:length(compiled_env)) {
      compiled_env[[v]][!hab_mask,ts] <- NA
    }
  }
  names(compiled_env) <- names(h3_list)

  gs <- gen3sis2::create_spaces(env=compiled_env,
                                type="h3",
                                duration=duration,
                                area=list(extent=NA,
                                          total_area=NA,
                                          n_sites=NA,
                                          unit="km2"),
                                geodynamic=geodynamic,
                                cost_function = list(cost_function),
                                ...
  )

  # filling spaces
  total_area <- h3jsr::cell_area(h3_cells, unit = "km2", simple = T) |> sum(na.rm = T)
  n_sites <- length(h3_cells)
  gs$meta$area$total_area <- total_area
  gs$meta$area$n_sites <- n_sites
  gs$meta$area$extent <- c("xmin" = min(compiled_env[[1]][["x"]]),
                           "xmax" = max(compiled_env[[1]][["x"]]),
                           "ymin" = min(compiled_env[[1]][["y"]]),
                           "ymax" = max(compiled_env[[1]][["y"]]))
  gs$meta$type_spec <- list(res = h3jsr::get_res(h3_cells[[1]]))

  gen3sis2::check_spaces(gs)

  # in case geodynamic is set to FALSE, double check env matrix
  if (!geodynamic){
    # if compiled_env is dynamic, reset it
    if (gen3sis2:::is_geodynamic(compiled_env)){ # get the geodynamic status
      warning("geodynamic is set to FALSE but environment says otherwise.
          changing geodynamic to TRUE")
      geodynamic <- TRUE
      gs$meta$geodynamic <- TRUE
    }
  }

  # save spaces.rds
  saveRDS(gs, file.path(output_directory, "spaces.rds"))

  # create local distances
  # iterate over times-teps

  # number of time-steps
  if (geodynamic){
    nts <- length(timesteps)
  } else {
    nts <- 1 # to only compute the first
  }

  coords <- h3_list[[1]][,c("x","y")]
  for( step in 1:nts ) {
    if (verbose) {
      cat(paste("starting distance calculations for timestep", step, '\n'))
    }

    var_step <- do.call(cbind, lapply(compiled_env, function(v){v[[timesteps[step]]]}))
    var_step <- cbind(coords,var_step)

    habitable_mask <- ts_habitabilty[[timesteps[step]]]

    distance_local <- get_h3_distances(var_step, h3_cells, habitable_mask, cost_function)

    file_name <- paste0("distances_local_", as.character(nts-step), ".rds")
    saveRDS(distance_local, file = file.path(output_directory, "distances_local", file_name))

    if(full_dists){
      # transpose to preserve src/dest relation for efficent local traversal in get_distance_matrices function
      distance_local <- t(distance_local)

      habitable_cells <- which(habitable_mask[])

      dist_matrix <- gen3sis2:::get_distance_matrix(
        habitable_cells,
        length(habitable_mask[]),
        distance_local@p,
        distance_local@i,
        distance_local@x,
        Inf
      )

      file_name <- paste0("distances_full_", as.character(nts-step), ".rds")
      saveRDS(dist_matrix, file = file.path(output_directory, "distances_full", file_name))
      rm(dist_matrix)
      gc()
    }
  }
}

#' Creates a list of H3 data.frame
#'
#' @param h3_address character. A set of H3 cell addresses. If provided
#'   together with \code{res}, the resolution is ignored.
#' @param res integer. H3 resolution at which to generate cells covering the
#'   globe. Only used if \code{h3_address} is \code{NULL}.
#' @param raster_list list of named list(s) of raster(s) or raster file(s) name(s). Starting from the past towards the present.
#' NOTE: the list names are important since these are the environmental names
#'
#' @returns A named list of data frames. Each data frame corresponds to one
#'   element of \code{raster_list} and contains:
#'   \itemize{
#'     \item the H3 cell identifier,
#'     \item extracted raster values for each layer,
#'     \item centroid coordinates (\code{x}, \code{y}).
#'   }
#'
#' @importFrom h3jsr get_children get_res0 cell_to_point
#' @importFrom sf st_coordinates st_set_geometry
#' @importFrom terra ext nlyr extract vect
#'
#' @export
#' @examples # TODO
data_raster_to_h3 <- function(
    h3_address = NULL,
    res = NULL,
    raster_list){
  if(!is.null(h3_address) && !is.null(res)){
    warning("Ignoring parameter res")
  } else if (is.null(h3_address) && !is.null(res)) {
    h3_address <- h3jsr::get_res0() |>
      h3jsr::get_children(res = res) |>
      unlist()
  } else if (is.null(h3_address) && is.null(res)) {
    stop("Either a set of cells or a resolution must be provided")
  }

  r_ext <- terra::ext(raster_list[[1]][[1]])

  points_sf <- h3jsr::cell_to_point(h3_address, simple = F)
  pts_before <- nrow(points_sf)

  coords <- sf::st_coordinates(points_sf)

  points_sf <- points_sf[
    coords[,1] >= r_ext[1] & coords[,1] <= r_ext[2] &  # longitude
      coords[,2] >= r_ext[3] & coords[,2] <= r_ext[4],   # latitude
  ]

  pts_after <- nrow(points_sf)

  if(pts_after < pts_before) {
    warning(pts_before-pts_after," cells outside the raster extent deleted.")
  }

  df_list <- lapply(names(raster_list), function(v){
    names(raster_list[[v]]) <- paste0("ts_",1:terra::nlyr(raster_list[[v]]))
    vals <- terra::extract(raster_list[[v]], terra::vect(points_sf))

    val_points <- cbind(points_sf, vals[,-1])
    val_coords <- sf::st_coordinates(val_points)
    colnames(val_coords) <- c("x","y")
    val_df <- val_points |>
      sf::st_set_geometry(NULL) |>
      cbind(val_coords)
    rownames(val_df) <- 1:nrow(val_df)

    val_df
  })
  names(df_list) <- names(raster_list)

  return(df_list)
}

#' Compute H3-based transition cost matrix
#'
#' This internal function calculates a sparse transition matrix between H3 cells
#' based on spatial coordinates, cell values, a habitable mask, and a user-defined
#' cost function. The function generates neighboring relationships, computes
#' geodesic distances, and applies corrections to account for local distances.
#'
#' @param var_step Data frame. Must contain columns \code{x} and \code{y} for
#'   coordinates, as well as other variables representing the state of each cell.
#' @param h3_cells Character vector. H3 cell IDs corresponding to the points in
#'   \code{var_step}.
#' @param habitable_mask Logical or numeric vector. Indicates which cells are
#'   considered habitable. Must have length equal to \code{nrow(var_step)}.
#' @param cost_function Function. A user-supplied function that takes two lists
#'   (\code{source_cell} and \code{destination_cell}) and returns a numeric cost.
#'   Each list contains \code{index}, \code{coordinates}, \code{value}, and
#'   \code{habitable}.
#'
#' @return A sparse matrix (\code{dgCMatrix}) representing the transition costs
#'   between cells. Rows correspond to destination cells and columns correspond
#'   to source cells. The matrix accounts for geodesic distances and user-defined
#'   costs, with zero values for non-traversable or infinitely costly transitions.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Generate neighbors for each H3 cell using a distance of 1 hexagon.
#'   \item Create a directed graph representing adjacency between cells.
#'   \item Compute geodesic distances between neighboring cells using \code{geosphere::distGeo}.
#'   \item Apply an inverse-distance correction to account for spatial separation.
#'   \item Iterate over all edges, computing transition costs using \code{cost_function}.
#'   \item Construct a sparse transition matrix combining costs and distance corrections.
#' }
#' This function is intended for internal use only and is not exported.
#'
#' @importFrom geosphere distGeo
#' @importFrom h3jsr get_disk
#' @importFrom igraph graph_from_edgelist V ecount as_edgelist E as_adjacency_matrix
#' @importFrom Matrix sparseMatrix drop0
#'
#' @noRd
#' @examples # TODO
get_h3_distances <- function(var_step, h3_cells, habitable_mask, cost_function){
  coords <- as.matrix(var_step[, c("x", "y")])

  neighs <- h3jsr::get_disk(h3_cells, 1)
  adj <- lapply(1:length(neighs), function(cell){
    neighborhood <- neighs[[cell]][neighs[[cell]] != h3_cells[cell]]
    neighborhood <- which(h3_cells %in% neighborhood)

    cell_adj <- expand.grid(from = cell, to = neighborhood)
  })
  adj <- do.call(rbind, adj) |> as.matrix()

  h3_graph <- igraph::graph_from_edgelist(as.matrix(adj), directed = TRUE)

  ##
  igraph::V(h3_graph)$name <- 1:nrow(var_step)
  edges_values <- numeric(igraph::ecount(h3_graph))

  edge_list_matrix <- igraph::as_edgelist(h3_graph, names = TRUE)
  edges_values <- edge_list_matrix[,2]
  igraph::E(h3_graph)$weight <- edges_values

  transition_matrix <- igraph::as_adjacency_matrix(
    graph = h3_graph,
    attr = "weight",
    sparse = TRUE
  ) |> t()

  #
  from_mtx <- var_step[adj[,1],c("x","y")] |> as.matrix()
  rownames(from_mtx) <- NULL

  to_mtx <- var_step[adj[,2],c("x","y")] |> as.matrix()
  rownames(to_mtx) <- NULL

  correction <- cbind(from_mtx, to_mtx)

  coords_A <- as.matrix(correction[, 1:2])
  coords_B <- as.matrix(correction[, 3:4])

  distances <- geosphere::distGeo(coords_A, coords_B)

  scaleValue <- 1
  correctionValues <- 1 / (distances / scaleValue)
  #

  i <- as.integer(adj[,2] - 1)
  j <- as.integer(adj[,1] - 1)
  xv <- as.vector(correctionValues) #check for Inf values!
  dims <- length(h3_cells)
  correctionMatrix <- new("dgTMatrix", i = i, j = j, x = xv, Dim = as.integer(c(dims,dims)))
  correctionMatrix <- (as(correctionMatrix,"sparseMatrix"))

  correctionMatrix@x <- 1 / correctionMatrix@x

  # return_list <- list(
  #   transitionMatrix = transition_matrix,
  #   correctionMatrix = correctionMatrix
  # )
  #
  # return(return_list)
  transition_cells <- summary(transition_matrix)

  tmp_cost <- numeric(nrow(transition_cells))

  habitable_mask <- as.logical(as.vector(habitable_mask))
  var_names <- names(var_step)[which(!names(var_step)%in%c("x","y"))]
  space_stack <- as.matrix(var_step[,var_names])
  for(k in 1:nrow(transition_cells)){
    ind_i <- transition_cells[k, "i"] # destination
    ind_j <- transition_cells[k, "j"] # origin

    coords_i <- coords[ind_i,] # destination coordinates
    coords_j <- coords[ind_j,] # origin coordinates

    cell_i <- space_stack[ind_i,] # destination values
    cell_j <- space_stack[ind_j,] # origin values

    habitable_i <- habitable_mask[ind_i] # is the destination habitable?
    habitable_j <- habitable_mask[ind_j] # is the origin habitable?

    # IMPORTANT
    # We want efficient iterations over sparse matrices
    # sparse matrices are column compressed, We therefore flip the indices and local distances are
    # indexed as local_distance[dest, src].
    # cost <- cost_function(cell_j, as.logical(habitable_j), cell_i, as.logical(habitable_i))
    source_cell <- list(
      index = ind_j,
      coordinates = coords_j,
      value = cell_j,
      habitable = habitable_j
    )

    destination_cell <- list(
      index = ind_i,
      coordinates = coords_i,
      value = cell_i,
      habitable = habitable_i
    )

    cost <- cost_function(source_cell, destination_cell)

    if(cost == Inf){
      cost <- 0
    }
    tmp_cost[k] <- cost
  }

  transition_matrix <- Matrix::sparseMatrix(
    i = transition_cells[, "i"],
    j = transition_cells[, "j"],
    x = tmp_cost
  )

  transition_matrix <- Matrix::drop0(transition_matrix) * correctionMatrix

  rownames(transition_matrix) <- 1:dim(space_stack)[1]
  colnames(transition_matrix) <- 1:dim(space_stack)[1]
  return(transition_matrix)
}

#' #' create empty gen3sis_spaces
#' #'
#' #' @param env named list of environmental variables in a with X,Y, and time-steps from
#' #' the most recent of future time-step
#' #' @param type string with type of gen3sis space. Accepted values are \code{check_space()$type}
#' #' @param duration list containing information on temporal dimension list(from, to, by, unit)
#' #' *from* is the oldest time-step; negative number if starting in the past, zero if starting in the present
#' #' *to* CAN ONLY BE smaller than *from*, since *to* is the latest time
#' #' *by* is the time interval increment. Note that this is constant and can only be positive;
#' #' *unit* is the time unit used. Accepted units are \code{check_space()$duration}
#' #' e.g. list(-20, 0, 1, "Ma") has a landscape that covers the last 20 Ma until the present, every 1 Ma.
#' #'      list(-800, 300, 10, "mil") has a landscape that covers the last 800 kyra or mil for millions of years
#' #'      and goes until the future 300 kya at every 100 mil years.
#' #' @param area list containing information on the 2D spacial dimension: list(total_area, n_sites, unit)
#' #' *total_area* the covered area by the points
#' #' Accepted units are square meter (m2) and square kilometer (km2) \code{check_space()$area}
#' #' @param crs Coordinate Reference Systems, as string and PROJ.4 format. Default is,
#' #' WGS 84 -- WGS84 - World Geodetic System 1984 crs="+proj=longlat +datum=WGS84 +no_defs"
#'
#' #' @param cost_function list of cost_function(s) used to calculate the cost distances between sites.
#' #' Depends on type and other methods used to calculate the cost distances.
#' #' @param geodynamic boolean for if sites change location or disappear over time
#' #' ,default is NULL and deduces from NA over time. If false, only one cost_distance is calculated and used
#' #' i.e. whencost_distances are the same, i.e. no geodynamic changes in the space
#' #' @param type_spec list containing type specific information, e.g. res (resolution) for raster and h3 spaces
#' #' @param author string with the name of the author. Default is NULL which gives the system user name obtained from \code{Sys.info()["user"]}
#' #' This is far from ideal, but it is better than nothing. Please fill this up and even consider leaving a contact information
#' #' @param source string with the source of the data, ideally should contain a publication with DOI and a valid URL
#' #' Default is list(env="missing", methods="missing")
#' #' @param description list with "environment" and "methods" containing information on the environmental data,
#' #' methods used and other relevant information such as source to raw data. Default is list("missing"
#' #' The "env" should describe the environmental data used, including it's units
#' #'
#' #' @return an informal and empty gen3sis_space object
#' #' @export
#' #'
#' #' @examples
#' create_space <- function(env=list(NA),
#'                          type="raster",
#'                          duration=list(from=NA, to=NA, by=NA, unit="Ma"),
#'                          area=list(total_area=NA, n_sites=NA, unit="km2"),
#'                          crs="+proj=longlat +datum=WGS84 +no_defs",
#'                          cost_function=list(NA),
#'                          geodynamic=NULL,
#'                          type_spec=list("res"=NA),
#'                          author="missing",
#'                          source="missing",
#'                          description=list(environment="missing", methods="missing")
#'                          ){
#'   space <- list()
#'   space$"env" <- env
#'   space$meta <- list()
#'   space$meta$"type" <- type
#'   space$meta$"type_spec" <- type_spec
#'   space$meta$"duration" <- duration
#'   space$meta$"area" <- area
#'   space$meta$"crs" <- crs
#'   space$meta$"cost_function" <- cost_function
#'   if(is.null(geodynamic)&!is.null(dim(env[[1]]))){
#'     geodynamic <- is_geodynamic(env)
#'   }
#'   space$meta$"geodynamic" <- geodynamic
#'   if (is.null(author)){
#'     author <- Sys.info()["user"]
#'     # remove attributes
#'     attributes(author) <- NULL
#'   }
#'   space$meta$"author" <- author
#'   space$meta$"source" <- source
#'   space$meta$"description" <- description
#'   return(invisible(space))
#' }
#'
#'
#'
#' #' Check gen3sis_space
#' #'
#' #' @param gen3sis_space either a gen3sis_space object to be checked, or NULL
#' #'
#' #' @return If a gen3sis_space object is provided, either a stop with printed error
#' #' report or a pass statement.
#' #' If gen3sis_space=NULL, this function returns lists of accepted values categories
#' #' according to \code{check_space()}
#' #' @export
#' #'
#' #' @examples
#' check_space <- function(gen3sis_space=NULL){
#'
#'   accepted <- list()
#'   accepted[["type"]] <- c("raster", "points", "h3")
#'   dur_units <- c("day", "wk", "mon", "yr", "dec", "cen", "mil", "Ma")
#'   accepted[["duration"]] <- dur_units[dur_units%in%measurements::conv_unit_options$duration]
#'   # area_units <- "km2"
#'   # accepted[["area"]]<- area_units[area_units%in%measurements::conv_unit_options$area]
#'
#'   if (is.null(gen3sis_space)){
#'     return(accepted)
#'   }
#'
#'   error_report <- NULL
#'   sp_ref <- create_space()
#'   #sp_ref <- gen3sis2::create_spaces()
#'
#'   error_report <- check_names(reference="env", datags=gen3sis_space, error_report)
#'   for (n_i in names(sp_ref$meta)){
#'     # n_i <- names(sp_ref$meta)[6]
#'     error_report <- check_names(reference=n_i, datags=gen3sis_space$meta, error_report)
#'     n_sub_e <- names(sp_ref$meta[[n_i]])
#'     if (length(n_sub_e)>1){
#'       for (s_i in n_sub_e){
#'         # s_i <- n_sub_e[1]
#'         error_report <- check_names(reference=s_i, datags=gen3sis_space$meta[[n_i]], error_report)
#'       }
#'     }
#'   }
#'   # stop in case of errors reported
#'   if (!is.null(error_report)){
#'     stop(error_report)
#'   }
#'   return("gen3sis_space [OK]")
#' }
#'
#' #' Check names of gen3sis_space
#' #'
#' #' @param reference string with the variable name to be tested, e.g. type, env
#' #' @param datags list of which \code{names(datags) is contrastet to reference}
#' #' @param error_report an error report that increases in case mismatches are found.
#' #' default is NULL
#' #'
#' #' @return
#' #' @export
#' check_names <- function(reference, datags, error_report=NULL){
#'   target=names(datags)
#'   if (!reference%in%target) { # if variables is missing
#'     error_report <- paste(
#'       error_report,
#'       (paste("The following space data is missing:", reference)),
#'       "\n")
#'   } # end if var name is missing
#'
#'   if (reference=="env"){ # check env
#'     if (!is.list(datags[[reference]])){ # if env is not a list
#'       error_report <- paste(
#'         error_report,
#'         (paste0("! >", reference, "< has to be a list of environmental variable(s)")),
#'         "\n")
#'     }
#'     # check if NAs are the same
#'     mask_NAs <- lapply(datags$env, function(x){
#'       is.na(x[,!colnames(x)%in%c("x","y"), drop=FALSE])
#'     })
#'     for (env_i in names(mask_NAs)){
#'       # env_i <- names(mask_NAs)[2]
#'       if (!identical(mask_NAs[[1]], mask_NAs[[env_i]])){
#'
#'         error_report <- paste(
#'           error_report,
#'           (paste0("! >please check your env data. NAs must match for all env's: \n
#'          i.e. environmental variables stored in env as a list.")),
#'           "\n")
#'       }
#'     } # end NA comparison loop
#'   } else if (any(is.na(datags[[reference]]))){ # if there is NA
#'     error_report <- paste(
#'       error_report,
#'       (paste0("! >", reference, "< can not be NA! please specify")),
#'       "\n")
#'   } # end if NA# end if env
#'
#'   return(error_report)
#' }
#'
#'
#' #' Prepare input and output directories
#' #'
#' #' @param dir_input "path to dir_input"
#' #' @param dir_output "path to dir_output"
#' #'
#' #' @return no value is returned
#' #' @export
#' #'
#' #' @examples
#' prepare_dirs <- function(dir_input, dir_output){
#'   if(!dir.exists(dir_input)){
#'     stop(paste("Input directory does not exist:", dir_input))
#'   }
#'   cat(paste0("Input directory found: \n [", dir_input, "] \n"))
#'   dir.create(dir_output, showWarnings = FALSE)
#'   cat(paste0("Output directory: \n [", dir_output, "] \n"))
#' }
#'
#'
