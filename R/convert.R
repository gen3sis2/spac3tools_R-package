#' @title spac3tools: Tools for gen3sis2 input manipulation including conversion,
#' compression and decompression.

#' @description spac3tools is a package that provides tools for manipulating
#' gen3sis input data, including conversion of landscapes to spaces,
#' compression and decompression of spaces, and conversion of spaces from raster
#' to h3 format and much more.
#' @references Development team
#' @details This package is meant to provide support for
#' gen3sis2 previous and future generations, facilitating storage, modification
#' and reproducibility of the input data used in gen3sis and gen3sis2 simulations.
#' @seealso \code{\link{landscape_to_space}}   \code{\link[gen3sis2]{check_spaces}}  \code{\link{compress_space}}  \code{\link{decompress_space}}
#' @concept spatial tools for inputs used by gen3sis2 eco-evolutionary modeling engine
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
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#'
"_PACKAGE"





#' Convert a gen3sis landscape to a gen3sis2 space object
#'
#' This function reads a `landscapes.rds` file (in gen3sis input format) from the specified input directory,
#' converts it to a `spaces.rds` object compatible with gen3sis2, and saves the result in the output directory.
#' It ensures that time steps are properly ordered and that spatial metadata is derived automatically.
#'
#' @param dir_input Location of landscape.rds to be converted
#' @param dir_output location to store the converted space.rds
#' @param duration see \code{?gen3sis2::create_spaces}
#' @param crs see \code{?gen3sis2::create_spaces}
#' @param cost_function Cost function for connectivity: see \code{?create_spaces}. Defaults to a cost of 2 for sites with missing data (NA).
#' @param ... see \code{?gen3sis2::create_spaces}
#'
#' @importFrom gen3sis2 create_spaces check_spaces
#' @importFrom terra rast crs cellSize ext res
#'
#' @return Saves a `spaces.rds` file in the specified output directory.
#' @export
#'
#' @example inst/examples/conv_landscape_to_space_help.R
landscape_to_space <- function(dir_input=NA,
                               dir_output=dir_input,
                               duration=list(from=NA, to=NA, by=NA, unit="Ma"),
                               crs="+proj=longlat +datum=WGS84 +no_defs",
                               cost_function=list(xx=gcf$XXHarderNA_dist_Km),...){

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
                     cost_function = list(cost_function),
                     geodynamic=NULL,
                     type_spec=list("res"=terra::res(ex_r)),
                     ...
                     )
  gen3sis2::check_spaces(gs)
  saveRDS(gs, file.path(dir_output, "spaces.rds"), compress=T)
  print(paste0("space.rds type=", gs$type , " saved to [", file.path(dir_output, "spaces.rds"), "]"))
}


#' Convert a space.rds type raster to type h3
#'
#' The function projects the raster data in any H3 resolution. In the process, values aggregation (e.g., values mean) is often necessary.
#' The original distance matrices are updated. If more than one raster cell is inside the same H3 cell, the closest raster cell to the H3 centroid is choosen as the representative.
#' The distance between two H3 cells is assumed as the distance calculated between their representative raster cells.
#' It is recommended to use this on top of the smallest available raster resolution.
#'
#'@param dir_input Location of space.rds type raster
#'@param dir_output location to store the converted space.rds of type h3
#'@param res resolution of the h3 grid. This parameter is passed directly to the
#' `res` parameter of `hrjsr::get_children()`. For details refer to h3 documentation
#'@param verbose integer, 0 for no messages, 1 for some messages, 2 for all messages
#'@param compute_distances logical, if TRUE, the cost distances are computed and stored (default TRUE)
#'@param agg_fun function. An aggregation function to resume raster values in a single H3 cell
#'when necessary. Examples include:
#'\itemize{
#'  \item \code{mean} - mean of all values (default)
#'  \item \code{median} - median of all values
#'  \item \code{min} or \code{max} - minimum or maximum value
#'  \item \code{summary} - summary statistics including min, max, mean, media, and quartiles
#'}
#'@param reports logical. Default is FALSE. If TRUE, saves four files:
#'\itemize{
#'  \item envs_raster.rds - the original enviromental variables in raster spaces
#'  \item raster_h3_dictionary.txt - a data.frame containing H3 cells and the indexes of raster cells resumed in it
#'  \item representative_cells.txt - a data.frame containing H3 cells and the index of the chosen representative raster cell
#'  \item conversion_error.txt - a data.frame containing the centroids for each raster and its respective H3 cell and the distance between them
#'}
#'
#'@return Saves a `spaces.rds` file in the specified output directory of type "h3".
#'
#'@importFrom geosphere distGeo
#'@importFrom h3jsr polygon_to_cells get_disk point_to_cell cell_to_point cell_area
#'@importFrom sf st_as_sf st_transform st_convex_hull st_union st_coordinates
#'@importFrom terra plot rast
#'@importFrom graphics hist
#'@importFrom stats aggregate na.omit as.formula
#'@importFrom utils write.csv
#'
#'@export
#'@example inst/examples/conv_space_raster_to_h3_help.R
space_raster_to_h3 <- function(
    dir_input,
    dir_output = file.path(dir_input, "h3_v2"),
    res = 0,
    verbose = 0,
    compute_distances = TRUE,
    agg_fun = mean,
    reports = FALSE
) {
  # ---------------------------------------------------------------------------
  # 1. Prepare output directories and load the source space
  # ---------------------------------------------------------------------------

  gen3sis2:::prepare_dirs(dir_input, dir_output)

  if (reports) {
    report_path <- file.path(dir_output, "reports")
    dir.create(report_path, recursive = TRUE, showWarnings = FALSE)
  }

  source_space <- readRDS(file.path(dir_input, "spaces.rds"))

  time_steps <- colnames(source_space$env[[1]])[
    !colnames(source_space$env[[1]]) %in% c("x", "y")
  ]
  n_ts <- length(time_steps)

  # Distance files are only required when distance conversion is requested.
  distance_files <- character(0)

  if (compute_distances) {
    distance_dir <- file.path(dir_input, "distances_full")
    distance_files <- list.files(
      distance_dir,
      pattern = "^distances_full_.*\\.rds$",
      full.names = FALSE
    )

    if (!length(distance_files)) {
      stop("No full distance matrices found in: ", distance_dir)
    }

    distance_numbers <- as.numeric(
      gsub("[^0-9.-]+", "", distance_files)
    )
    distance_files <- distance_files[order(distance_numbers)]

    if (
      isTRUE(source_space$meta$geodynamic) &&
      n_ts != length(distance_files)
    ) {
      stop("Mismatch between environmental timesteps and distance files.")
    }

    dir.create(
      file.path(dir_output, "distances_full"),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Convert raster-cell coordinates to H3 membership
  # ---------------------------------------------------------------------------

  all_pts <- source_space$env[[1]][, c("x", "y"), drop = FALSE]
  raster_ids <- rownames(all_pts)
  # in case their are no row names (although unlikely)
  if (is.null(raster_ids)) {
    raster_ids <- as.character(seq_len(nrow(all_pts)))
    rownames(all_pts) <- raster_ids
  }

  all_pts_sf <- sf::st_as_sf(
    as.data.frame(all_pts),
    coords = c("x", "y"),
    crs = source_space$meta$crs
  )
  # the crs may not equal 4326, so convert to be sure
  all_pts_4326 <- sf::st_transform(all_pts_sf, 4326)
  all_pts_4326_xy <- sf::st_coordinates(all_pts_4326)
  rownames(all_pts_4326_xy) <- raster_ids

  h3_cell_idx <- h3jsr::point_to_cell(all_pts_4326, res = res)

  # Build the H3 domain once from the union of raster cells that contain habitat
  # in at least one timestep. The NA mask of the first environmental variable is
  # sufficient because gen3sis2 requires matching NA masks across variables.
  habitat_any_raster <- apply(
    source_space$env[[1]][, time_steps, drop = FALSE],
    1,
    function(values) any(!is.na(values))
  )

  # Sorting the H3 indexes makes the assignment deterministic: the same input
  # always produces the same numeric site ID. These IDs remain fixed even if
  # habitat appears or disappears through time.
  h3_ids <- sort(unique(
    as.character(h3_cell_idx[habitat_any_raster])
  ))

  if (!length(h3_ids)) {
    stop("No H3 cells contain habitat in any timestep.")
  }

  h3_coordinates <- h3jsr::cell_to_point(h3_ids) |>
    sf::st_coordinates()

  h3_cell_dictionary <- data.frame(
    cell_id = seq_along(h3_ids),
    h3_cell = h3_ids,
    x = h3_coordinates[, 1],
    y = h3_coordinates[, 2],
    stringsAsFactors = FALSE
  )

  # Character numeric IDs are used because gen3sis2 indexes matrices by row
  # names and some dispersal code explicitly coerces those names to numeric.
  h3_cell_dictionary$cell_id <- as.character(
    h3_cell_dictionary$cell_id
  )

  h3_to_cell_id <- stats::setNames(
    h3_cell_dictionary$cell_id,
    h3_cell_dictionary$h3_cell
  )

  cell_id_to_h3 <- stats::setNames(
    h3_cell_dictionary$h3_cell,
    h3_cell_dictionary$cell_id
  )

  # Save the permanent lookup regardless of reports. It is required to convert
  # model output back to real H3 indexes.
  utils::write.csv(
    h3_cell_dictionary,
    file.path(dir_output, "h3_cell_dictionary.csv"),
    row.names = FALSE
  )

  # Keep every raster member of each H3 cell. Representatives are selected from
  # these candidates separately for every distance timestep.
  h3_cell_members <- split(raster_ids, as.character(h3_cell_idx))
  h3_cell_members <- h3_cell_members[h3_ids]

  if (reports) {
    raster_h3_dictionary <- do.call(
      rbind,
      lapply(h3_ids, function(h3_id) {
        data.frame(
          cell_id = unname(h3_to_cell_id[h3_id]),
          h3_cell = h3_id,
          raster_cell = h3_cell_members[[h3_id]],
          stringsAsFactors = FALSE
        )
      })
    )

    utils::write.csv(
      raster_h3_dictionary,
      file.path(report_path, "raster_h3_dictionary.txt"),
      row.names = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 3. Rank candidate raster representatives by distance to each H3 centroid
  # ---------------------------------------------------------------------------

  h3_centroids <- h3_coordinates
  rownames(h3_centroids) <- h3_ids

  # CHANGED:
  # The original function compared WGS84 H3 centroids with all_pts in the
  # source CRS. Here both coordinate sets are explicitly in EPSG:4326.
  #
  # We store a ranked list rather than one representative. At each timestep,
  # the first candidate that exists in that timestep's distance matrix is used.
  h3_representative_rank <- lapply(names(h3_cell_members), function(h3_id) {
    candidates <- h3_cell_members[[h3_id]]

    candidate_xy <- all_pts_4326_xy[
      candidates,
      c("X", "Y"),
      drop = FALSE
    ]

    centroid_xy <- matrix(
      h3_centroids[h3_id, c("X", "Y")],
      nrow = length(candidates),
      ncol = 2,
      byrow = TRUE
    )

    candidate_distance <- geosphere::distGeo(
      centroid_xy,
      candidate_xy
    )

    candidates[order(candidate_distance)]
  })

  names(h3_representative_rank) <- names(h3_cell_members)

  if (reports) {
    static_representatives <- data.frame(
      cell_id = unname(h3_to_cell_id[names(h3_representative_rank)]),
      h3_cell = names(h3_representative_rank),
      raster_cell = vapply(
        h3_representative_rank,
        `[`,
        character(1),
        1
      ),
      stringsAsFactors = FALSE
    )

    utils::write.csv(
      static_representatives,
      file.path(report_path, "representative_cells_static.txt"),
      row.names = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 4. Aggregate every environmental variable from raster cells to H3 cells
  # ---------------------------------------------------------------------------

  environmental_superlist <- list()

  for (variable_name in names(source_space$env)) {
    variable_df <- as.data.frame(source_space$env[[variable_name]])
    variable_df$h3_cell <- h3_cell_idx

    aggregated_timesteps <- vector("list", length(time_steps))
    names(aggregated_timesteps) <- time_steps

    for (time_name in time_steps) {
      timestep_df <- variable_df[
        ,
        c("x", "y", "h3_cell", time_name),
        drop = FALSE
      ]

      names(timestep_df)[ncol(timestep_df)] <-
        paste0("place_holder_", time_name)

      aggregation_formula <- stats::as.formula(
        paste0("`place_holder_", time_name, "` ~ h3_cell")
      )

      aggregated_value <- stats::aggregate(
        aggregation_formula,
        data = timestep_df,
        FUN = agg_fun
      )

      aggregated_timesteps[[time_name]] <- merge(
        timestep_df,
        aggregated_value,
        by = "h3_cell",
        all.x = TRUE,
        suffixes = c("_raster", "_h3")
      )
    }

    merged_variable <- base::Reduce(
      function(left, right) {
        merge(
          left,
          right,
          by = c("h3_cell", "x", "y"),
          all.x = TRUE
        )
      },
      aggregated_timesteps
    )

    original_coordinates <- variable_df[, c("x", "y"), drop = FALSE]
    original_coordinates$original_order <- seq_len(nrow(original_coordinates))

    merged_variable <- merge(
      original_coordinates,
      merged_variable,
      by = c("x", "y"),
      all.x = TRUE
    )

    merged_variable <- merged_variable[
      order(merged_variable$original_order),
      ,
      drop = FALSE
    ]
    merged_variable$original_order <- NULL
    rownames(merged_variable) <- NULL

    value_columns <- !names(merged_variable) %in% c("x", "y", "h3_cell")
    names(merged_variable)[value_columns] <- gsub(
      "place_holder_",
      "",
      names(merged_variable)[value_columns]
    )

    # Preserve support for aggregation functions that return multiple values.
    if (any(vapply(merged_variable, is.matrix, logical(1)))) {
      matrix_columns <- names(merged_variable)[
        vapply(merged_variable, is.matrix, logical(1))
      ]

      for (matrix_column in matrix_columns) {
        expanded_column <- as.data.frame(merged_variable[, matrix_column])
        names(expanded_column) <- paste0(
          names(expanded_column),
          "--",
          matrix_column
        )

        merged_variable[[matrix_column]] <- NULL
        merged_variable <- cbind(merged_variable, expanded_column)
      }
    }

    environmental_superlist[[variable_name]] <- merged_variable
  }

  # ---------------------------------------------------------------------------
  # 5. Build the H3 environmental tables
  # ---------------------------------------------------------------------------

  envs_raster <- list()
  envs_h3 <- list()

  for (variable_name in names(environmental_superlist)) {
    variable_data <- environmental_superlist[[variable_name]]

    raster_columns <- grep(
      "_raster$",
      names(variable_data),
      value = TRUE
    )

    raster_df <- variable_data[
      ,
      c("x", "y", raster_columns),
      drop = FALSE
    ]
    names(raster_df) <- gsub("_raster$", "", names(raster_df))
    envs_raster[[variable_name]] <- raster_df

    h3_columns <- grep(
      "_h3$",
      names(variable_data),
      value = TRUE
    )

    h3_df <- variable_data[
      ,
      c("h3_cell", h3_columns),
      drop = FALSE
    ]
    h3_df <- h3_df[!duplicated(h3_df$h3_cell), , drop = FALSE]

    # Keep only H3 cells that contain habitat in at least one timestep. This is
    # the same permanent domain used by h3_cell_dictionary.csv.
    h3_df <- h3_df[
      as.character(h3_df$h3_cell) %in% h3_ids,
      ,
      drop = FALSE
    ]

    # Put every environmental table into permanent numeric-ID order.
    h3_df <- h3_df[
      match(h3_ids, as.character(h3_df$h3_cell)),
      ,
      drop = FALSE
    ]

    names(h3_df) <- gsub("_h3$", "", names(h3_df))

    h3_coordinates <- h3jsr::cell_to_point(h3_df$h3_cell) |>
      sf::st_coordinates() |>
      as.data.frame()

    names(h3_coordinates)[1:2] <- c("x", "y")
    h3_df <- cbind(h3_coordinates, h3_df)

    # CHANGED: assign the permanent numeric ID rather than the H3 string or a
    # representative raster ID. The H3 string remains recoverable through
    # h3_cell_dictionary.csv.
    rownames(h3_df) <- unname(h3_to_cell_id[h3_df$h3_cell])

    if (ncol(h3_df) != n_ts + 3) {
      variable_columns <- setdiff(
        names(h3_df),
        c("x", "y", "h3_cell")
      )
      versions <- unique(sub("--.*", "", variable_columns))

      version_list <- lapply(versions, function(version_name) {
        version_columns <- grep(
          paste0("^", version_name, "--"),
          names(h3_df),
          value = TRUE
        )

        version_df <- h3_df[
          ,
          c("x", "y", "h3_cell", version_columns),
          drop = FALSE
        ]

        names(version_df)[-(1:3)] <- sub(
          paste0("^", version_name, "--"),
          "",
          names(version_df)[-(1:3)]
        )

        rownames(version_df) <- unname(
          h3_to_cell_id[version_df$h3_cell]
        )
        version_df$h3_cell <- NULL
        version_df
      })

      names(version_list) <- paste0(variable_name, "-", versions)
      envs_h3 <- append(envs_h3, version_list)
    } else {
      rownames(h3_df) <- unname(h3_to_cell_id[h3_df$h3_cell])
      h3_df$h3_cell <- NULL
      envs_h3[[variable_name]] <- h3_df
    }
  }


  # Validate that all environmental variables use the same permanent site IDs.
  expected_cell_ids <- h3_cell_dictionary$cell_id
  env_id_ok <- vapply(
    envs_h3,
    function(environment) {
      identical(rownames(environment), expected_cell_ids)
    },
    logical(1)
  )

  if (!all(env_id_ok)) {
    stop(
      "Converted H3 environmental tables do not share the permanent ",
      "numeric cell-ID ordering."
    )
  }

  if (reports) {
    envs_raster <- lapply(envs_raster, function(x) {
      x$raster_cell_idx <- rownames(x)
      x
    })

    saveRDS(
      envs_raster,
      file.path(report_path, "envs_raster.rds")
    )
  }

  rm(environmental_superlist, envs_raster)
  gc()

  # ---------------------------------------------------------------------------
  # 6. Report raster-to-H3 centroid displacement
  # ---------------------------------------------------------------------------

  point_h3_centroids <- h3jsr::cell_to_point(h3_cell_idx) |>
    sf::st_coordinates()

  conversion_error_m <- vapply(
    seq_len(nrow(point_h3_centroids)),
    function(i) {
      geosphere::distGeo(
        point_h3_centroids[i, ],
        all_pts_4326_xy[i, ]
      )
    },
    numeric(1)
  )

  if (verbose > 1) {
    graphics::hist(
      conversion_error_m / 1000,
      main = "Raster point to H3 centroid error (km)"
    )
  }

  if (reports) {
    conversion_error <- cbind(
      as.data.frame(all_pts_4326_xy),
      as.data.frame(point_h3_centroids),
      error_m = conversion_error_m
    )

    names(conversion_error) <- c(
      "raster_lon",
      "raster_lat",
      "h3_lon",
      "h3_lat",
      "error_m"
    )

    utils::write.csv(
      conversion_error,
      file.path(report_path, "conversion_error.txt"),
      row.names = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 7. Convert each raster distance matrix
  # ---------------------------------------------------------------------------

  representative_log <- list()

  if (compute_distances) {
    for (ti in seq_len(n_ts)) {
      if (
        !isTRUE(source_space$meta$geodynamic) &&
        ti > 1
      ) {
        next
      }

      if (isTRUE(source_space$meta$geodynamic)) {
        distance_file <- rev(distance_files)[ti]
      } else {
        distance_file <- "distances_full_0.rds"
      }

      if (verbose > 0) {
        cat(
          "--\n",
          ti,
          "of",
          n_ts,
          "timesteps:",
          time_steps[ti],
          "\n"
        )
      }

      source_distance <- readRDS(
        file.path(
          dir_input,
          "distances_full",
          distance_file
        )
      )

      if (is.null(rownames(source_distance)) ||
          is.null(colnames(source_distance))) {
        stop(
          "Distance matrix has no row/column names: ",
          distance_file
        )
      }

      available_raster_ids <- intersect(
        rownames(source_distance),
        colnames(source_distance)
      )

      # CHANGED: derive the distance-matrix domain from actual H3 habitat at
      # this timestep, but retain the permanent numeric IDs assigned above.
      time_name <- time_steps[ti]

      if (!time_name %in% names(envs_h3[[1]])) {
        stop("Environmental timestep not found: ", time_name)
      }

      active_cell_ids <- rownames(envs_h3[[1]])[
        !is.na(envs_h3[[1]][[time_name]])
      ]

      active_h3 <- unname(cell_id_to_h3[active_cell_ids])

      # Select the closest member of each active H3 cell that is actually
      # available in the original raster distance matrix for this timestep.
      representatives <- vapply(
        active_h3,
        function(h3_id) {
          ranked_candidates <- h3_representative_rank[[h3_id]]
          available_candidates <- ranked_candidates[
            ranked_candidates %in% available_raster_ids
          ]

          if (!length(available_candidates)) {
            stop(
              "H3 habitat cell ", h3_id,
              " (static cell ID ", h3_to_cell_id[[h3_id]], ") has no ",
              "source raster cell in ", distance_file
            )
          }

          available_candidates[1]
        },
        character(1)
      )

      # Explicit indexing preserves the environment-table order.
      converted_distance <- source_distance[
        representatives,
        representatives,
        drop = FALSE
      ]

      # Rename rows and columns to the permanent numeric site IDs. These names
      # are identical to the corresponding environment and coordinate row names.
      rownames(converted_distance) <- active_cell_ids
      colnames(converted_distance) <- active_cell_ids

      if (!identical(rownames(converted_distance), active_cell_ids) ||
          !identical(colnames(converted_distance), active_cell_ids)) {
        stop(
          "Converted distance matrix IDs do not match environmental habitat ",
          "IDs for timestep ", time_name
        )
      }

      saveRDS(
        converted_distance,
        file.path(
          dir_output,
          "distances_full",
          distance_file
        )
      )

      representative_log[[distance_file]] <- data.frame(
        timestep = time_name,
        distance_file = distance_file,
        cell_id = active_cell_ids,
        h3_cell = active_h3,
        raster_cell = unname(representatives),
        stringsAsFactors = FALSE
      )

      if (verbose > 0) {
        cat(
          "Saved:",
          file.path(
            dir_output,
            "distances_full",
            distance_file
          ),
          "\n"
        )
      }
    }
  }

  if (reports && length(representative_log)) {
    utils::write.csv(
      do.call(rbind, representative_log),
      file.path(
        report_path,
        "representative_cells_by_timestep.txt"
      ),
      row.names = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 8. Save the final H3 space
  # ---------------------------------------------------------------------------

  final_space <- source_space
  final_space$env <- envs_h3
  final_space$meta$type <- "h3"
  final_space$meta$type_spec <- list(
    res = res,
    site_id = "numeric",
    h3_dictionary = "h3_cell_dictionary.csv"
  )

  # H3 cell centroids are returned in longitude/latitude. Keep metadata and
  # coordinates in the same CRS and recalculate the extent accordingly.
  final_space$meta$crs <- "EPSG:4326"
  final_space$meta$area$extent <- c(
    xmin = min(h3_cell_dictionary$x),
    xmax = max(h3_cell_dictionary$x),
    ymin = min(h3_cell_dictionary$y),
    ymax = max(h3_cell_dictionary$y)
  )

  # Metadata describes the permanent H3 domain. Habitat availability through
  # time remains represented by NA values in the environmental tables.
  final_space$meta$area$total_area <- sum(
    h3jsr::cell_area(
      h3_ids,
      final_space$meta$area$unit,
      simple = TRUE
    )
  )
  final_space$meta$area$n_sites <- length(h3_ids)

  saveRDS(
    final_space,
    file.path(dir_output, "spaces.rds")
  )

  if (verbose > 0) {
    cat(
      "Space converted to H3 and saved to [",
      dir_output,
      "]\n"
    )
  }

  invisible(dir_output)
}


#' Do local conversion of a space.rds type h3 to type points locally or not
#'
#' The function simply changes the type marker, leading \code{gen3sis2} to properly handle the spaces.rds
#'
#' @param dir_input Location of spaces.rds type h3 to be converted to points
#' @param dir_output Where to save the spaces.rds type points. Not necessary if inplace = TRUE. Default is NULL
#' @param inplace logical. If TRUE, conversion is local, overwriting the input spaces.rds, if FALSE, a new directory is created and all files copied to it prior conversion. If FALSE, a dir_output must be provided. Default is TRUE
#' @export
#' @example inst/examples/conv_space_h3_to_points_help.R
space_h3_to_points <- function(dir_input, dir_output=NULL, inplace = TRUE){
  if (inplace){
    space <- readRDS(file.path(dir_input, "spaces.rds"))
    space$meta$type <- "points"
    saveRDS(space, file.path(dir_input, "spaces.rds"))
  } else if (!inplace & !is.null(dir_output)){
    input_files <- list.files(dir_input, recursive = T, full.names = F)

    for (f in input_files) {
      dir_target <- file.path(dir_output, dirname(f))
      if (!dir.exists(dir_target)){
        dir.create(dir_target, recursive = TRUE)
      }

      # Copy the file
      file.copy(from = file.path(dir_input, f),
                to   = file.path(dir_output, f),
                overwrite = TRUE)
    }
    space <- readRDS(file.path(dir_output, "spaces.rds"))
    space$meta$type <- "points"
    saveRDS(space, file.path(dir_output, "spaces.rds"))
  }
}

