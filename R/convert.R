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
space_raster_to_h3 <- function(dir_input,
                               dir_output=file.path(dir_input, "h3_v2"),
                               res=0,
                               verbose=0,
                               compute_distances=TRUE,
                               agg_fun = mean,
                               reports = FALSE
){
  # variables lexicon
  # o_gls_all: original raster-spaces.rds
  # o_cd_fl: full distances files
  # numb_ts: number of timesteps

  # Set up dirs
  gen3sis2:::prepare_dirs(dir_input, dir_output)
  if (reports) {
    report_path <- file.path(dir_output,"reports")
    dir.create(report_path)
  }

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
  time_steps <- colnames(o_gls_all$env[[1]][,-c(1,2)])
  n_ts <- length(time_steps)

  if (o_gls_all$meta$geodynamic & n_ts!=length(o_cd_fl)){
    #stop function and show error message
    stop("Mismatch on number of time steps and full distances.rds")
  }

  # get points and coordinates
  all_pts <- o_gls_all$env[[1]][,c("x", "y")]
  all_pts_sf <- sf::st_as_sf(as.data.frame(all_pts), coords = c("x","y"), crs = o_gls_all$meta$crs)
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

  # Creates a dictionary of which raster cells are inside each h3 cell
  h3_cell_dict <- list()
  for (h3_cll in unique(h3_cell_idx)) {
    point_cell <- cbind(all_pts,h3_cell_idx) |> as.data.frame()
    l <- list(row.names(point_cell[point_cell$h3_cell_idx==h3_cll,]))
    names(l) <- h3_cll
    h3_cell_dict <- append(h3_cell_dict, l)
  }

  # saves the dictionary as data.frame
  if (reports) {
    temp_df <- do.call(rbind, lapply(names(h3_cell_dict), function(h3) {
      data.frame(h3_cell = h3, raster_cell = h3_cell_dict[[h3]], stringsAsFactors = FALSE)
    }))

    write.csv(temp_df,file.path(report_path,"raster_h3_dictionary.txt"), row.names = FALSE)
  }

  # reduces the dict to a single representative cell based on distance to the H3 centroid
  for (h3_cll in names(h3_cell_dict)) {
    cent_h3 <- h3jsr::cell_to_point(h3_cll) |>
      sf::st_coordinates()

    dist_from_cent <- lapply(h3_cell_dict[[h3_cll]], function(x){
      geosphere::distGeo(cent_h3, all_pts[as.numeric(x),])
    }) |>
      unlist()

    h3_cell_dict[[h3_cll]] <- h3_cell_dict[[h3_cll]][which(dist_from_cent == min(dist_from_cent))]
  }

  # saves the representative cell
  if (reports) {
    temp_df <- do.call(rbind, lapply(names(h3_cell_dict), function(h3) {
      data.frame(h3_cell = h3, raster_cell = h3_cell_dict[[h3]], stringsAsFactors = FALSE)
    }))

    write.csv(temp_df, file.path(report_path,"representative_cells.txt"), row.names = FALSE)
  }

  # Project values from raster to h3 aggregating when necessary
  envar <- names(o_gls_all$env)
  envar_superlist <- list()
  for (vari in envar) {
    # vari <- envar[[1]]
    vari_df <- o_gls_all$env[[vari]] |> as.data.frame()
    vari_df$h3_cell <- h3_cell_idx

    agg_df <- list()
    for (t_s in time_steps) {
      # t_s <- time_steps[[1]]
      ts_df <- vari_df[,c("x","y","h3_cell",t_s)]
      colnames(ts_df)[ncol(ts_df)] <- paste0("place_holder_",colnames(ts_df)[ncol(ts_df)])
      agg_formula <- paste(paste0("place_holder_",t_s), "~", "h3_cell") |> stats::as.formula()
      agg_df_var <- stats::aggregate(agg_formula, data = ts_df, FUN = agg_fun)
      agg_df_var <- merge(ts_df, agg_df_var, by = "h3_cell", all.x = TRUE, suffixes = c("_raster","_h3"))
      agg_df <- append(agg_df,list(agg_df_var))
    }

    merged_final_left <- base::Reduce(function(df1, df2) {
      base::merge(df1, df2, by = c("h3_cell","x","y"), all.x = T)
    }, agg_df)

    vari_coords <- vari_df[,c("x","y")]
    vari_coords$original_order <- 1:nrow(vari_coords)

    vari_df <- merge(vari_coords, merged_final_left, by = c("x","y"), all.x = TRUE)

    vari_df <- vari_df[order(vari_df$original_order),]
    vari_df <- vari_df[,-c(which(colnames(vari_df)=="original_order"))]
    row.names(vari_df) <- NULL

    new_colnames <- colnames(vari_df[,-c(which(colnames(vari_df)%in%c("x","y","h3_cell")))])
    new_colnames <- gsub("place_holder_", "", new_colnames)
    names(vari_df)[which(!names(vari_df)%in%c("x","y","h3_cell"))] <- new_colnames

    # Create multiple list if agg_fun returns more than 1 value
    if (any(sapply(vari_df, is.matrix))){
      mtx_columns <- names(vari_df)[sapply(vari_df, is.matrix)]
      for (mtx_col in mtx_columns) {
        # mtx_col <- mtx_columns[[1]]
        temp_df <- as.data.frame(vari_df[,mtx_col])
        names(temp_df) <- paste0(names(temp_df),"--",mtx_col)
        vari_df <- vari_df[,-c(which(names(vari_df)==mtx_col))]
        vari_df <- cbind(vari_df,temp_df)
      }
    }

    l <- list(vari_df)
    names(l) <- vari

    envar_superlist <- append(envar_superlist, l)
  }

  # organize and plot everything
  envs_raster <- list()
  envs_h3 <- list()
  for (envar in names(envar_superlist)) {
    # envar <- names(envar_superlist)[[1]]
    raster_columns <- names(envar_superlist[[envar]])[grep("_raster$",names(envar_superlist[[envar]]))]
    raster_df <- envar_superlist[[envar]][,c("x","y",raster_columns)]
    names(raster_df) <- gsub("_raster$","",names(raster_df))
    l <- list(raster_df)
    names(l) <- envar
    envs_raster <- append(envs_raster, l)

    h3_columns <- names(envar_superlist[[envar]])[grep("_h3$",names(envar_superlist[[envar]]))]
    h3_df <- envar_superlist[[envar]][,c("h3_cell",h3_columns)]
    h3_df <- h3_df[!duplicated(h3_df$h3_cell),]

    names(h3_df) <- gsub("_h3$","",names(h3_df))

    h3_df <- cbind(
      h3jsr::cell_to_point(h3_df$h3_cell) |>
        sf::st_coordinates() |>
        as.data.frame(),
      h3_df
    )

    names(h3_df)[c(1,2)] <- c("x","y")

    new_rownames <- c()
    for (cll in h3_df$h3_cell) {
      new_rownames <- c(new_rownames, h3_cell_dict[[cll]])
    }

    rownames(h3_df) <- new_rownames

    if (ncol(h3_df) != (n_ts+3)) {
      var_columns <- setdiff(names(h3_df), c("x","y","h3_cell"))
      ts_versions <- unique(sub("--.*", "", var_columns))

      version_list <- lapply(ts_versions, function(v) {
        cols <- grep(paste0("^", v, "--"), names(h3_df), value = TRUE)

        df_sub <- h3_df[,c(c("x","y","h3_cell"), cols)]

        names(df_sub)[-1] <- sub(paste0("^", v, "--"), "", names(df_sub)[-1])

        return(df_sub)
      })

      if (verbose > 2 & envar == names(envar_superlist[1])) {
        plot(terra::rast(raster_df[,c("x","y",t_s)], type = "xyz"), main=paste("Raster in",t_s))

        h3jsr::cell_to_polygon(stats::na.omit(version_list[[1]][,c("h3_cell",t_s)])$h3_cell) |>
          plot(main = paste("H3 cells in", t_s))
      }

      names(version_list) <- paste0(envar,"-",ts_versions)

      version_list <- lapply(version_list, function(v){
        v[,names(v)!="h3_cell"]
      })

      envs_h3 <- append(envs_h3, version_list)
    } else {
      if (verbose > 2 & envar == names(envar_superlist[1])) {
        plot(terra::rast(raster_df[,c("x","y",t_s)], type = "xyz"), main=paste("Raster in",t_s))

        h3jsr::cell_to_polygon(na.omit(h3_df[,c("h3_cell",t_s)])$h3_cell) |>
          plot(main = paste("H3 cells in", t_s))
      }

      h3_df <- h3_df[,c(!names(h3_df)%in%c("h3_cell"))]

      l <- list(h3_df)
      names(l) <- envar
      envs_h3 <- append(envs_h3, l)
    }
  }

  # saves the original raster envs for comparision
  if (reports) {
    envs_raster <- lapply(envs_raster, function(x){
      x$raster_cell_idx <- rownames(x)
      x
    })

    saveRDS(envs_raster, file.path(report_path,"envs_raster.rds"))
  }

  # garbage colector to save ram
  rm(vari_df, envar_superlist, envs_raster)
  gc()

  # calculate the conversion error based on centroids distance
  h3_centroids <- h3jsr::cell_to_point(h3_cell_idx) |>
    sf::st_coordinates()

  error_m <- c()
  for (i in 1:nrow(h3_centroids)) {
    distance <- geosphere::distGeo(h3_centroids[i,], all_pts[i,])
    error_m[i] <- distance
  }

  if (verbose>1){
    hist(error_m/1000, main="Error in Km")
  }

  if (reports) {
    temp_df <- cbind(all_pts,h3_centroids,error_m)
    colnames(temp_df) <- c("raster_x","raster_y","h3_x","h3_y","error_m")

    utils::write.csv(temp_df, file.path(report_path,"conversion_error.txt"),row.names = FALSE)
  }

  # gc to save ram
  rm(h3_df, raster_df)
  gc()

  # Loop over timesteps to update distance matrix
  for (ti in 1:n_ts){
    # ti <- 1
    num_ts <- (n_ts-1):0
    if(verbose>0){
      cat("--\n")
      cat(paste(ti,"of ", n_ts, "time_steps:", time_steps[ti],"\n"))
    }

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
      cat(paste("Loaded:", rev(o_cd_fl)[ti],"\n"))

      new_cost_dist_full <- lfd_ti[row.names(lfd_ti) %in% row.names(envs_h3[[1]]), colnames(lfd_ti) %in% row.names(envs_h3[[1]])]

      # save cost function ti
      saveRDS(new_cost_dist_full, file.path(dir_output, "distances_full", tiis))
      if (verbose>0){
        cat(paste("Saved:", file.path(dir_output, "distances_full", tiis),"\n"))
      }
    }
  }

  # save final landscape
  final_space <- o_gls_all # copy old landscape as ref.
  final_space$env <- envs_h3
  final_space$meta$type <- "h3"
  final_space$meta$type_spec <- list("res"=res)
  final_space$meta$area$total_area <- sum(h3jsr::cell_area(h3_fill, final_space$area$unit, simple=TRUE))
  final_space$meta$area$n_sites <- as.numeric(length(unique(h3_cell_idx)))
  saveRDS(final_space, file.path(dir_output, "spaces.rds"))
  return(cat("Space converted to h3 and saved to [", dir_output,"]","\n"))
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

      # Copia o arquivo
      file.copy(from = file.path(dir_input, f),
                to   = file.path(dir_output, f),
                overwrite = TRUE)
    }
    space <- readRDS(file.path(dir_output, "spaces.rds"))
    space$meta$type <- "points"
    saveRDS(space, file.path(dir_output, "spaces.rds"))
  }
}

