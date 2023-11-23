#' Create empty gen3sis_space
#'
#' @param type string with type of gen3sis space. Accepted values are \code{check_space()$type}
#' @param duration list containing information on temporal dimension list(from, to, by, unit)
#' *from* is the oldest time-step; negative number if starting in the past, zero if starting in the present
#' *to* CAN ONLY BE smaller than *from*, since *to* is the latest time
#' *by* is the time interval increment. Note that this is constant and can only be positive;
#' *unit* is the time unit used. Accepted units are \code{check_space()$duration}
#' e.g. list(-20, 0, 1, "Ma") has a landscape that covers the last 20 Ma until the present, every 1 Ma.
#'      list(-800, 300, 10, "mil") has a landscape that covers the last 800 kyra or mil for millions of years
#'      and goes until the future 300 kya at every 100 mil years.
#' @param area list containing information on the 2D spacial dimension: list(total_area, n_sites, unit)
#' *total_area* the covered area by the points
#' unit is the time unit used. Accepted units are \code{check_space()$area}
#' @param crs Coordinate Reference Systems, as string
#' @param env Named list of environmental variables with X,Y, and time-steps from
#' the most recent of future time-step
#' @param cost_function List of cost_function(s) used
#' @param geo_dynamic Boolean, default is TRUE, should only be false
#' if cost_distances are the same, i.e. no geodynamic changes in the space
#'
#' @return a gen3sis_space object
#' @export
#'
#' @examples
create_space <- function(type="raster",
                         duration=list(from=NA, to=NA, by=NA, unit="Ma"),
                         area=list(total_area=NA, n_sites=NA, unit="km2"),
                         crs="+proj=longlat +datum=WGS84 +no_defs",
                         env=list(NA),
                         cost_function=list(NA),
                         geo_dynamic=TRUE){
  space <- list()
  # see convert_units from measurements, i.e.conv_units ?
  space[["type"]] <- type
  space[["duration"]] <- duration
  space[["area"]] <- area
  space[["crs"]] <- crs
  space[["env"]] <- env
  space[["cost_function"]] <- cost_function
  space[["geo_dynamic"]] <- geo_dynamic
  class(space) <- "gen3sis_space"
  return(invisible(space))
}



#' Check gen3sis_space
#'
#' @param gen3sis_space either a gen3sis_space object to be checked, or NULL
#'
#' @return If a gen3sis_space object is provided, either a stop with printed error
#' report or a pass statement.
#' If gen3sis_space=NULL, this function returns lists of accepted values cathegories
#' according to \code{check_space()}
#' @export
#'
#' @examples
check_space <- function(gen3sis_space=NULL){

  accepted <- list()
  accepted[["type"]] <- c("raster", "points", "h3")
  accepted[["duration"]] <- measurements::conv_unit_options$duration
  accepted[["area"]]<- measurements::conv_unit_options$area

  if (is.null(gen3sis_space)){
    return(accepted)
  }

  error_report <- NULL
  sp_ref <- create_space()

  for (n_i in names(sp_ref)){
    # n_i <- names(sp_ref)[5]
    error_report <- check_names(reference=n_i, datags=gen3sis_space, error_report)
    n_sub_e <- names(sp_ref[[n_i]])
    if (length(n_sub_e)>1){
      for (s_i in n_sub_e){
        # s_i <- n_sub_e[1]
        error_report <- check_names(reference=s_i, datags=gen3sis_space[[n_i]], error_report)
      }
    }
  }
  # stop in case of errors reported
  if (!is.null(error_report)){
    stop(error_report)
  }
  return("gen3sis_space [OK]")
}

#' Title
#'
#' @param reference string with the variable name to be tested, e.g. type, env
#' @param datags list of which \code{names(datags) is contrastet to reference}
#' @param error_report an error report that increases in case mismatches are found.
#' default is NULL
#'
#' @return
#' @export
check_names <- function(reference, datags, error_report=NULL){
  target=names(datags)
  if (!reference%in%target) { # if variables is missing
    error_report <- paste(
      error_report,
      (paste("The following space data is missing:", reference)),
      "\n")
  } # end if var name is missing

  if (reference=="env"){ # check env
    if (!is.list(datags[[reference]])){
      error_report <- paste(
        error_report,
        (paste0("! >", reference, "< has to be a list of environmental variable(s)")),
        "\n")
    }
    mask_NAs <- lapply(datags$env, function(x){
      is.na(x[,!colnames(x)%in%c("x","y"), drop=FALSE])
    })
    for (env_i in names(mask_NAs)[-1]){
      # env_i <- names(mask_NAs)[2]
      if (!identical(mask_NAs[[1]], mask_NAs[[env_i]])){

        error_report <- paste(
          error_report,
          (paste0("! >please check your env data. NAs must match for all env's: \n
         i.e. environmental variables stored in env as a list.")),
          "\n")
      }
    } # end NA comparison loop
  } else if (is.na(datags[[reference]])){ # if there is NA
    error_report <- paste(
      error_report,
      (paste0("! >", reference, "< can not be NA! please specify")),
      "\n")
  } # end if NA# end if env

  return(error_report)
}


#' Prepare input and output directories
#'
#' @param dir_input "path to dir_input"
#' @param dir_output "path to dir_output"
#'
#' @return no value is returned
#' @export
#'
#' @examples
prepare_dirs <- function(dir_input, dir_output){
  if(!dir.exists(dir_input)){
    stop(paste("Input directory does not exist:", dir_input))
  }
  cat(paste0("Input directory found: \n [", dir_input, "] \n"))
  dir.create(dir_output, showWarnings = FALSE)
  cat(paste0("Output directory: \n [", dir_output, "] \n"))
}

