# data

#' A H3 list containing one dataset from temperature values over 5 time-steps for a set of h3 cells
#'
#' @format A list with one element, a data.frame with 40 row and 8 columns:
#' \describe{
#'   \item{h3_address}{H3 cells addresses}
#'   \item{ts_1}{temperature value for the given H3 cell at the first time-step}
#'   \item{ts_2}{temperature value for the given H3 cell at the second time-step}
#'   \item{ts_3}{temperature value for the given H3 cell at the third time-step}
#'   \item{ts_4}{temperature value for the given H3 cell at the fourth time-step}
#'   \item{ts_5}{temperature value for the given H3 cell at the fifth time-step}
#'   \item{x}{the x coordinate for cells' centroid}
#'   \item{y}{the y coordinate for cells' centroid}
#' }
#' @examples
#' \dontrun{
#' require(spac3tools)
#' data("h3_list")
#' h3_list
#' }
"h3_list"
