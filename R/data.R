#' Time series of hospitalisations over time with current strategy.
#'
#' A dataset containing the dates and hospitalisation numbers for the current strategy.
#'
#' @format A data frame with 306 rows and 2 variables:
#' \describe{
#'   \item{dates}{\code{string}, date}
#'   \item{hospitalisations}{\code{int}, hopitalisation cases}
#' }
"cases_current_strategy"

#' Time series of hospitalisations over time with current strategy including infection-curve flattening.
#'
#' A dataset containing the dates and hospitalisation numbers for the current strategy with infection-curve flattening.
#'
#' @format A data frame with 306 rows and 2 variables:
#' \describe{
#'   \item{dates}{\code{string}, date}
#'   \item{hospitalisations}{\code{int}, hopitalisation cases}
#' }
"cases_current_strategy_flattened"

#' Time series of hospitalisations over time with a 'do nothing' scenario.
#'
#' A dataset containing the dates and hospitalisation numbers for a 'do nothing' scenario.
#'
#' @format A data frame with 306 rows and 2 variables:
#' \describe{
#'   \item{dates}{\code{string}, date}
#'   \item{hospitalisations}{\code{int}, hopitalisation cases}
#' }
"cases_do_nothing"
