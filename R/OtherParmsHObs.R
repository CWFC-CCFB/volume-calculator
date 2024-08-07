#'
#' A data.frame with the Other Parameters
#'
#' @format a data table
#' \describe{
#' \item{species_code}{A 6-character code for the species (See getSpeciesList function)}
#' \item{fixed}{The fixed effect term of the model}
#' \item{std_prov}{Standard deviation of the province random effect}
#' \item{std_plot}{Standard deviation of the plot random effect}
#' \item{std_tree}{Standard deviation of the tree random effect}
#' \item{std_res}{Standard deviation of the residual error term}
#' \item{corr}{The correlation parameter estimate}
#' \item{movaverage}{The moving average parameter estimate}
#' \item{varFunc}{The type of variance function}
#' \item{AB}{Alberta specific variance parameter estimate}
#' \item{BC}{British Columbia specific variance parameter estimate}
#' \item{MB}{Manitoba specific variance parameter estimate}
#' \item{ON}{Ontario specific variance parameter estimate}
#' \item{QC}{Quebec specific variance parameter estimate}
#' \item{SK}{Saskatchewan specific variance parameter estimate}
#' \item{YT}{Yukon specific variance parameter estimate}
#' \item{all}{All-provinces variance parameter estimate}
#' }
"OtherParmsHObs"
