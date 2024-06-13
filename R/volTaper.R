#' Wood volume calculation using taper models
#'
#' @param data a data frame that includes at least 4 columns: species_code, treeID, dbh_cm, Htot_m
#' @param measuredH logic, TRUE if Htot is available .
#' @param limit.dob numeric, the limit of dob (cm) .
#' @param limit.height numeric, the limit of height (m) .
#'
#' @import dplyr
#' @import tidyr
#' @importFrom plyr join
#' @rawNamespace import(data.table, except = c(last, first, between))
#'
#'
#' @examples
#' sample_data <- data.table::data.table(
#'   species_code = c("PICEMAR","PICEGLA"),
#'   treeID = c(1,2),
#'   dbh_cm = c(30,30),
#'   Htot_m = c(18,18)
#' )
#' V_ObsH.R <- volTaper(data = sample_data, measuredH = TRUE, limit.dob = 9, limit.height = NULL)
#' V_EstH.R <- volTaper(data = sample_data, measuredH = FALSE, limit.dob = 9, limit.height = NULL)
#' @export
#'
#' @author Xiao Jing Guo \email{xiaojing.guo@NRCan-RNCan.gc.ca}
#' @references Ung, C.-H.; Guo, X.J.; Fortin, M. 2013. Canadian national taper models. For. Chron. 89(2):211-224
volTaper <- function(data , forceHeightEstimation = FALSE, limit.dob = NULL, limit.height = NULL){
  species_code <- dob2P_top <- h_top <- Vi <- dob2P <- NULL
  treeID <- c <- dbh_cm <- h <- cf <- Dob2P0 <- NULL
  measuredH <- !forceHeightEstimation & "Htot_m" %in% colnames(data) & all(!is.na(data[,"Htot_m"]))

  var_b <- var_a2 <- Htot_m <- Der2_b <- Der2_a2 <- NULL
  stddev_prov_b <- stddev_idPlot_b <- stddev_idtree_b <- stddev_prov_a2 <- stddev_idPlot_a2 <- a1 <- a2 <- b <- NULL

  if (measuredH) {
    parms <- VolumeTaper::Parms.ObsH
    parms[is.na(parms)] <- 0
    parms[,var_b := stddev_prov_b^2 + stddev_idPlot_b^2 + stddev_idtree_b^2 ]
    parmsToKeep <- parms[,c("species_code", "fixed_b", "var_b")]
  } else {
    parms <- VolumeTaper::Parms.EstH
    colnames(parms)[which(colnames(parms) == "b")] <- "fixed_b" ### TODO this is a patch the file should be resaved with the proper filename MF20240613
    parms[is.na(parms)] <- 0
    parms[,var_b := stddev_prov_b^2 + stddev_idPlot_b^2 + stddev_idtree_b^2 ]
    parms[,var_a2 := stddev_prov_a2^2 + stddev_idPlot_a2^2 ]
    parmsToKeep <- parms[,c("species_code", "a1", "a2", "fixed_b" , "var_a2", "var_b")]
    data <- data[parms[, c("species_code", "a1", "a2")], on = c("species_code"), nomatch = 0]
    data[, Htot_m := round(a1*dbh_cm^a2, 1)]
    data[, c("a1","a2") := NULL]
  }

  # expand the data
  df <- data %>%
    rowwise() %>%
    transmute(species_code,treeID, dbh_cm, Htot_m, h = list(seq_a(0.1, Htot_m, 0.1))) %>%
    unnest_longer(h)

  df2 <- data.table::setDT(df)[parmsToKeep, on = c("species_code"), nomatch = 0]
  df2 [, c:= (Htot_m - h)/(Htot_m - 1.3)]
  df2[  , Dob2P0 := dbh_cm^2* c * (h/1.3)^(2-fixed_b) ] # prediction conditional on the mode of the random effects
  df2[  , Der2_b := Dob2P0*(log(h/1.3))^2] # second-order derivatives with respect to random effects
  if (measuredH) {
    df2[  , cf := 0.5*var_b*Der2_b] # correction factor (see Eq.5 in Ung et al. 2013)
  } else {
    df2[,Der2_a2 := -1*dbh_cm^2*(h/1.3)^(2-fixed_b)*(log(dbh_cm))^2*(h-1.3)*a1*dbh_cm^a2*( a1*dbh_cm^a2+1.3)/( a1*dbh_cm^a2 - 1.3)^3]
    df2[  , cf := 0.5*var_b*Der2_b + 0.5*var_a2*Der2_a2] # correction factor (see Eq.5 in Ung et al. 2013)
  }
  df2[  , dob2P := Dob2P0+cf]
  df2 <- df2[dob2P >= 0]

  df2 <- df2[, c("species_code", "treeID", "h", "dob2P")]

  # volume
  data.table::setorder(df2, species_code, treeID, -h)
  df2[, dob2P_top:= data.table::shift(dob2P), by = .(species_code, treeID)]
  df2[, h_top:= data.table::shift(h), by = .(species_code, treeID)]
  # h criteria
  if(!(is.null(limit.height))) df2 <- df2[h_top <= limit.height]
  # dob criteria
  if(!is.null(limit.dob)) df2 <- df2[dob2P_top >= limit.dob^2]

  df2 <- df2[!is.na(dob2P_top)]
  df2[, Vi := (dob2P + dob2P_top)*0.1*acos(0)*2/80000]
  v_tree <- df2[, .(Vol_m3 = sum(Vi)), by = .(species_code, treeID)]
  if (measuredH == TRUE) equation <- "ObsH" else equation <- "EstH"
  # setnames(v_tree, "Vol_m3", paste0("vol_", equation))
  if (is.null(limit.dob)) limit.dob <- NA_real_
  if (is.null(limit.height)) limit.height <- NA_real_
  v_tree[, c("limit.dob", "limit.height", "Equation") := .(limit.dob, limit.height, equation)]
  data.v <- join(data,v_tree, by = c("species_code", "treeID"))
  return(data.v)
}
#
# getDiameters <- function(df, parms) {
#   df2 <- data.table::setDT(df)[Parms.ObsH [, c("species_code", "fixed_b", "var_b")], on = c("species_code"), nomatch = 0]
#
#   df2 [, c:= (Htot_m - h)/(Htot_m - 1.3)]
#
#   df2[  , Dob2P0 := dbh_cm^2* c * (h/1.3)^(2-fixed_b) ] # prediction conditional on the mode of the random effects
#   df2[  , Der2_b := Dob2P0*(log(h/1.3))^2] # second-order derivatives with respect to random effects
#
#   df2[  , cf := 0.5*var_b*Der2_b] # correction factor (see Eq.5 in Ung et al. 2013)
#
#   df2[  , dob2P := Dob2P0+cf]
#
# }


seq_a <- function(ifirst, ilast, iby){
  if (ilast*10 > floor(ilast*10))
    return(c(seq(ifirst, ilast, iby), ilast)) else
      return(seq(ifirst, ilast, iby))

}

