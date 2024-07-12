
possibleProvinces <- c("AB", "BC", "MB", "ON", "QC", "SK", "YT")

#'
#' Wood volume calculation using taper models
#'
#' @param data a data frame that includes at least 4 columns: species_code, treeID, dbh_cm, Htot_m
#' @param forceHeightEstimation logic, TRUE to force the height estimation even if Htot is available.
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
#' V_ObsH.R <- volTaper(data = sample_data, limit.dob = 9, limit.height = NULL)
#' V_EstH.R <- volTaper(data = sample_data, forceHeightEstimation = TRUE, limit.dob = 9, limit.height = NULL)
#' @export
#'
#' @author Xiao Jing Guo \email{xiaojing.guo@NRCan-RNCan.gc.ca}
#' @references Ung, C.-H.; Guo, X.J.; Fortin, M. 2013. Canadian national taper models. For. Chron. 89(2):211-224
volTaper <- function(data, forceHeightEstimation = FALSE, limit.dob = NULL, limit.height = NULL) {
  species_code <- dob2P_top <- h_top <- Vi <- dob2P <- NULL
  treeID <- c <- dbh_cm <- h <- cf <- Dob2P0 <- NULL
  measuredH <- !forceHeightEstimation & "Htot_m" %in% colnames(data) & all(!is.na(data[,"Htot_m"]))

  var_b <- var_a2 <- Htot_m <- Der2_b <- Der2_a2 <- NULL
  stddev_prov_b <- stddev_idPlot_b <- stddev_idtree_b <- stddev_prov_a2 <- stddev_idPlot_a2 <- a1 <- a2 <- b <- NULL

  o <- prepareDataAndParms(data, measuredH)
  data <- o$data
  parms <- o$parms

  # expand the data
  df <- data %>%
    rowwise() %>%
    transmute(species_code,treeID, dbh_cm, Htot_m, h = list(seq_a(0.1, Htot_m, 0.1))) %>%
    unnest_longer(h)

  df2 <- getSquaredDiametersOverBarkCm2(df, parms, measuredH)

  # volume
  data.table::setorder(df2, species_code, treeID, -h)
  df2[, dob2P_top:= data.table::shift(dob2P), by = .(species_code, treeID)]
  df2[, h_top:= data.table::shift(h), by = .(species_code, treeID)]
  # h criteria
  if(!(is.null(limit.height))) df2 <- df2[h_top <= limit.height]
  # dob criteria
  if(!is.null(limit.dob)) df2 <- df2[dob2P_top >= limit.dob^2]

  df2 <- df2[!is.na(dob2P_top)]
  df2[, Vi := (dob2P + dob2P_top)*0.1*pi/80000] # 0.1 is the section length
  v_tree <- df2[, .(Vol_m3 = sum(Vi)), by = .(species_code, treeID)]
  if (measuredH == TRUE) equation <- "ObsH" else equation <- "EstH"
  # setnames(v_tree, "Vol_m3", paste0("vol_", equation))
  if (is.null(limit.dob)) limit.dob <- NA_real_
  if (is.null(limit.height)) limit.height <- NA_real_
  v_tree[, c("limit.dob", "limit.height", "Equation") := .(limit.dob, limit.height, equation)]
  data.v <- join(data,v_tree, by = c("species_code", "treeID"))
  return(data.v)
}

#'
#' Provide the Species List
#'
#' The values in the fiels `species_code` must be used to select
#' the appropriate species.
#'
#' @export
getSpeciesList <- function() {
  return(VolumeTaper::Species.List)
}

#
# A private function to select the appropriate set of parameter estimates and
# eventually add two fields in the data set in case of missing height.
#
prepareDataAndParms <- function(data, measuredH) {
  if (measuredH) {
    parms <- VolumeTaper::Parms.ObsH
    parms[is.na(parms)] <- 0
    parms[,var_b := stddev_prov_b^2 + stddev_idPlot_b^2 + stddev_idtree_b^2 ]
#    parmsToKeep <- parms[,c("species_code", "fixed_b", "var_b")]
  } else {
    parms <- VolumeTaper::Parms.EstH
    colnames(parms)[which(colnames(parms) == "b")] <- "fixed_b" ### TODO this is a patch the file should be resaved with the proper filename MF20240613
    parms[is.na(parms)] <- 0
    parms[,var_b := stddev_prov_b^2 + stddev_idPlot_b^2 + stddev_idtree_b^2 ]
    parms[,var_a2 := stddev_prov_a2^2 + stddev_idPlot_a2^2 ]
#    parmsToKeep <- parms[,c("species_code", "a1", "a2", "fixed_b" , "var_a2", "var_b")]
    data <- data[parms[, c("species_code", "a1", "a2")], on = c("species_code"), nomatch = 0]
    data[, Htot_m := round(a1*dbh_cm^a2, 1)]
    data[, c("a1","a2") := NULL]
  }
  return(list(parms = parms, data = data))
}



#
# Private function to produce the squared dbh prediction conditional on the mode of the random effects.
#
# @param df expanded data
# @param parms parameter estimates
# @param includeFirstDerivative a logical to indicate whether the first derivative should be
#  included in the output data.frame (false by default)
#
getConditionalSquaredDiameterOverBarkCm2 <- function(df, parms, includeFirstDerivative = F) {
  df2 <- data.table::setDT(df)[parms, on = c("species_code"), nomatch = 0]
  df2[, c:= (Htot_m - h)/(Htot_m - 1.3)]
  df2[, Dob2P0 := dbh_cm^2* c * (h/1.3)^(2-fixed_b) ] # prediction conditional on the mode of the random effects
  if (includeFirstDerivative) {
    df2[, Der1_b := -Dob2P0*log(h/1.3)] # first-order derivative with respect to random effects
  }
  return(df2)
}






#
# Private function to produce the corrected squared dbh prediction.
#
# @param df expanded data
# @param parms parameter estimates
# @param measuredHAvailable a logical to indicate whether the observed height is available
# @param applyCorrection a logical to indicate whether the correction should be applied (true by default)
#
# @references Eq. 5 in Ung et al. (2013)
#
getSquaredDiametersOverBarkCm2 <- function(df, parms, measuredHAvailable, applyCorrection = T) {
  df2 <- getConditionalSquaredDiameterOverBarkCm2(df, parms)
  # df2 <- data.table::setDT(df)[parms, on = c("species_code"), nomatch = 0]
  # df2[, c:= (Htot_m - h)/(Htot_m - 1.3)]
  # df2[, Dob2P0 := dbh_cm^2* c * (h/1.3)^(2-fixed_b) ] # prediction conditional on the mode of the random effects
  if (applyCorrection) {
    df2[, Der2_b := Dob2P0*(log(h/1.3))^2] # second-order derivatives with respect to random effects
    if (measuredHAvailable) {
      df2[, cf := 0.5*var_b*Der2_b] # correction factor (see Eq.5 in Ung et al. 2013)
    } else {          # else we need to account for the height-related parameters in the correction
      df2[,Der2_a2 := -1*dbh_cm^2*(h/1.3)^(2-fixed_b)*(log(dbh_cm))^2*(h-1.3)*a1*dbh_cm^a2*( a1*dbh_cm^a2+1.3)/( a1*dbh_cm^a2 - 1.3)^3]
      df2[, cf := 0.5*var_b*Der2_b + 0.5*var_a2*Der2_a2] # correction factor (see Eq.5 in Ung et al. 2013)
    }
    df2[, dob2P := Dob2P0+cf]
  } else {
    df2[, dob2P := Dob2P0]
  }
  df2 <- df2[dob2P >= 0]

  df2 <- df2[, c("species_code", "treeID", "h", "dob2P")]
  return(df2)
}


seq_a <- function(ifirst, ilast, iby){
  if (ilast*10 > floor(ilast*10))
    return(c(seq(ifirst, ilast, iby), ilast)) else
      return(seq(ifirst, ilast, iby))

}


#'
#' Predict the random effects
#'
#' @param data a data.frame object properly formatted
#'
#' @export
#'
predictRandomEffects <- function(data) {
  requiredFields <- colnames(VolumeTaper::exampleRandomEffectPrediction)
  for (f in requiredFields) {
    if (!(f %in% colnames(data))) {
      stop(paste("At least one field is missing in the data.",
                 paste("The data should contain the following fields:", paste(requiredFields, collapse = ", ")),
                 "See the `exampleRandomEffectPrediction` data.frame for an example.", sep="\n"))
    }
  }

  for (prov in unique(data$provinceID)) {
    if (!(prov %in% possibleProvinces)) {
      stop("The province ", prov, " is not recognized!\n Accepted province codes are: ", paste(possibleProvinces, collapse = ", "))
    }
  }

  screenedData <- NULL
  speciesInData <- unique(data$species_code)
  for (s in speciesInData) {
    if (!(s %in% getSpeciesList()$species_code)) {
      warning("Species ",s, " is not recognized and will be discarded!")
    } else {
      screenedData <- rbind(screenedData, data[which(data$species_code == s),])
    }
  }

  o <- prepareDataAndParms(screenedData, measuredH = T)
  parms <- o$parms
  screenedData <- o$data
  screenedData <- getConditionalSquaredDiameterOverBarkCm2(screenedData, parms, includeFirstDerivative = T)
  outputList <- list()
  for (s in unique(screenedData$species_code)) {
    message("Processing species ", s)
    data.s <- screenedData[which(screenedData$species_code == s),]
    data.s <- data.s[order(data.s$provinceID, data.s$plotID, data.s$treeID, data.s$h),]
    otherParms <- VolumeTaper::OtherParmsHObs
    otherParms <- as.data.frame(otherParms[which(otherParms$species_code == s),])
#    print(otherParms)
    species.blups <- NULL
    first <- F # change to true to have the first rMat and zGzt matrices displayed
    for (p in unique(data.s$plotID)) {
      data.ps <- data.s[which(data.s$plotID == p),]
      data.ps <- as.data.frame(data.ps)
      treeList <- unique(data.ps$treeID)
      zList <- list()
      zList[["plot"]] <- constructZMatrix(data.ps, unique(data.ps$plotID), data.ps$plotID, data.ps[1,]$stddev_idPlot_b + data.ps[1,]$stddev_prov_b)
      zList[["tree"]] <- constructZMatrix(data.ps, treeList, data.ps$treeID, data.ps[1,]$stddev_idtree_b)
      variances <- c(data.s[1,]$stddev_idPlot_b^2 + data.ps[1,]$stddev_prov_b^2, data.s[1,]$stddev_idtree_b^2)
      levels <- c("plot", "tree")

      gMat <- NULL
      for (i in 1:length(levels)) {
        l <- levels[i]
        zLevel <- zList[[l]]
        if (!is.null(zLevel)) {
          colnames(zList[[l]]) <- paste(l,colnames(zList[[l]]),sep="_")
          gMat <- c(gMat, rep(variances[i], ncol(zLevel)))
        }
      }

      zMat <- cbind(zList[["plot"]], zList[["tree"]])
      gMat <- diag(gMat, ncol = length(gMat), nrow = length(gMat))

      stdMat <- createStdMatrix(otherParms, data.ps, treeList)
      corrMat <- createCorrMatrix(otherParms, data.ps, treeList)
      rMat <- stdMat %*% corrMat %*% stdMat
      if (first) {
        message("rMat")
        print(rMat)
      }
      zGzt <- zMat %*% gMat %*% t(zMat)
      if (first) {
        message("zGzt")
        print(zGzt)
      }
      first <- F
      vMat <- zGzt + rMat
      invV <- solve(vMat)
      r <- data.ps$d_cm^2 - data.ps$Dob2P0
      blups <- gMat %*% t(zMat) %*% invV %*% r
      species.blups <- rbind(species.blups, data.frame(subject = colnames(zMat), blups))
    }

    outputList[[s]] <- species.blups
  }

  return(outputList)
}



#
# Private function to produce the correlation part of the R matrix
# @param otherParms the other parameter estimates for observed height models
# @param data the data to be processed
# @param a vector containing the individual tree id
#
createCorrMatrix <- function(otherParms, data, treeList) {
  if (is.na(otherParms$movaverage)) {
    corrFunctionType <- "CAR1"   # continuous first-order autoregressive structure
  } else {
    corrFunctionType <- "ARMA11" # ARMA structure
  }

  corrMatrix <- matrix(0, ncol = nrow(data), nrow = nrow(data))
  for (t in treeList) {
    index <- which(data$treeID == t)
    nObs <- length(index)
    if (corrFunctionType == "CAR1") {
      distancesM <- data[index,]$h
      dMat <- matrix(distancesM, nrow = nObs, ncol = nObs)
      dMat <- abs(dMat - t(dMat))
      corrMatrix[index, index] <- otherParms$corr^dMat
    } else {
      distancesM <- index
      dMat <- matrix(distancesM, nrow = nObs, ncol = nObs)
      dMat <- abs(dMat - t(dMat)) - 1
      diag(dMat) <- 0
      ar1Part <- otherParms$corr^dMat
      ma1Part <- matrix(otherParms$movaverage, nrow = nObs, ncol = nObs)
      diag(ma1Part) <- 1
      corrMatrix[index, index] <- ma1Part * ar1Part
    }
  }
  return(corrMatrix)
}



#
# Private function to produce the std part of the R matrix
# @param otherParms the other parameter estimates for observed height models
# @param data the data to be processed
# @param a vector containing the individual tree id
#
createStdMatrix <- function(otherParms, data, treeList) {
  varianceFunctionType <- otherParms$varFunc
  if (!is.na(otherParms$all)) {
    varParm <- otherParms$all
  } else {
    provID <- data[1,"provinceID"]
    if (!(provID %in% possibleProvinces)) {
      stop("The province id ", provID, " does not have a variance parameter!")
    }
    varParm <- otherParms[,provID]
    if (is.na(varParm)) {
      stop("The variance parameter for province ", provID, " is not available!")
    }
  }
  stdVect <- rep(0,nrow(data))
  for (t in treeList) {
    index <- which(data$treeID == t)
    if (varianceFunctionType == "Power") {
      covar <- data[index,]$dbh_cm  ### MF20240712 Presumably this is the good covariate. According to Jing's email this might vary across the species.
      stdVect[index] <- otherParms$std_res * covar ^ varParm
    } else if (varianceFunctionType == "Exp") {
      covar <- data[index,]$dbh_cm
      stdVect[index] <- otherParms$std_res * exp(covar * varParm)
    }
  }
  return(diag(stdVect))
}


#
# Private function to construct the Z matrix for the different levels
#
# If the standard deviation of the random effect is 0, the function
# returns NULL. The column names are also set to the levels of this
# random effect.
#
# @param data the data properly formatted
# @param levelList a vector with the different levels
# @param field the field that contains the value of this level
# @param stdRandomEffect the standard deviation of this random effect
#
constructZMatrix <- function(data, levelList, field, stdRandomEffect) {
  if (stdRandomEffect != 0) {
    z_mat <- matrix(0, nrow = nrow(data), ncol = length(levelList))
    names <- c()
    for (i in 1:length(levelList)) {
      names <- c(names, as.character(levelList[i]))
      index <- which(field == levelList[i])
      z_mat[index, i] <- data[index,]$Der1_b
    }
    colnames(z_mat) <- names
    return(z_mat)
  } else {
    return(NULL)
  }
}

