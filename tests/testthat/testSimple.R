#'
#' Simple unit tests
#'

library(VolumeTaper)

sample_data <- data.table::data.table(
   species_code = c("PICEMAR","PICEGLA"),
   treeID = c(1,2),
   dbh_cm = c(30,30),
   Htot_m = c(18,18)
)

V_ObsH.R <- volTaper(data = sample_data, limit.dob = 9, limit.height = NULL)
test_that("Testing predictions with observed heights", {
  expect_equal(nrow(V_ObsH.R), 2)
  expect_equal(V_ObsH.R$Vol_m3[1], 0.5635172, tolerance = 1E-6)
  expect_equal(V_ObsH.R$Vol_m3[2], 0.5741721, tolerance = 1E-6)
})


V_EstH.R <- volTaper(data = sample_data, forceHeightEstimation = TRUE, limit.dob = 9, limit.height = NULL)
test_that("Testing predictions without observed heights", {
  expect_equal(nrow(V_EstH.R), 2)
  expect_equal(V_EstH.R$Vol_m3[1], 0.6079802, tolerance = 1E-6)
  expect_equal(V_EstH.R$Vol_m3[2], 0.6741257, tolerance = 1E-6)
})

speciesList <- getSpeciesList()
test_that("Testing species list", {
  expect_equal(nrow(speciesList), 34)
})

test_that("Testing example for random effect predictions", {
  expect_equal(nrow(VolumeTaper::exampleRandomEffectPrediction), 1868)
})

o <- suppressWarnings(predictRandomEffects(VolumeTaper::exampleRandomEffectPrediction))
# acerSac <- o$ACERSAH
# hist(acerSac[which(startsWith(acerSac$subject,"tree")),"blups"])
# mean(acerSac[which(startsWith(acerSac$subject,"tree")),"blups"]^2)^.5
# acerRub <- o$ACERRUB
# hist(acerRub[which(startsWith(acerRub$subject,"tree")),"blups"])
# mean(acerRub[which(startsWith(acerRub$subject,"tree")),"blups"]^2)^.5
# betupap <- o$BETUPAP
# hist(betupap[which(startsWith(betupap$subject,"tree")),"blups"])
# mean(betupap[which(startsWith(betupap$subject,"tree")),"blups"]^2)^.5
# betuall <- o$BETUALL
# hist(betuall[which(startsWith(betuall$subject,"tree")),"blups"])
# mean(betuall[which(startsWith(betuall$subject,"tree")),"blups"]^2)^.5

test_that("Testing blups", {
  expect_equal(length(o), 4)
  expect_equal(o$ACERSAH[1,"blups"],-0.02873648, tolerance=1e-6)
  expect_equal(o$ACERRUB[1,"blups"],0.0005296406, tolerance=1e-6)
  expect_equal(o$BETUPAP[1,"blups"],-0.08022338, tolerance=1e-6)
  expect_equal(o$BETUALL[1,"blups"],-0.07736193, tolerance=1e-6)
})

