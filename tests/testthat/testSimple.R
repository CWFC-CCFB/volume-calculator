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
