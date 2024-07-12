#
# Compiling other parms for covariance modelling with observed height
# @author Mathieu Fortin - July 2024
#

rm(list = ls())

OtherParms <- read.csv("./data-raw/OtherParmsHObs.csv")
OtherParms2 <- read.csv("./data-raw/OtherParmsHObs2.csv")
OtherParms <- merge(OtherParms, OtherParms2, by="species_code")
colnames(OtherParms)[1] <- "ENGLISH_NAME"
speciesList <- as.data.frame(VolumeTaper::getSpeciesList())
OtherParms <- merge(speciesList[,c("species_code", "ENGLISH_NAME")], OtherParms, by="ENGLISH_NAME")
OtherParmsHObs <- OtherParms[,colnames(OtherParms)[which(colnames(OtherParms) != "ENGLISH_NAME")]]
save(OtherParmsHObs, file = "./data/OtherParmsHObs.rda")
