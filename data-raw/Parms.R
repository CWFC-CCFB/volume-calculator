
library(data.table)
library(readxl)
library(stringr)
library(readr)


Parms.ObsH <- read_excel("P:\\Jing 2\\2007-07\\Taper Data\\online calculator\\Tbl Eq2.xls")
Parms.EstH <- read_excel("P:\\Jing 2\\2007-07\\Taper Data\\online calculator\\Tbl Eq3.xls")
# prepare species list
# vcdExtra::datasets("VolumeTaper")
spp1 <- read_excel("P:\\Jing 2\\2007-07\\Taper Data\\J data\\Presentation\\species name list.XLS",sheet = "UNFORMAT")
spp <- read_excel("P:\\Jing 2\\2007-07\\Taper Data\\J data\\Presentation\\species name list.XLS",sheet = "formatted", skip = 3, col_names = F)
names(spp1)
setDT(spp1)
setDT(spp)
names(spp) <- c("ENGLISH_NAME",    "FRENCH_NAME" ,    "SCIENTIFIC_NAME")

speciescanfi <- read_excel( 'P:\\Jing\\2006-04 ---\\_Project\\Huor\\GY age\\Ref tables\\intolerant species Ung.xls')
setDT(speciescanfi)

spp1[, name_eng.a:= paste0( toupper(str_sub(name_eng, 1,1)), tolower(str_sub(name_eng, 2,-1)))]
spp1[name_eng == "Tamarack", name_eng.a:= "Tamarack larch"]
spp <- spp[!is.na(ENGLISH_NAME)]
spp[, ENGLISH_NAME.a:= paste0( toupper(str_sub(ENGLISH_NAME, 1,1)), tolower(str_sub(ENGLISH_NAME, 2,-1)))]
spp <- merge(spp1[,-c("ENGLISH_NAME",    "FRENCH_NAME" ,    "SCIENTIFIC_NAME")] , spp, by.x = "name_eng.a", by.y = "ENGLISH_NAME.a", all = T)

spp[name_eng.a != ENGLISH_NAME]
setdiff(spp$name_eng, Parms.EstH$name_eng)
setdiff(Parms.EstH$name_eng, spp$name_eng)
names(spp)

speciescanfi[CANFI_CODE==3403]
speciescanfi[str_detect(toupper(speciescanfi$SCIENTIFIC_NAME), "POP")]
speciescanfi[CANFI_CODE==1208]
spp.ra <- data.table(name_eng.a = "Red ash", canfi_code3 = 3403, name_eng = "Red Ash", ENGLISH_NAME= "Red ash", FRENCH_NAME = "Frêne rouge" , SCIENTIFIC_NAME = "Fraxinus pennsylvanica Marsh.")
spp <- rbind(spp, spp.ra)
spp <- merge(spp,Parms.EstH[, "name_eng"], by = "name_eng")
spp[canfi_code3==9991, canfi_code3:= 1208]
spp <- merge(spp, speciescanfi[, c("CANFI_CODE", "NFI_CODE" )], by.x = "canfi_code3", by.y = "CANFI_CODE", all.x = T)
spp[, species_code:= paste0(str_sub(NFI_CODE,1,4), str_sub(NFI_CODE,6,8))]


names(speciescanfi)
speciescanfi[NFI_Genus=="CARY"]
spp[, .N, by = .(species_code)][N>1]

Parms.EstH<- merge(spp[, c("name_eng", "species_code")], Parms.EstH)
Parms.EstH[, name_eng:=NULL]
Parms.ObsH<- merge(spp[, c("name_eng", "species_code")], Parms.ObsH)
Parms.ObsH[, name_eng:=NULL]
Parms.ObsH[, c("ar1", "ma1"):=NULL]
names(Parms.EstH)
Parms.ObsH[, .N, by = .(species_code)][N>1]
Parms.EstH[, .N, by = .(species_code)][N>1]

names(spp)
Species.List<- spp[, c("species_code", "ENGLISH_NAME", "FRENCH_NAME", "SCIENTIFIC_NAME")]
setorder(Species, species_code)
rm(spp, spp1, speciescanfi, spp.ra)



readr::write_csv(Species.List, "data-raw/Species.List.csv")
usethis::use_data(Species.List)
usethis::use_data(Parms.EstH)
usethis::use_data(Parms.ObsH)
