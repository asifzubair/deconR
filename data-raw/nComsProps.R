## code to prepare `nComsProps` dataset goes here

theDataDir <- "~/projects/decon/data/"
nComsProps <- read.csv(paste(theDataDir, "ncomms9971-s2.csv", sep=""), as.is=T)
usethis::use_data(nComsProps)
