.onAttach = function(libname, pkgname){
   if (interactive()) {
     packageStartupMessage(
#"-----------------------------------------------------------------------      
#Genotyping Uncertainty with Sequencing data for linkage MAPping (GUSMap)
#Copyright (C) 2017-2020 Timothy P. Bilton
#This program comes with ABSOLUTELY NO WARRANTY.
#This is free software, and you are welcome to redistribute it
#under certain conditions.
#------------------------------------------------------------------------
"Welcome to GUSMap v2.2.0"
     )}
}

.onUnload <- function (libpath) {
  library.dynam.unload("GUSMap", libpath)
}
