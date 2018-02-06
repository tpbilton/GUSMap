.onLoad = function(libname, pkgname){
  if (interactive()) {
    packageStartupMessage(
      "Genotyping Uncertainty with Sequening data for linkage MAPping (GUSMap)
      Copyright (C) 2018 Timothy P. Bilton
      This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
      This is free software, and you are welcome to redistribute it
      under certain conditions; type `show c' for details."
  )}
}
