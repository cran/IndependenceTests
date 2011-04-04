.First.lib <- function(lib, pkg){
  library.dynam("IndependenceTests", pkg, lib)
  x <- read.dcf(file = system.file("DESCRIPTION", package = "IndependenceTests"))
  cat("\n")
  write.dcf(x)
  cat("\n")
}
