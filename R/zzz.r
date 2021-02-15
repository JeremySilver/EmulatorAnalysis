.First.lib <- function(lib, pkg) {
##   library.dynam("EmulatorAnalysis", pkg, lib)
##   packageStartupMessage("Use help(EmulatorAnalysis) for an overview of this library")
  cat('Running .First.lib \n')
}

.Last.lib <- function(lib, pkg) {
##   library.dynam("EmulatorAnalysis", pkg, lib)
##   packageStartupMessage("Use help(EmulatorAnalysis) for an overview of this library")
  cat('Running .Last.lib \n')
}

## .onLoad <- function(lib, pkg) {
##   library.dynam("EmulatorAnalysis", pkg, lib)
##   packageStartupMessage("Use help(EmulatorAnalysis) for an overview of this library")
##   cat('Running onLoad \n')
## }

.onUnload <- function(lib, pkg) {
##   library.dynam("EmulatorAnalysis", pkg, lib)
##   packageStartupMessage("Use help(EmulatorAnalysis) for an overview of this library")
  cat('Running onLoad \n')
}

.onAttach <- function(lib, pkg) {
    print(c(pkg,lib))
    # library.dynam("EmulatorAnalysis", pkg, lib)
    packageStartupMessage("Use help(EmulatorAnalysis) for an overview of this library")
    cat('Running onAttach \n')
}

## .onLoad<-function(libname, pkgname){
##   library.dynam("EmulatorAnalysis", pkg, lib)
##   packageStartupMessage("Use help(EmulatorAnalysis) for an overview of this library")
## }

## .onUnload<-function(libpath){
##   library.dynam.unload('EmulatorAnalysis', libpath)
## }
