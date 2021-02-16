.onDetach <- function(libpath) {
}

.onAttach <- function(lib, pkg) {
}

.onLoad<-function(libname, pkgname){
    library.dynam("EmulatorAnalysis", pkgname, libname)
}

.onUnload<-function(libpath){
   library.dynam.unload('EmulatorAnalysis', libpath)
}
