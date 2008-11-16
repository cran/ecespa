.First.lib <- function (lib, pkg)
{
    library.dynam("ecespa", pkg, lib)
}
.Last.lib <- function (libpath){
    .isMethodsDispatchOn(FALSE)
}
