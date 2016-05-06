.onAttach <- function(libname, pkgname){

    v <- packageVersion(sprintf("%s", pkgname))
    msg0 <- "\t********************************************************\n"
    msg1 <- sprintf("\tCurrent version: %s\n", v)
    msg2 <- "\tThis version may contain important changes.\n"
    msg3 <- sprintf("\tUse news(Version == '%s', package = 'rCGH').\n", v)

    packageStartupMessage(c(msg0, msg1, msg2, msg3, msg0))

    invisible()

}
