glibc_version <- function() {
    if (.Platform$OS.type!="unix" || !grepl("gnu",R.version$platform))
        return(NA)
    ss <- system("ldd --version",intern=TRUE)
    lnum <- as.package_version(gsub("^.* ([0-9.]+)$","\\1",ss[1]))
    llist <- list(lnum)
    class(llist) <- c("package_version","numeric_version")
    llist
}

check_rforge <- function()
    !is.na(gg <- glibc_version()) && gg < as.package_version("2.14")
