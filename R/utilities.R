glibc_version <- function() {
    if (.Platform$OS.type!="unix" || !grepl("gnu",R.version$platform))
        return(NA)
    ss <- system("ldd --version",intern=TRUE)
    lnum <- as.package_version(gsub("^.* ([0-9.]+)$","\\1",ss[1]))
    class(lnum) <- c("package_version","numeric_version")
    lnum
}

check_rforge <- function()
    !is.na(gg <- glibc_version()) && gg < "2.14"
