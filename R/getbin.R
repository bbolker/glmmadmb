get_bbot_versions <- function(os="all",rev="latest",bits="all",
                              bburl="http://admb-project.org/buildbot/glmmadmb/") {
    OSvals <- c("fedora","macos","ubuntu","windows")
    z <- suppressWarnings(readLines(url(bburl)))
    desc <- z[grep("glmmadmb",z)]
    sizes <- z[grep("[0-9][BKM]",z)]
    desc <- desc[-1]  ## header line
    desc <- gsub(".+\"(([a-z0-9.]|-)+).*","\\1",desc)
    ## assume 0B = 0 -- everything else will be in M
    sizes <- as.numeric(gsub(".+>([0-9.]+).*","\\1",sizes))
    d <- data.frame(desc,sizes)
    d <- subset(d,sizes>0)
    if (os!="all") d <- d[grep(os,d$desc),]
    if (bits!="all") d <- d[grep(paste(d$bits,"bit",sep="")),]
    if (rev=="latest") {
        ## d$type <- gsub("-r[0-9]+.*$","",d$desc)
        d$type <- gsub("\\.bin$","",gsub("(glmmadmb-|r[0-9]+-)","",d$desc))
        d$ver <- as.numeric(gsub(".*-r([0-9]+).*$","\\1",d$desc))
        d <- droplevels(ddply(d, .(type), function(x) x[which.max(x$ver),]))
    } else if (ver!="all") {
        d <- d[d$ver==ver,]
    }
    OSopts <- paste(OSvals,collapse="|")
    ## suppress false-positive R CMD check warnings
    type <- ver <- sdesc <- . <- fullOS <- NULL
    d <- mutate(d,
                desc=as.character(desc),
                type=as.character(type),
                sdesc=gsub("(glmmadmb|r[0-9]+|\\.bin|\\.exe)","",desc),
                glmmadmbver=as.numeric(gsub("[-r]","",
                str_extract(desc,"-r[0-9]+-"))),
                fullOS=str_extract(sdesc,paste0("(",OSopts,")[.0-9]+")),
                OSver=str_extract(fullOS,"[0-9.]+"),
                OSstr=str_extract(fullOS,"[[:alpha:]]+"))
    d <- mutate(d,sdesc=gsub("^-+","",gsub("--","-",sdesc)))
    return(d)
}

download_bin <- function(x,quiet=FALSE,
                         bburl="http://admb-project.org/buildbot/glmmadmb/") {
  download.file(paste(bburl,x,sep=""),x)
}

