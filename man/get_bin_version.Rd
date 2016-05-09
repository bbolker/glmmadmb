\name{get_bin_version}
\alias{get_bin_version}
\title{
  Report ADMB version
}
\description{
   Report AD Model Builder version for compiled executable being used 
}
\usage{
get_bin_version(file_name = "glmmadmb")
}
\arguments{
  \item{file_name}{  Name of executable to search for
}
}
\details{
There are three version numbers that could be
relevant to debugging glmmADMB problems.  The first
is the glmmADMB version (retrievable via \code{help(package="glmmADMB")}
or \code{\link{sessionInfo}}); the second, for development
versions of glmmADMB, is the SVN (Subversion) revision number;
the third is the version of the compiled ADMB binary included
in the package.  By running \code{glmmadmb -version} on this
binary, this function retrieves the third (binary) version.
}
\value{
   Prints the version information and returns it, invisibly,
as a character vector.
}
\author{
Ben Bolker
}
\examples{
  si <- sessionInfo()
  sig <- si$otherPkgs$glmmADMB
  sig$Version       ## glmmADMB version
  if (!check_rforge()) {
        get_bin_version() ## ADMB binary version
      }
}
\keyword{misc}
