\name{getdata}
\alias{getdata}
\title{
  Retrieve package data
}
\description{
  Pulls in auxiliary data associated with the package
}
\usage{
getdata(s)
}
\arguments{
  \item{s}{
    The name of a data file (only .rda for now)
}
}
\details{
  Details of the way that R tests and saves data means that
  data that depend on the package being loaded for their definition
  (in particular, fitted \code{glmmADMB} model objects) cannot
  be treated as standard data.  This is a workaround.
}
\value{
  As with \code{\link{load}}, the values returned are the
  \emph{names} of the objects loaded.  As a side effect, the
  objects in the associated files will be loaded in the
  global workspace.
}
\author{
  Ben Bolker
}
\examples{
getdata("OwlModel")
## regular data:
data(package="stats")
}
\keyword{misc}
