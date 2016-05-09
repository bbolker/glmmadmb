\name{admbControl}
\alias{admbControl}
\title{
  Control options for ADMB
}
\description{
  Determines options for running AD Model Builder
}
\usage{
admbControl(impSamp=0,maxfn=500,imaxfn=500,maxph=5,noinit=TRUE,shess=TRUE,
  run=TRUE, ZI_kluge=FALSE, poiss_prob_bound=TRUE)
}
\arguments{
  \item{impSamp}{
    Number of importance sampling steps
  }
  \item{maxfn}{
    Maximum number of evaluation steps
  }
  \item{imaxfn}{number of function evaluations used in intermediate
    optimization steps.}
  \item{maxph}{maximum phase}
  \item{noinit}{
    Don't initialize random effects before inner optimization: see ADMB-RE reference manual
  }
  \item{shess}{
    Use sparse Hessian: see ADMB-RE reference manual
  }
  \item{ZI_kluge}{retain a zero-inflation term of 0.001 even when \code{zeroInflation=FALSE} (for backward
    compatibility/testing; deprecated)}
  \item{poiss_prob_bound}{add a small term to prevent estimated Poisson probabilities of zero}
  \item{run}{Whether to run ADMB at all. This will write out files to the appropriate working directory, \emph{or} read existing
    output files in the working directory
    (\code{glmmadmb.\{std,bar,par,cor\}}) into a \code{glmmadmb} object.}
}
\section{Troubleshooting}{
  Some useful tricks for diagnosing and overcoming fitting problems:
  \itemize{
    \item{Diagnostics
      \itemize{
      \item to find the version of \code{glmmADMB} you are running, use
      \code{\link{sessionInfo}}
      \item to find the specific version of the ADMB binaries built into the package, use \code{\link{get_bin_version}}
      \item to see exactly where \code{glmmADMB} is looking for the
      executable on your system, and to see what command-line arguments
      it is passing to the ADMB executable, use \code{debug=TRUE}
      \item to see lots more output from ADMB, use \code{verbose=TRUE}
    }
    \item{Fixing problems
      \itemize{
	\item if the fit fails with "couldn't find STD file",
	try setting \code{noinit=FALSE} and/or \code{shess=FALSE}
	\item if you want to try running the ADMB model outside of R, use
	\code{run=FALSE} and set \code{save.dir} in the \code{glmmADMB}
	arguments.  This will result in a saved directory containing the
	\code{glmmadmb} executable (on systems other than Windows) as well as \code{glmmadmb.pin} and
	\code{glmmadmb.dat} files.  You can then run the \code{glmmadmb}
	executable with appropriate command-line arguments in a command
	shell or terminal window.  If you run the same R command again
	(i.e. still with \code{run=FALSE} and with the same \code{save.dir})
	\code{glmmADMB} will find any output files that have been produced
	and read them into a \code{glmmadmb} object.  You can also try sending
	the \code{glmmadmb.pin} and \code{glmmadmb.dat} files to an ADMB
	expert for diagnosis.
	\item if the ADMB output shows a statement like "Use the [] command line option ...",
        use the \code{extra.args} argument to \code{\link{glmmadmb}} to set command line
        options for ADMB. ADMB's diagnostic statements are usually fairly specific about
        which command line options to use, but see the ADMB reference manuals (linked from Details below)
        for more information. 
      }
    }
}
}
}
\details{
  See the AD Model Builder reference manuals at \url{http://admb-project.org/documentation/manuals/admb-user-manuals}
}
\value{
  Returns a list of options suitable for passing as the \code{admb.opts}
  argument to \code{\link{glmmadmb}}
}
\author{
  Ben Bolker
}
\examples{
admbControl(impSamp=200)
\dontrun{
## running to produce saved files in directory:
data(bacteria,package="MASS")
bacteria$present <- as.numeric(bacteria$y)-1
## run to generate files
glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria,
         save.dir="tmp",admb.opts=admbControl(run=FALSE),
         debug=TRUE)
## now run glmmadmb outside of R ...
## ./glmmadmb -maxfn 500 -maxph 5 -noinit -shess
## run to read in data
result <- glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                     family = "binomial", data = bacteria,
         save.dir="tmp",admb.opts=admbControl(run=FALSE))
}
}
\keyword{misc}
