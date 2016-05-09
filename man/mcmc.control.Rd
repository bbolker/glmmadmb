\name{mcmcControl}
\alias{mcmcControl}
\title{
  Control options for MCMC after ADMB fitting
}
\description{
  Determines the options (number of steps, save interval, etc.)
  for running MCMC based on the estimated mode (maximum likelihood
  estimate) and parameter variance-covariance matrix
}
\usage{
mcmcControl(mcmc = 1000, mcmc2=0, mcsave, mcnoscale = FALSE, mcgrope = FALSE, mcmult = 1)
}
\arguments{
  \item{mcmc}{
    Total number of MCMC steps
  }
  \item{mcmc2}{
    MCMC2 steps (see ADMB-RE manual)
  }
  \item{mcsave}{
    Thinning interval.  Default is \code{pmax(1,floor(mcmc/1000))},
    i.e. aim to save 1000 steps 
  }
  \item{mcnoscale}{
    
  }
  \item{mcgrope}{
    
}
  \item{mcmult}{
    
  }
}
\details{
  See the AD Model Builder reference manual
}
\value{
  Returns a list of options suitable for passing as the \code{mcmc.opts}
  argument to \code{\link{glmmadmb}}
}
\author{
  Ben Bolker
}
\note{
  Some options (\code{mcmc2}, etc.) that can be used in AD Model Builder
  and ADMB-RE may not be available
}
\note{The functions in the \code{coda} and \code{plotMCMC} packages
  are useful for diagnosing convergence and other problems, summarizing,
  and displaying the results of MCMC runs.
  
  The tools for working with MCMC output in \code{glmmADMB} are 
limited at the moment (sorry!) --- in particular, the fixed-effect
parameter estimates are given in terms of the internally fitted
variables (which uses an orthogonalized version of the original
design matrix), not the original coefficients. If you need
to use the MCMC output, please contact the maintainers and
encourage them to work on them some more.
}
\examples{
mcmcControl(mcmc=2000)
}
\keyword{misc}
