\name{glmmadmb}
\alias{glmmadmb}
\alias{glmm.admb}
\title{Generalized Linear Mixed Models Using AD Model Builder}
\description{
  Fits (extended)
  generalized linear mixed-effects models to data using a variety of
  distributions and link functions, including zero-inflated models.
}
\usage{
glmmadmb(formula, data, family="poisson", link,start,
          random, corStruct="diag", easyFlag=TRUE,
          zeroInflation=FALSE, 
          admb.opts=admbControl(),
          mcmc=FALSE,mcmc.opts=mcmcControl(),save.dir,
          verbose=FALSE,extra.args="",
          bin_loc=NULL,debug=FALSE)
glmm.admb(\dots)  ## deprecated
}
\arguments{
  \item{formula}{a two-sided linear formula describing the 
    part of the model. The response, which may be a numeric vector
    (for most models) or a matrix (for binomial models)
    is on the left of a \code{~} operator and the terms,
    separated by \code{+} operators, are on the
    right.  Random effects terms should be enclosed in parentheses
    \code{(x|g)}, where \code{x} represents the formula for the
    parameter(s) that vary across grouping variables and \code{g}
    represents the grouping variables; for example, \code{(1|g)} indicates
    a model where only intercepts vary across groups.
    \emph{Grouping variables must be factors.}
    The model may include an offset term, which must be specified on the
    scale of the linear predictor and included in \code{data}:
    see \code{\link{offset}}. Alternatively, the random effects terms may be specified
  separately using the \code{random} argument.}
  \item{data}{a data frame containing the variables named in \code{formula}.}
  \item{family}{a string determining the response distribution:
    currently implemented options (and default link functions) are
    \itemize{
      \item exponential family distributions: \code{"poisson"} (log),
      \code{"binomial"} (logit), \code{"gamma"} (log), \code{"gaussian"}
      (identity)
      \item negative binomial (log): \code{"nbinom"},
      \code{"nbinom1"} (negative binomial ``type 1'': see Details),
      \code{"nbinom2"} (a synonym for \code{"nbinom"})
      \item extended distributions: \code{"beta"} (logit),
      \code{"logistic"} (identity), \code{"betabinomial"} (logit)
      \item truncated distributions (log):
      \code{"truncpoiss"}, \code{"truncnbinom"}, \code{"truncnbinom1"}
      (these can be used to fit hurdle models)
    }
  }
  \item{link}{a string specifying the link function:
    currently implemented options are
    \code{"log"}, \code{"logit"}, \code{"probit"}, \code{"cloglog"},
    and \code{"identity"}.}
  \item{start}{an (optional) named list specifying starting values, with
    (optional) components
    \itemize{
    \item \code{fixed} (fixed-effect parameters: default 0);
    \item \code{pz} (zero-inflation probability: default 0.02);
    \item \code{log_alpha} (log of dispersion parameter: default 1);
    \item \code{RE_sd} (log-standard deviations of random effects: default 0.25);
    \item \code{RE_cor} (off-diagonal elements of Cholesky decomposition of
    random-effect correlation matrix: default 0.0001);
    \item \code{u} (random effects [unscaled, i.e. not multiplied by
    corresponding RE standard deviation]: default 0)}
  }
  \item{random}{a formula specifying the random effects.
    A single random effect may be specified as \code{~x|g},
    but parentheses must be used if more than one random
    effect is specified (e.g., \code{~(x1|g1)+(x2|g2)}).}
  \item{corStruct}{a string specifying the covariance structure of the
    random effects vector. Two types of covariance matrices are
    currently implemented: \code{"diag"} (diagonal matrix) and
    \code{"full"} (positive definite matrix with all elements being
    estimated).  At present the same covariance structure must be
    specified for all grouping variables in the model; at some point
    in the future this requirement may be relaxed.}
  \item{easyFlag}{(\emph{currently inactive}: whether a faster but less robust optimization
    algorithm should be employed (only \code{"poisson"} and
    \code{"nbinom"} families).}
  \item{zeroInflation}{whether a zero-inflated model should be fitted
    (only \code{"poisson"} and \code{"nbinom"} families).}
  \item{admb.opts}{options for AD Model Builder: see \code{\link{admbControl}}}
  \item{mcmc}{(logical) run MCMC after fitting, based on estimated mode
    and parameter variance-covariance matrix}
  \item{mcmc.opts}{list of options for MCMC run: see \code{\link{mcmcControl}}}
  \item{save.dir}{Optionally, a directory name where all ADMB output
    files are saved. If \code{save.dir} is not specified, a temporary
    directory is created (with \code{\link{tempfile}}) and all ADMB input
    and output files are deleted after the ADMB run finishes and the results are read into R.}
  \item{verbose}{(logical) Print interim output from AD Model Builder
    run?}
  \item{extra.args}{(character) extra arguments to AD Model Builder (see
    ADMB manual)}
  \item{bin_loc}{location of \code{glmmadmb} binary file (normally determined automatically by
              \code{glmmADMB:::get_bin_loc})}
  \item{debug}{(logical) print debugging statements?}
  \item{\dots}{arguments to pass through to \code{glmmadmb}}
}
\details{
  \itemize{
    \item For advice on troubleshooting, see \code{\link{admbControl}}.
    \item Parameterization of the negative binomial distribution:
    \itemize{
      \item for \code{family=="nbinom"}:
      \deqn{\mathrm{Var}(Y) = E(Y) \left[1 +
	\frac{E(Y)}{\alpha}\right]}{Var(Y) = E(Y) * (1 + E(Y)/alpha)}
      (i.e. Hardin and Hilbe NB2; \code{alpha} corresponds to the \code{size} parameter in
      \code{\link{dnbinom}})
      \item for \code{family=="nbinom1"}:
      \deqn{\mathrm{Var}(Y) = E(Y) \cdot \alpha}{Var(Y) = E(Y) * alpha}
      (i.e. Hardin and Hilbe NB1: \code{alpha} corresponds to the scale
      parameter in a quasi-Poisson model)
    }
    \item Parameterization of the beta-binomial distribution:
    as in Morris (1997), in terms of \eqn{N} (number of trials),
    \eqn{p} (per-trial probability, determined by the inverse-link
    function applied to the linear predictor), \eqn{alpha}
    (overdispersion parameter: as \eqn{alpha} becomes large,
    the distribution converges to binomial). The density function is
    \deqn{p(x) = \frac{C(N,x) \mbox{Beta}(N-x+\theta(1-p),x+\theta p)}%
      {\mbox{Beta}(\theta(1-p),\theta p)}%
    }{p(x) = (C(N,x)*Beta(N-x+theta*(1-p),x+theta*p))/Beta(theta*(1-p),theta*p)}.

    \item Zero-inflation: With probability \eqn{1-pz}, \eqn{Y} comes from a
    Poisson (or negative binomial) distribution, and with probability
    \eqn{pz}, \eqn{Y} is zero (Bohning et al. 1999).
    
    \item Parameters are estimated by maximum likelihood, using the Laplace
    approximation to evaluate the marginal likelihood.
    When \code{impSamp >0}, 
    importance sampling is used to improve the Laplace approximation
    (Skaug and Fournier 2006).
    \item Due to technical limitations, offset expressions can
    contain \emph{at most} one internal set of parentheses; i.e.
    \code{offset(log(x))} works, but \code{offset((log(x)))} will not.
    If necessary, define a new variable in the data frame, e.g.
    \code{mydata$myoffset <- with(mydata,log(3*exp(x)))}.
    \item Internally, glmmADMB works with an orthogonalized version of
    the fixed-effect design matrix.  The \code{phi} matrix returned in the results can
    be used to convert from orthogonalized to 'real' fixed parameters,
    i.e. \code{b_real == b_orth \%*\% phi} or \code{b_orth = b_real
      \%*\% solve(phi).}  Fixed effects specified in the starting values
    are automatically orthogonalized.
    \item If you are specifying starting values from a different GLMM
    fit for more than the fixed-effect parameters, you may want to use
    \code{extra.args="-phase 5"} to skip over the estimation phase where
    ADMB fits the model without random effects.
  }
}
\value{
  An object of class \code{"glmmadmb"} representing the model fit,
  including (among others) components:
  \item{b}{vector of fixed effects}
  \item{S}{covariance matrix of random effects}
  \item{alpha}{scale/overdispersion parameter
    (negative binomial, Gamma, beta)}
  \item{pz}{Zero-inflation parameter (only when
    \code{zeroInflation=TRUE})}
  \item{phi}{Matrix for converting from 'orthogonalized' to 'real'
    parameters: see Details}
}
\note{
  \itemize{
    \item \code{glmm.admb} is defined for backward compatibility; it will be removed in a future version.
    \item Tools for working with MCMC output are preliminary; see \code{\link{mcmcControl}} or the vignette for more information.
  }
}
\references{
  Bohning, D., E. Dietz, and P. Schlattmann. 1999. The zero-inflated
  Poisson model and the decayed, missing and filled teeth index in
  dental epidemiology. Journal of the Royal Statistical Society. Series
  A (Statistics in Society) 162:195--209.

  Skaug, H.J. and D.A. Fournier. 2006. Automatic approximation of the
  marginal likelihood in non-Gaussian hierarchical models. Computational
  Statistics & Data Analysis 51:699--709.

  Morris, W. 1997. American Naturalist 150:299-327
}
\seealso{
  \code{\link{glmmADMB-package}}.
}
\examples{
## Australian H. influenzae data from MASS
data(bacteria,package="MASS")
bacteria$present <- as.numeric(bacteria$y)-1
if (!check_rforge())
    (bfit <-  glmmadmb(present ~ trt + I(week > 2), random = ~ 1 | ID,
                       family = "binomial", data = bacteria))
## simple simulated zero-inflated Poisson example
### simulate values
set.seed(101)
d <- data.frame(f=factor(rep(LETTERS[1:10],each=10)),x=runif(100))
u <- rnorm(10,sd=2)
d$eta <- with(d,u[f]+1+4*x)
pz <- 0.3
zi <- rbinom(100,size=1,prob=pz)
d$y <- ifelse(zi,0,rpois(100,lambda=exp(d$eta)))
## fit
\dontrun{
    zipmodel <- glmmadmb(y~x+(1|f),data=d,family="poisson",zeroInflation=TRUE)
}

\dontrun{
## Epilepsy data
## all grouping variables must be factors:
  epil2$subject <- factor(epil2$subject)
  (fm <- glmmadmb(y~Base*trt+Age+Visit+(Visit|subject),
                  data=epil2, family="nbinom"))
## Owls data
 om <- glmmadmb(SiblingNegotiation~FoodTreatment*SexParent+
           (1|Nest)+offset(log(BroodSize)),
          zeroInflation=TRUE,family="nbinom",data=Owls)
}
}
\keyword{models}
