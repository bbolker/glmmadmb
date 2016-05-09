\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{glmmADMB News}
\encoding{UTF-8}

\section{Changes in version 0.8.2 (18-08-2015)}{
  \itemize{
    \item update binaries (at least for Windows, MacOS) to version
    r2862
  }
}
\section{Changes in version 0.8.0 (23-05-2014)}{
  \itemize{
    \item changed suggested package 'scapeMCMC' to 'plotMCMC'
  }
}
\section{Changes in version 0.7.6 (21-05-2013)}{
  \itemize{
    \item added Poisson-zeroprob hack to TPL (only in Linux binary so far)
    \item stubbed out possibility of corStruct>1; added comments
  }
}

\section{Changes in version 0.7.5 (7-05-2013)}{
  \itemize{
    \item vignette tweaks
    \item fixed corStruct="full"
  }
}

\section{Changes in version 0.7.4 (20-02-2013)}{
 \itemize{
\item updated glmmADMB binaries
\item beta-binomial models (NOT YET WORKING)
\item stubs for model matrices for alpha, zi effects
\item output tweaks
}
}

\section{Changes in version 0.7.4 (20-02-2013)}{
 \itemize{
\item updated glmmADMB binaries
\item beta-binomial models (NOT YET WORKING)
\item stubs for model matrices for alpha, zi effects
\item output tweaks
}
}

\section{Changes in version 0.7.3 (17-11-2012)}{
 \itemize{
\item fixed bug in fitted/residuals calculation (failed to scale random effects)
}
}
\section{Changes in version 0.7.2.12}{
\itemize{
  \item added citation info
  \item MCMC results transformed automatically (but *not* correlation
  parameters or random effects; the function warns of this)
  \item fixed buglet with empty fixed formula
  \item drop1() works
  \item fixed tiny glitch in glmmADMB model printing (alpha printed even
  if not relevant)
  \item now uses results from .bar file for greater precision in fixed
  effects (mostly useful for cross-package comparisons)
  \item step, stepAIC **DO NOT WORK**; added masking functions
  \item fixed MacOS X.5 detection code (throws warning)
}
}

\section{Changes in version 0.7.2.9}{
\itemize{
  \item  version bump, new version of saved OwlModels
  \item added df.residuals and model.frame methods, car::Anova now works
  \item added car, lme4.0 to Suggests:
  \item vignette tweaks/additions
}
}

\section{Changes in version 0.7.2.8}{
\itemize{
  \item  added information to vignette about parameter conversions from MCMC to original scales
}
}

\section{Changes in version 0.7.2.7}{
\itemize{
  \item anova() is a little more flexible, tiny format changes
  \item removed dat_write (use version from R2admb)
}
}

\section{Changes in version 0.7.2.6}{
\itemize{
  \item fix scaling of RE standard deviations; add scale option to ranef
}
}

\section{Changes in version 0.7.2.5}{
\itemize{
  \item new binaries (ADMB r229)
  \item R calculation of phi, pre-transformation of fixed effect parameters
}
}

\section{Changes in version 0.7.2.4}{
 \itemize{
   \item Bug fixes in sys.result handling
   \item Retrieve and document information about orthogonalized
   fixed-effect parameters
   \item re-instituted informative error message for non-factor grouping
   variables
   \item added maxphase (-maxph) to admbControl
   \item further robustification of TPL file
   \item documentation cleanup/improvement
 }
}
\section{Changes in version 0.7.2 (2012-01-06)}{
  \itemize{
    \item updated TPL file; improved robustness, reduced warnings
    \item added get_bin_version
    \item improvements to reporting (drop unused levels in nested
  groups, better labeling of variance components and random effects, return 'alpha'
  (variance/scale parameter etc.) for all appropriate distributions);
  better formatting of RE variance-covariance matrices
  \item documented default link functions
  }
}

\section{Changes in version 0.7.1 (2011-12-28)}{
  \itemize{
    \item added logistic distribution, Gaussian (fixed stupid bug)!
    \item fixed bugs in predict function with offsets, documented predict.glmmadmb
  }
}

\section{Changes in version 0.7 (2011-12-18)}{
  \itemize{
    \item added R2admb dependency, removed copies of R2admb files
    \item improved vignette (hurdle and MCMC)
    \item improved labeling of MCMC output
    \item added test for MacOS <10.6
    \item fixed predict bug; enhanced predict (still needs documentation)
  }
}
  
\section{Changes in version 0.6.7.1 (2011-12-05)}{
  \itemize{
    \item SVN rev 158
    \item added debug flag, debugging statements
    \item attempt to avoid temporary-directory collisions
    \item TPL file adds truncated Poisson and negative binomial
    \item update saved OwlModel
    \item added 'run' option to admbControl
    \item fixed predict bug (remove LHS)
    \item fix problems with calculating sd.est for models other than
  Bernoulli, Poisson, NB
  }
}

\section{Changes in version 0.6.6 (2011-11-23)}{
  \itemize{
    \item SVN rev 147
    \item add identity, cloglog links
    \item stub for gaussian family (not yet working!)
    \item use tolower() on family
    \item add (experimental) predict() function
    \item comment on easyFlag
    \item move ZI_kluge to admbControl
    \item fix name assignment bug for random effect levels
    
  }
}
  
\section{Changes in version 0.6.5 (2011-10-22)}{
  \itemize{
    \item SVN rev 142
    \item BRANCH MERGE with alpha version
    \item added mcmcControl, admbControl
    \item added shess, noinit options (poorly documented)
    \item added scapeMCMC, ggplot2 to Suggests: (the latter for vignette use)
  }
}

\section{Changes in version 0.6.4 (2011-08-30)}{
  \itemize{
    \item new binaries
  }
}

\section{Changes in version 0.6.3 (2011-08-08)}{
  \itemize{
    \item include more buildbot binaries
    \item added nbinom1 option
    \item fixed phase bugs in nbinom
  }
}

\section{Changes in version 0.6.2}{
  \itemize{
    \item include sd of zero-inflation in output and summary
    \item include buildbot binaries
    \item more information when failing to find binaries
  }
}

\section{Changes in version 0.6.1 (2011-07-08)}{
  \itemize{
    \item Allow random effects to be specified separately in 'random' (i.e. as in 'old' glmmADMB)
    \item Add fixef() as synonym for coef(), for lme4-compatibility
    \item Add 'start' argument (poorly tested)
}}

\section{Changes in version 0.6 (2011-05-04)}{
  \itemize{
    \item Changed formula to lme4-style: random effects, and offsets,
    included within the formula
    \item Allow multiple random effects
    \item Turn off minimal zero-inflation in non-ZI models (pz was
    always 0.0001)
    \item Require grouping variables to be factors
    \item New families: Binomial with N>1, Gamma, stub for Beta
    \item Allow alternative link functions (poorly tested): logit,
    probit, log
    \item new accessors: VarCorr(); ranef() [returns scaled random
    effects]; updated summary method (includes Wald tests); vcov();
    stdEr()
    \item Basic access to MCMC functionality
    \item function and class name changed from glmm.admb to glmmadmb
}}
\section{CHANGES IN VERSION 0.5-2 (2011-02-09)}{
  \itemize{
    \item SVN rev 54
    \item Added ranef(), stdEr(), and vcov() methods.
    \item Added 'type' argument to residuals().
    \item Improved summary() method.
    \item Added 'append' argument to dat_write().
    \item Added element 'corMat' to glmm.admb() output.
    \item Added attributes 'df' and 'nobs' to logLik() output.
    \item Added AIC to print() output.
    \item Added column 'logBroodSize' to Owls data.
    \item Added OwlModel example object.
    }}

\section{CHANGES IN VERSION 0.5-1 (2011-01-18)}{
  \itemize{
    \item SVN rev 32
    \item Fixed scrambled order of fitted values and predictions.
    }}

\section{CHANGES IN VERSION 0.5-0 (2011-01-17)}{
\itemize{
  \item SVN rev 30
  \item Added coef(), residuals(), and fitted() methods.
  \item Added bvprobit.tpl and nbmm.tpl ADMB source files.
  \item Added Owls example data.
  }}

\section{CHANGES IN VERSION 0.4-2 (2010-10-26)}{
  \itemize{
    \item SVN rev 23 [2010-11-05]
    \item Added MacOS executables.
    \item Added 'verbose' argument to glmm.admb().
    }}

\section{Changes in version 0.4-1 (2010-10-26)}{
  \itemize{
    \item Fixed executable bits in Linux.
    \item Improved backward compatibility with old versions of R.
}}
\section{Changes in version 0.4-0 (2010-10-17)}{
  \itemize{
    \item Removed unnecessary internal functions.
    \item Removed unnecessary warnings.
}}
\section{Changes in version 0.3-0 (2010-07-28)}{
  \itemize{
    \item Legacy version, now on R-Forge.
}}
