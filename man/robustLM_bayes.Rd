% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustLM_bayes.R
\name{robustLM_bayes}
\alias{robustLM_bayes}
\title{Fit a robust model, with fully Bayesian model comparisons}
\usage{
robustLM_bayes(formula, data, ..., iter = 10000, showProgress = T, quiet = T)
}
\arguments{
\item{formula}{Formula, as in \code{lm()}, \code{brm()} or \code{lmer()}}

\item{data}{Data, as in \code{lm()} or \code{lmer()}}

\item{...}{Additional arguments to pass to \code{\link[brms]{brm}}. For example, if multiple cores are available, \code{cores=2} or more should drastically speed up performance.}

\item{iter}{Number of iterations. More is better, within reasonable time constraints}

\item{showProgress}{Show a small progress bar?}

\item{quiet}{Logical. To show model progress, or hide it.}
}
\description{
Fits a median linear regression, with random effects (i.e., mixed-effects model) if desired,
using \code{\link[brms]{brm}} and the asymmetric laplace distribution (fitting quantile = .5, i.e., median).
Then fits a set of models dropping each fixed-effect term in turn, and uses
bridge sampling to estimate the Bayes Factors of the full model relative to each
sub-model. Returns the full model, with a data frame \code{model$bayes_factors} with
each fixed effect's log (base 3) Bayes Factor.
}
\details{
Bridge sampling is run 5 times for each effect, and the most conservative (i.e., equivocal)
of the resulting BF is returned.

These models are computationally expensive and require quite a bit of memory. Without enough
memory they will likely crash the R session. Other sources of errors or crashing may be 
out-of-date packages, such as \code{RcppParallel}, \code{StanHeaders}, \code{brms}, \code{bridgesampling}.

While using \code{algorighm = 'fullrank'} speeds up the model (over the default \code{algorithm  = 'sampling'} by a factor of over 3.5 in preliminary tests), it would be advised
to fit and compare several fitting runs to assess the robustness of the resulting coefficients and Bayes Factors. In brief tests, some
values have been quite consistent across fitting runs and others have not when using \code{fullrank} or \code{meanfield}; in general, waiting for
\code{sampling} provides the most robust estimates.
}
\examples{
\dontrun{
# This fairly simple example takes about 160 seconds on a fairly good computer.
m <- robustLM_bayes(mot ~ ufov * enum, dat_cochraneEtAl_2019_PLOSOne) 

m$bayes_factors
}
}
\seealso{
\code{\link[brms]{brm}}, \code{\link[brms]{bayes_factor}}, and \code{\link[brms]{brmsfamily}}.
}
