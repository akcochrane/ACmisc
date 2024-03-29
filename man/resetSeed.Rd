% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resetSeed.R
\name{resetSeed}
\alias{resetSeed}
\title{Re-randomize seed}
\usage{
resetSeed()
}
\description{
If the seed has been set, this will un-set it.
}
\examples{
set.seed(13) # set an arbitrary seed
rnorm(5)     # see the resulting "random" values
set.seed(13) # re-set the same arbitrary seed
rnorm(5)     # see the same "random" values
set.seed(13) # re-set the same arbitrary seed
resetSeed()  # un-set the current seed
rnorm(5)     # see different "random" values

}
