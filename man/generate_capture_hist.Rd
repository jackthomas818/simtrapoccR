% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_capture_hist.R
\name{generate_capture_hist}
\alias{generate_capture_hist}
\title{Generate Capture History Data}
\usage{
generate_capture_hist(ntraps, nsample_cap, tau, p, homeranges)
}
\arguments{
\item{ntraps}{number of camera traps to place in study area}

\item{nsample_cap}{number of capture-recapture sampling occasions}

\item{tau}{movement parameter of species}

\item{p}{the species capture probability}

\item{homeranges}{array of homerange centers for all individuals}
}
\value{
capture_hist a 2D array (n rows by nsample_cap columns) of
the capture history for each individual, where 1 indicates a
capture.
}
\description{
Function that generates the capture history array. Uses camera trapping as
the capture-recapture method. NOTE: the study
area is assumed to be a unit square (1x1).
}
\examples{
# number of individuals in capture history
n <- 10

#number of sampling occasions in capture-recapture
nsample_cap <- 5

#amount of movement for the species
tau <- 0.25

# number of camera traps
ntraps <- 3

# species capture probability
p <- 0.6

# create homeranges for the n individuals
homeranges <- create_homeranges(n)

capture_hist <- generate_capture_hist(ntraps, nsample_cap, tau, p, homeranges)
}
