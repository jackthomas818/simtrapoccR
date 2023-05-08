# simtrapoccR

<!-- badges: start -->
<!-- badges: end -->

Simulate Camera Trap Capture-Recapture and Presence-Absence Occupancy Data

The 'simtrapoccR' package allows you to generate capture-recapture history matrices and
presence-absence history matrices for feeding into population estimation statistical models.

For information about how the data is generated, see [here](https://github.com/jackthomas818/masters-project).

## Installation

You can install the development version of simtrapoccR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jackthomas818/simtrapoccR")
```

## Usage

``` r
# number of sites
nsites <- 4**2

# number of individuals in capture history
n <- 10

# number of sampling occasions in capture-recapture
nsample_cap <- 5

# number of sampling occasions in presence-absence
nsample_pa <- 7

# amount of movement for the species
tau <- 0.25

# number of camera traps
ntraps <- 3

# species capture probability
p <- 0.6

# marker deposition rate
sigma <- 1

# time between sampling occasions (presence-absence)
delta <- 2

# homeranges for the n individuals
homeranges <- create_homeranges(n)

capture_hist <- generate_capture_hist(ntraps, nsample_cap, tau, p, homeranges)

presence_absence <- generate_presence_absence(nsites, nsample_pa, tau, sigma, delta, homeranges)
```
