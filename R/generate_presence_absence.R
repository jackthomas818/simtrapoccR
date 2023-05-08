#' Setup Grid Area
#' function that divides the study area (a unit square) into a grid of equally
#' divided regions.
#' @param nsites number of grid cells required
#'
#' @return regions, a 2-d array with grid cells defined as x0,x1,y0,y1
#'                  bounding boxes
#'
#' @export
#'
#' @examples
#'
#' # number of individuals
#' n <- 5
#' homeranges <- create_homeranges(n)
#' nsites <- 16
#' regions <- setup_grid(nsites)
#' #
#' # plot grid cells on study area with homerange centers and individuals labeled
#' par(pty = "s")
#' plot(homeranges, xlim = c(0, 1), ylim = c(0, 1))
#' m <- sqrt(nsites)
#' for (i in 0:m) {
#'   abline(v = i / m)
#'   abline(h = i / m)
#' }
#' # individuals labeled from 1 to n
#' text(homeranges, labels = rownames(homeranges), cex = 1.5, font = 2, pos = 3)
setup_grid <- function(nsites) {
  # We now have a m x m grid for the unit square as our divided area

  # number of grid cells per row
  m <- sqrt(nsites)

  # x coordinates
  regions_x <- array(data = NA, dim = c(m * m, 2))

  for (k in 0:m - 1) {
    for (i in 1:m) {
      regions_x[i + k * m, 1] <- (i - 1) / m
      regions_x[i + k * m, 2] <- regions_x[i + k * m, 1] + 1 / m
    }
  }

  # y coordinates
  regions_y <- array(data = NA, dim = c(m * m, 2))
  j <- 0
  for (i in 1:(m * m)) {
    regions_y[i, 1] <- j
    regions_y[i, 2] <- j + 1 / m

    if (i %% m == 0) {
      j <- j + 1 / m
    }
  }

  # regions are in x0,x1,y0,y1 form
  regions <- cbind(regions_x, regions_y)

  return(regions)
}

#' Calculate Proportion of Time Spent in Grid Cell
#'
#' @param x The x coordinate of the homerange for an individual
#' @param y The y coordinate of the homerange for an individual
#' @param x0 The starting x value for the grid cell
#' @param x1 The ending x value for the grid cell
#' @param y0 The starting y value for the grid cell
#' @param y1 The ending y value for the grid cell
#' @param tau The movement parameter for the species, lower is less movement
#'
#' @return A value between 0 and 1 indicating the proportion of time that the
#'         individual spends in the given grid cell
#' @export
#'
#' @import OpenMx
#' @examples
#' # show high proportion when homerange in center of grid cell
#' calc_prop(0.25 / 2, 0.25 / 2, 0, 0.25, 0, 0.25, 0.1)
calc_prop <- function(x, y, x0, x1, y0, y1, tau) {
  covariance <- matrix(c(tau^2, 0, 0, tau^2), nrow = 2, ncol = 2)
  means <- c(x, y)
  lbound <- c(x0, y0)
  ubound <- c(x1, y1)

  prop <- OpenMx::omxMnor(covariance, means, lbound, ubound)[1]
  return(prop)
}


#' Simulate Presence-Absence Data
#'
#' Function that simulates presence-absence data for each region and for
#' each period.
#'
#' Assumes that the amount of time spent by each individual
#' in each region is proportional to the probability of detection.
#'
#' @param nsites number of regions in the study area
#' @param nsample_pa number of sampling occasions
#' @param homeranges the homeranges of the individuals as 2D array of coordinates
#' @param regions the subareas of the study area, as 2D array of (x0,x1,y0,y1)
#'                bounding boxes
#' @param tau the movement parameter for the species, higher is more movement
#' @param sigma the marker deposition rate (integer)
#' @param delta the time between sampling occasions (integer)
#'
#' @return presence_absence a 2D array (nsites rows by nsample_pa columns)
#'         of 0s and 1s where 1 indicates a detection for a site and a period.
#'
sim_presence_absence <- function(nsites, nsample_pa, homeranges, regions, tau, sigma, delta) {
  # simulate presence-absence data for each region for each period nsample_occ

  m <- sqrt(nsites)
  # 2D matrix with rows being regions and columns being time periods
  presence_absence <- array(data = 0, dim = c(m * m, nsample_pa))

  indi_props <- array(data = NA, dim = c(nrow(regions), nrow(homeranges)))
  for (i in 1:nrow(regions)) {
    # compute proportion of time each individual spends in this region
    for (j in 1:nrow(homeranges)) {
      x <- homeranges[j, 1]
      y <- homeranges[j, 2]

      x0 <- regions[i, 1]
      x1 <- regions[i, 2]

      y0 <- regions[i, 3]
      y1 <- regions[i, 4]

      prop <- calc_prop(x, y, x0, x1, y0, y1, tau)
      # P(presence detection) = P(Z_{ij}>0)
      indi_props[i, j] <- 1 - exp(-sigma * delta * prop)
    }

    # run simulation for presence-absence data
    # in this region
    for (j in 1:nrow(homeranges)) {
      for (k in 1:nsample_pa) {
        # assume detection probability is proportional to time spent in region
        detection_prob <- indi_props[i, j]
        rv <- stats::runif(1)

        if (detection_prob >= rv) {
          presence_absence[i, k] <- 1
        }
      }
    }
  }
  # columns do not necessarily sum to 1, because individual may spend
  # time outside the study area
  # print(indi_props)
  return(presence_absence)
}

#' Generate Presence-Absence Data
#'
#' Function that generates the presence-absence array. NOTE: the study
#' area is assumed to be a unit square (1x1).
#'
#' @param nsites number of sites in study area
#' @param nsample_pa number of presence-absence sampling occasions
#' @param homeranges array of homerange centers for all individuals
#' @param tau movement parameter of species
#' @param sigma the marker deposition rate
#' @param delta the time between sampling occasions
#'
#' @return presence_absence a 2D array (nsites rows by nsample_pa columns)
#'         of 0s and 1s where 1 indicates a detection for a site on a period.
#' @export
#'
#' @examples
#'
#' # number of sites
#' nsites <- 4**2
#'
#' # number of individuals in study area
#' n <- 10
#'
#' # amount of movement for the species
#' tau <- 0.25
#'
#' # marker deposition rate
#' sigma <- 1
#'
#' # time between sampling occasions (presence-absence)
#' delta <- 2
#'
#' # number of sampling occasions in presence-absence
#' nsample_pa <- 7
#'
#' # homeranges for the n individuals
#' homeranges <- create_homeranges(n)
#'
#' presence_absence <- generate_presence_absence(
#'   nsites, nsample_pa, tau, sigma,
#'   delta, homeranges
#' )
generate_presence_absence <- function(nsites, nsample_pa, tau, sigma, delta, homeranges) {
  # setup the study area grid
  regions <- setup_grid(nsites)

  # simulate the presence-absence data
  presence_absence <- sim_presence_absence(nsites, nsample_pa, homeranges, regions, tau, sigma, delta)

  return(presence_absence)
}
