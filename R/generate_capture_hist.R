#' Define Trap Locations
#'
#' Function that places camera traps in the study area (a unit square).
#'
#' Uses a Poisson process where x, y coordinates of each trap
#' are independent uniform random variables between 0 and 1.
#'
#' @param ntraps the number of camera traps to place
#'
#' @return traps a 2D array of trap locations in x,y coordinates
#' @export
#'
#' @examples
#' n <- 5
#' homeranges <- create_homeranges(n)
#' ntraps <- 3
#' traps <- define_traps(ntraps)
#' # plot traps on study area with homeranges included
#' par(pty = "s", mar = c(5, 4, 4, 8), xpd = TRUE)
#'
#' plot(homeranges, xlim = c(0, 1), ylim = c(0, 1))
#' points(traps, col = "red", pch = 15)
#' legend("topright",
#'   inset = c(-.5, 0), c("homerange", "camera trap"),
#'   cex = .8, col = c("black", "red"), pch = c(1, 15)
#' )
define_traps <- function(ntraps) {
  traps <- array(data = NA, dim = c(ntraps, 2), dimnames = list(c(1:ntraps), c("x", "y")))

  # randomly place camera traps
  for (i in 1:ntraps) {
    x <- stats::runif(1)
    y <- stats::runif(1)
    traps[i, 1] <- x
    traps[i, 2] <- y
  }

  return(traps)
}

#' Calculate Trap Detection Probabilities
#'
#' Function that calculates the detection probabilities for each
#' individual i and each trap j
#'
#' @param ntraps number of traps
#' @param homeranges homeranges for all individuals
#' @param traps locations of all traps
#' @param tau movement parameter for species
#' @param p capture probability for species
#'
#' @return trap_detec 2D array of detection probabilities for each individual i
#'                    and each trap j.
#' @export
#'
#' @examples
#'
#' # plot individuals and their camera trap detection probabilities
#' # labelled for trap j. Blue label is individual i number.
#'
#' # number of individuals in population
#' n <- 5
#' homeranges <- create_homeranges(n)
#'
#' # number of camera traps
#' ntraps <- 3
#'
#' # movement parameter
#' tau <- 0.55
#'
#' # species capture probability
#' p <- 0.6
#'
#' traps <- define_traps(ntraps)
#' trap_detec <- calc_detection(ntraps, homeranges, traps, tau, p)
#'
#' for (j in 1:ntraps) {
#'   par(pty = "s", mar = c(5, 4, 4, 8), xpd = TRUE)
#'   plot(homeranges, xlim = c(0, 1), ylim = c(0, 1))
#'   points(x = traps[j, 1], y = traps[j, 2], col = "red", pch = 15)
#'   text(
#'     x = traps[j, 1], y = traps[j, 2], labels = paste(c("Trap", j), collapse = " "),
#'     cex = 1.5, font = 2, pos = 3
#'   )
#'   text(homeranges, labels = round(trap_detec[, j], 2), cex = 1, font = 2, pos = 3)
#'   text(homeranges, labels = rownames(homeranges), cex = 1, font = 2, pos = 2, col = "blue")
#'   legend("topright",
#'     inset = c(-.5, 0), c("homerange", "camera trap"), cex = .8,
#'     col = c("black", "red"), pch = c(1, 15)
#'   )
#' }
calc_detection <- function(ntraps, homeranges, traps, tau, p) {
  n <- nrow(homeranges)
  # 2D matrix of detection probabilities for each
  # individual i and each trap j
  trap_detec <- array(data = 0, dim = c(n, ntraps))

  # calculate camera trap detection probabilities

  for (i in 1:nrow(homeranges)) {
    for (j in 1:ntraps) {
      # calculate distance of homerange to camera trap j
      d <- distance(homeranges[i, ], traps[j, ])
      # calculate probability of detection
      detec_prob <- half_norm_detec(d, tau)

      trap_detec[i, j] <- p * detec_prob
    }
  }
  return(trap_detec)
}


#' Simulate Trap Interactions
#'
#' Function that simulates interactions of individuals with the camera traps.
#'
#' @param homeranges homeranges of the individuals
#' @param ntraps number of camera traps
#' @param nsample_cap number of capture-recapture sampling occasions
#' @param trap_detec  2D array of detection probabilities for each individual i
#'                    and each trap j.
#'
#' @return camera_trap_full a 3D array of whether each individual i was seen by
#'                          trap j at period k
#' @export
#' @examples
#'
#' # movement parameter
#' tau <- 0.55
#'
#' # number of individuals in population
#' n <- 5
#'
#' # species capture probability
#' p <- 0.6
#'
#' homeranges <- create_homeranges(n)
#' ntraps <- 3
#'
#' traps <- define_traps(ntraps)
#' trap_detec <- calc_detection(ntraps, homeranges, traps, tau, p)
#'
#' # number of capture-recapture sampling occasions
#' nsample_cap <- 5
#'
#' camera_trap_full <- sim_trap_interac(homeranges, ntraps, nsample_cap, trap_detec)
#'
#' # plot number of times each individual was seen by each trap
#' trap_counts <- rowSums(camera_trap_full, dims = 2)
#'
#' for (j in 1:ntraps) {
#'   barplot(trap_counts[, j],
#'     names.arg = rownames(homeranges),
#'     main = paste(c("Trap", j), collapse = " "), xlab = "individual", ylab = "count"
#'   )
#' }
sim_trap_interac <- function(homeranges, ntraps, nsample_cap, trap_detec) {
  # simulate interactions of individuals with the traps
  camera_trap_full <- array(data = 0, dim = c(nrow(homeranges), ntraps, nsample_cap))

  # for each individual i
  for (i in 1:nrow(homeranges)) {
    # for each trap j
    for (j in 1:ntraps) {
      # for each period k
      for (k in 1:nsample_cap) {
        # if individual detected by trap j, put 1 in matrix at i,j,k
        rv <- stats::runif(1)
        if (trap_detec[i, j] >= rv) {
          camera_trap_full[i, j, k] <- 1
        }
      }
    }
  }
  return(camera_trap_full)
}


#' Create Capture History
#'
#' Function that collapses the 3D array camera_trap_full (individual i, trap j, period k)
#' on the jth dimension to obtain the capture history for each individual i.
#'
#' @param homeranges homeranges for each individual
#' @param nsample_cap number of capture-recapture sampling occasions
#' @param camera_trap_full a 3D array of whether each individual i was seen by
#'                         trap j at period k
#'
#' @return capture_hist a 2D array of capture histories for each individual
#'
create_cap_hist <- function(homeranges, nsample_cap, camera_trap_full) {
  capture_hist <- array(data = 0, dim = c(nrow(homeranges), nsample_cap))

  # for each individual
  for (i in 1:nrow(homeranges)) {
    # for each period k
    for (k in 1:nsample_cap) {
      # if individual was detected by any trap on period k, add 1 to capture history
      if (sum(camera_trap_full[i, , k]) >= 1) {
        capture_hist[i, k] <- 1
      }
    }
  }
  return(capture_hist)
}

#' Generate Capture History Data
#'
#' Function that generates the capture history array. Uses camera trapping as
#' the capture-recapture method.
#'
#' @param ntraps number of camera traps to place in study area
#' @param nsample_cap number of capture-recapture sampling occasions
#' @param tau movement parameter of species
#' @param homeranges array of homerange centers for all individuals
#' @param p the species capture probability
#'
#' @return capture_hist a 2D array (n rows by nsample_cap columns) of
#'         the capture history for each individual, where 1 indicates a
#'         capture.
#' @export
#'
#' @examples
#' # number of individuals in capture history
#' n <- 10
#'
#' #number of sampling occasions in capture-recapture
#' nsample_cap <- 5
#'
#' #amount of movement for the species
#' tau <- 0.25
#'
#' # number of camera traps
#' ntraps <- 3
#'
#' # species capture probability
#' p <- 0.6
#'
#' # create homeranges for the n individuals
#' homeranges <- create_homeranges(n)
#'
#' capture_hist <- generate_capture_hist(ntraps, nsample_cap, tau, p, homeranges)
generate_capture_hist <- function(ntraps, nsample_cap, tau, p, homeranges) {
  # place camera traps in study area
  traps <- define_traps(ntraps)
  # calculate the probabilities of detection for each individual by each trap
  trap_detec <- calc_detection(ntraps, homeranges, traps, tau, p)
  # generate the 3D array of which individual (i) was seen by which trap (j) at each period (k)
  camera_trap_full <- sim_trap_interac(homeranges, ntraps, nsample_cap, trap_detec)
  # collapse the jth dimension to obtain the capture history for each individual
  capture_hist <- create_cap_hist(homeranges, nsample_cap, camera_trap_full)

  return(capture_hist)
}
