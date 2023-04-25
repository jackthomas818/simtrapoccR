#' Create Homeranges
#'
#' Function that assigns homeranges to n individuals using a Poisson
#' random process
#'
#' Assumes a study area that is a unit square (0,1) x (0,1)
#'
#' @param n number of individuals
#'
#' @return a 2D array of the homeranges for each individual
#' @export
#'
#' @examples
#' n <- 10
#' homeranges <- create_homeranges(n)
#' # plotting homeranges for validation
#' par(pty = "s")
#' plot(homeranges, xlim = c(0, 1), ylim = c(0, 1))
#' text(homeranges, labels = rownames(homeranges), cex = 1.5, font = 2, pos = 3)
create_homeranges <- function(n) {


  # our study area is the unit square
  # define homerange centers

  # Poisson process, with x, y coordinates as independent uniform
  # random variables

  homeranges <- array(data = NA, dim = c(n, 2), dimnames = list(c(1:n), c("x", "y")))

  for (i in 1:n) {
    x <- stats::runif(1)
    y <- stats::runif(1)
    homeranges[i, 1] <- x
    homeranges[i, 2] <- y
  }

  return(homeranges)
}
