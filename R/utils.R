#' Distance function
#'
#' @param a point a in vector form
#' @param b point b in vector form
#'
#' @return euclidean distance between point a and point b
#' @export
#'
#' @examples
#' a <- c(2, 6)
#' b <- c(3, 5)
#' distance(a, b)
distance <- function(a, b) {
  return(sqrt(sum((a - b)**2)))
}


#' Half-Normal Detection Function
#'
#' A half-normal detection function models the capture probability of individual i by trap j:
#' <p>p<sub>ij</sub> = p e<sup>(-d<sup>2</sup><sub>ij</sub>/&tau;<sup>2</sup>)</sup></p>
#' <p>where <em>p</em> is the capture probability for a camera trap at distance 0 from the home range
#' centre, <em>&tau;</em> is the movement parameter for the species, and <em>d<sub>ij</sub></em>
#' is the distance between the home-range centre <em>i</em> and the camera location <em>j</em>.</p>

#'
#' @param d distance of homerange center to trap
#' @param tau movement parameter for species
#'
#' @return probability of detection for individual by camera trap
#' @export
#'
#' @examples
#' # example on species with tau=0.1
#' trap <- c(0.25, 0.5)
#' homerange <- c(0.26, 0.6)
#' tau <- 0.1
#'
#' d <- distance(trap, homerange)
#' detect_prob <- half_norm_detec(d, tau)
half_norm_detec <- function(d, tau) {
  return(exp((-d**2) / (tau**2)))
}