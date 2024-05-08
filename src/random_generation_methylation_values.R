
#' Random generation for beta distribution simulating methylation data
#' 
#' @param n An \code{integer}: number of observations.
#' @return A \code{numeric} vector of length n following a beta distribution
#'         simulating the distribution of methylation data.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Random generation of 1000 methylation values following a bimodal distribution
#' sim_values <- rbimodal_meth(n = 1000)
#' hist(sim_values, breaks = 20)
#' @references
#' \href{https://stats.stackexchange.com/questions/355344/simulating-a-bimodal-distribution-in-the-range-of-15-in-r}{Simulating a bimodal distribution}

rbimodal_meth <- function(n){
  bimodal_distrib <- rbeta(n = n, shape1 = 0.2, shape2 = 0.2, ncp = 0)
  return(bimodal_distrib)
}

#' Random generation for trimodal distribution simulating bulk methylation data
#' 
#' @param n An \code{integer}: number of observations.
#' @return A \code{numeric} vector of length n following a trimodal distribution
#'         simulating the distribution of bulk methylation data.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' # Random generation of 1000 methylation values following a trimodal distribution
#' sim_values <- rtrimodal_meth(n = 1000)
#' hist(sim_values, breaks = 20)
#' @references
#' \href{https://CRAN.R-project.org/package=truncnorm}{Mersmann O, Trautmann H, Steuer D, Bornkamp B (2023). _truncnorm: Truncated Normal Distribution_. R package version 1.0-9}

rtrimodal_meth <- function(n){
  low_distrib <- truncnorm::rtruncnorm(n = n, a = 0, b = 1, mean = 0, sd = 0.05)
  high_distrib <- truncnorm::rtruncnorm(n = n, a = 0, b = 1, mean = 1, sd = 0.05)
  mid_distrib <- truncnorm::rtruncnorm(n = n, a = 0, b = 1, mean = 0.5, sd = 0.05)
  trimodal_distrib <- sample(
    x = c(low_distrib, mid_distrib, high_distrib), size = n)
  return(trimodal_distrib)
}
