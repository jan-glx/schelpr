% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_it.R
\name{bin_it}
\alias{bin_it}
\title{Bin a vector into a predefined number of as-equally-sized-as-possible intervals}
\usage{
bin_it(x, nbins = 5L)
}
\arguments{
\item{x}{Numeric vector to be binned.}

\item{nbins}{Integer indicating the number of bins to create. Defaults to 5.
Must be a positive integer with a length of 1.}
}
\value{
An ordered factor with levels corresponding to the bins.
}
\description{
This function bins a vector of comparable values into a specified number of as-equally-sized intervals.
It uses the empirical cumulative distribution function (ECDF) and an optimal (in terms of entropy) partioning approach
to determine the bin edges, ensuring a more informative binning process compared to simple range division.
It can be considered a better version of base:cut
}
\details{
The function first checks for input validity, ensuring that `x` is non-empty
and `nbins` is a valid integer greater than or equal to 1. It then computes the ECDF
of the input vector and applies an optimal binning algorithm (using find_best_partioning)
to find the best bin edges based on the desired number of bins. If `nbins` is set to 1,
the function returns an empty matrix.
}
\examples{
set.seed(123)
data <- rpois(100, rexp(100))
binned_data <- bin_it(data, nbins = 5)
binned_data
table(binned_data)

}
