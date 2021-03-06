% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiple_testing.R
\name{adjust_q_vals}
\alias{adjust_q_vals}
\alias{adjust_q_vals.default}
\alias{adjust_q_vals.list}
\title{Combine sets of FDR controlled tests and adjust for this second level of multiple testing
the returned adjusted q values control the average FDR for the sets with rejections (sets with no rejections always control FDR)}
\usage{
adjust_q_vals(q_vals, ...)

\method{adjust_q_vals}{default}(q_vals, sets)

\method{adjust_q_vals}{list}(q_vals, ...)
}
\description{
Combine sets of FDR controlled tests and adjust for this second level of multiple testing
the returned adjusted q values control the average FDR for the sets with rejections (sets with no rejections always control FDR)
}
