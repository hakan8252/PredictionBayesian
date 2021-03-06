% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Information.R
\name{entropy}
\alias{entropy}
\title{Entropy of Variable}
\usage{
entropy(x)
}
\arguments{
\item{x}{a numerical vector}
}
\value{
Total entropy of numerical vector.
}
\description{
Calculate total entropy of each variable in a numerical vector. However, if a vector contains value of `0`.
\code{correction()} function is applied.
}
\examples{
vector.1 <- c(1,4,5,0)
entropy(vector.1)
}
\references{
B. McMillan and D. Slepian, “Information Theory,” Proc. IRE, vol. 50, no. 5, pp. 1151–1157, 1962.
}
