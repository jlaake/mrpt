\name{delta}
\alias{delta}
\title{Delta function for Mark-Recapture Distance Sampling}
\description{
Computes delta value at given x as described in Buckland et al. (2010)
}
\usage{
delta(x, p1, p2, gamma, PI , use.offset, posdep)
}
\arguments{
  \item{x}{Design matrix for delta function}
  \item{p1}{Value of detection probabilities for observer 1}
  \item{p2}{Value of detection probabilities for observer 2}
  \item{gamma}{Parameters for delta function}
  \item{PI}{If TRUE, uses point independence assumption}
  \item{use.offset}{If TRUE, uses offset in delta function so PI and FI models can be fit}
  \item{posdep}{if TRUE, enforce positive dependence}
}
\value{Vector of delta function values}
\references{
Buckland, S. T., J. Laake and D. Borchers (2010). Double observer line transect methods: levels of independence. Biometrics 66: 169-177.
}
\author{Jeff Laake}
