\name{lnl.io}
\alias{lnl.io}
\title{Mark-Recapture Distance Sampling: log-likelihood for IO configuration}
\description{Computes negative log-likelihood of independent observer configuration for double-observer point/line sampling}
\usage{
lnl.io(par, x1, x2, nfreq, models, width, cutp, debug, indep, PI, use.offset, posdep, point)
}
\arguments{
  \item{par}{values for parameters}
  \item{x1}{dataframe for observer 1}
  \item{x2}{dataframe for observer 2}
  \item{nfreq}{frequency of observations for each row in x}
  \item{models}{list of formulae for p(y) (p.formula) and delta(y) (delta.formula)}
  \item{width}{radius of point transect or half-width of line transect}
  \item{cutp}{cut points for distances}
  \item{debug}{if TRUE, output given during iteration}
  \item{indep}{if TRUE, full independence assumed}
  \item{PI}{if TRUE, point independence assumed; if !indep & !PI -- limiting independence }
  \item{use.offset}{should always be default of TRUE}
  \item{posdep}{if TRUE enforces positive dependence}
  \item{point}{if TRUE for point counts; otherwise, line transect}
}
\details{}
\value{}
\author{Jeff Laake}
