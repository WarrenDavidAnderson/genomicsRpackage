% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_functions.R
\name{odemodel}
\alias{odemodel}
\title{Dynamic model}
\usage{
odemodel(time, y, params)
}
\arguments{
\item{time}{Vector of simulation times.}

\item{y}{Dynamic variables.}

\item{params}{Parameter matrix, the Jacobian.}
}
\value{
Model object for input into the ode() function from the deSolve library.
}
\description{
This function generates the mathematical formulation for the dynamic model.
}
\examples{
out <- ode(y=ic, times=time, parms=J, func=odemodel, method=lsoda)
}
\references{
\url{https://www.jstatsoft.org/article/view/v033i09}
}
