\name{MLE of the Dirichlet distribution via Newton-Rapshon}
\alias{diria0.est}
\title{
MLE of the Dirichlet distribution via Newton-Rapshon
}

\description{
MLE of the Dirichlet distribution via Newton-Rapshon.
}

\usage{
diria0.est(x, tol = 1e-06)
}

\arguments{
\item{x}{
A matrix containing compositional data. Zeros are not allowed.
}
\item{tol}{
The tolerance level indicating no further increase in the log-likelihood.
}
}

\details{
Maximum likelihood estimation of the Dirichlet distribution, assuming a common parameter, is performed via Newton-Raphson.
}

\value{
A list including:
\item{iter}{
The number of iterations required.
}
\item{loglik}{
The value of the log-likelihood.
}
\item{a0}{
The estimated common parameter.
}
}

\references{
Thomas P. Minka (2003). Estimating a Dirichlet distribution.
http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{diri.est}, \link{diri.contour} \link{rdiri}, \link{ddiri}, \link{dda}
}
}

\examples{
x <- rdiri( 100, c(5, 5, 5, 5) )
diri.nr(x)
diria0.est(x)
}

