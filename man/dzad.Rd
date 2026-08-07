\name{Density values of the zero adjusted Dirichlet distribution}
\alias{dzad}
\title{
Density values of the zero adjusted Dirichlet distribution
}

\description{
Density values of the zero adjusted Dirichlet distribution.
}

\usage{
dzad(y, phi, mu, logged = TRUE)
}

\arguments{
\item{y}{
A matrix with compositional data with zero values.
}
\item{phi}{
The concentration parameter.
}
\item{mu}{
The mean vector.
}
\item{logged}{
A boolean variable specifying whether the logarithm of the density values to
be returned. It is set to TRUE by default.
}
}

\details{
The density values of the zero adjusted Dirichlet distribution are computed.
}

\value{
A vector with the density values.
}

\references{
Tsagris M. and Stewart C. (2018). A Dirichlet regression model for compositional data with zeros.
Lobachevskii Journal of Mathematics, 39(3): 398--412.

Preprint available from https://arxiv.org/pdf/1410.5011.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{zad.est}, \link{zadr}
}
}

\examples{
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
y[sample(1:450, 15) ] <- 0
mod2 <- zad.est(y)
f <- dzad(y, mod2$phi, mod2$mu)
}

