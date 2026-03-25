\name{Imputation of missing values using the Jensen-Shannon based k--NN algorithm}
\alias{knnimp}
\title{
Imputation of missing values using the Jensen-Shannon based k--NN algorithm
}

\description{
Imputation of missing values using the Jensen-Shannon based k--NN algorithm.
}

\usage{
knnimp(x, k = 2)
}

\arguments{
\item{x}{
A matrix with the compositional data, zeros are allowed.
}
\item{k}{
The number of nearest neighbours to use. This can be a single number or a numerical vector.
}
}

\details{
Missing value imputation is performed using the k--NN algorithm with the Jensen-Shanon divergence.
}

\value{
A list including:
\item{pou}{
A vector with the rows (indexes) that contain missing values.
}
\item{y}{
A list with as many elements as the number of neighbours (length of k basically).
Each element of the list containts a matrix of the initial compositions, including the imputed missing values.
These correspond to each value of k.
}
}

\references{
Tsagris M., Alenazi A. and Stewart C. (2026). A Jensen--Shannon divergence based k--NN algorithm
for missing value imputation in compositional data.
Journal of Applied Statistics (Accepted for publication).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{acor}, \link{acor.tune}, \link{alfa}, \link{alfa.profile}
}
}

\examples{
x <- w <- rdiri(100, runif(7) )
dm <- dim(x)    ;    n <- dm[1]    ;    p <- dm[2]
xam <- matrix(nrow = 10, ncol =  3)
nu <- sample(1:n, 10)
for ( vim in 1:10)  xam[vim, ] <- sample(1:p, 3)
for ( vim in 1:10)  x[ nu[vim], xam[vim, ] ] <- NA
mod <- knnimp(x, k = 2)
}

