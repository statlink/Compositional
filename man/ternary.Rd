\name{Ternary diagram}
\alias{ternary}

\title{
Ternary diagram
}

\description{
Ternary diagram.
}

\usage{
ternary(x, dg = FALSE, hg = FALSE, means = TRUE, pca = FALSE, colour = NULL)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{dg}{
Do you want diagonal grid lines to appear? If yes, set this TRUE.
}
\item{hg}{
Do you want horizontal grid lines to appear? If yes, set this TRUE.
}
\item{means}{
A boolean variable. Should the closed geometric mean and the arithmetic mean appear
(TRUE) or not (FALSE)?.
}
\item{pca}{
Should the first PCA calculated Aitchison (1983) described appear? If yes, then this
should be TRUE, or FALSE otherwise.
}
\item{colour}{
If you want the points to appear in different colour put a vector with the colour
numbers or colours.
}
}

\details{
There are two ways to create a ternary graph. We used here that one where each edge is equal to 1
and it is what Aitchison (1986) uses. For every given point, the sum of the distances from the edges
is equal to 1. Horizontal and or diagonal grid lines can appear, so as the closed geometric and the
simple arithmetic mean. The first PCA is calculated using the centred log-ratio transformation
as Aitchison (1983, 1986) suggested. If the data contain zero values, the first PCA will not be
plotted. Zeros in the data appear with green circles in the triangle and you will also see NaN in
the closed geometric mean.
}

\value{
The ternary plot and a 2-row matrix with the means. The closed geometric and the simple arithmetic
mean vector and or the first principal component will appear as well if the user has asked for them.
Additionally, horizontal or diagonal grid lines can appear as well.
}

\references{
Aitchison, J. (1983). Principal component analysis of compositional data. Biometrika 70(1): 57--65.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ternary.mcr}, \link{ternary.reg}, \link{diri.contour}
}
}

\examples{
x <- as.matrix(iris[, 1:3])
x <- x / rowSums(x)
ternary(x, means = TRUE, pca = TRUE)
}

\keyword{ ternary plot }
\keyword{ visualisation }
