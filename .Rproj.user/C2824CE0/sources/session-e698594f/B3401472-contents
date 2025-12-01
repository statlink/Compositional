\name{Pyramid diagram}
\alias{pyramid}

\title{
Pyramid diagram
}

\description{
Pyramid diagram.
}

\usage{
pyramid(x, colour = "red", point_size = 4, face_alpha = 0.15)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{colour}{
If you want the points to appear in different colour put a vector with the colour
numbers or colours.
}
\item{point_size}{
The size of the plotted points. Note that the points appear as boxes.
}
\item{face_alpha}{
This controls the transparency of the tetrahedron faces:
When face_alpha= 0, this is completely transparent (invisible faces). When face_alpha=1, this is completely opaque (solid faces). When face_alpha = 0.15 (the default value) it is very transparent (you can see through the faces easily).
}
}

\details{
This is the pyramid that visualizes compositional data with 4 components in the 3D space.
}

\value{
The pyramid plot.
}

\references{
Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ternary}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
pyramid(x)
}
