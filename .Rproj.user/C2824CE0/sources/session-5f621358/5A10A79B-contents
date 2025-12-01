pyramid <- function(x, colour = "red", point_size = 4, face_alpha = 0.15) {
    
  labels <- colnames(x)
  if ( is.null(labels) )  labels <- LETTERS[1:4]
  vertices <- matrix(c( 0, 0, 0, 1, 0, 0, 0.5, sqrt(3)/2, 0, 0.5, sqrt(3)/6, sqrt(2/3)), ncol = 3, byrow = TRUE)
  coords_3d <- x %*% vertices
  
  # Plot
  rgl::open3d()
  rgl::bg3d("white")
  
  # Edges
  edges <- rbind( c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4) )
  for ( i in 1:6 ) rgl::segments3d( rbind( vertices[edges[i, 1], ], vertices[edges[i, 2], ]), col = "black", lwd = 2)
  # Faces
  faces <- rbind( c(1,2,3), c(1,2,4), c(1,3,4), c(2,3,4) )
  for ( i in 1:4 )  rgl::triangles3d( vertices[faces[i,],], col = "lightblue", alpha = face_alpha, lit = FALSE)
    
  rgl::points3d(coords_3d, col = point_col, size = point_size)
  # Vertices
  rgl::points3d(vertices, col = "blue", size = 1)
  rgl::text3d(vertices, texts = labels, col = "blue", cex = 1.5, adj = c(0.5, -0.5))
  rgl::view3d(theta = 30, phi = 20, zoom = 0.8)
  invisible(coords_3d)
}


# Example
set.seed(123)
x<- matrix(runif(50 * 4), ncol = 4)
x[1, 3] <- 0
x <- x / rowSums(x)
colnames(x) <- c("A", "B", "C", "D")

pyramid(x)