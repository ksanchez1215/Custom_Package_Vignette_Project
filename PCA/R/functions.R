#' @title pca_tab
#' @description Performs principal components analysis on a data matrix and returns the result in tabular form
#' @param x a numeric or complex data frame that provides the data for the analysis
#' @param retx a logical value indicating whether rotated variables should be returned
#' @param center a logical value indicating if variables should be shifted to be zero centered
#' @param scale. a logical value indicating whether variables should be scaled to have unit variance before analysis
#' @param tol a value indicating the magnitude below which components should be omitted
#' @keywords PCA components
#' @export
#' @examples
#' pca_tab(x = iris[1:4], scale. = TRUE)

pca_tab <- function(x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, ...)
{
  chkDots(...)
  x <- as.matrix(x)
  x <- scale(x, center = center, scale = scale.)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")
  if(any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance")
  s <- svd(x, nu = 0)
  s$d <- s$d / sqrt(max(1, nrow(x) - 1))
  if (!is.null(tol)) {
    rank <- sum(s$d > (s$d[1L]*tol))
    if (rank < ncol(x)) {
      s$v <- s$v[, 1L:rank, drop = FALSE]
      s$d <- s$d[1L:rank]
    }
  }
  dimnames(s$v) <-
    list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
  r <- list(sdev = s$d, rotation = s$v,
            center = if(is.null(cen)) FALSE else cen,
            scale = if(is.null(sc)) FALSE else sc)
  if (retx) r$x <- x %*% s$v
  class(r) <- "prcomp"
  r
}

#' @title pca_plot
#' @description Plots the results of a principal components analysis is 2 dimensions
#' @param d a numeric or complex data frame that provides the data for the analysis
#' @param dataframe the entire dataframe that d came from, that can have variables as factors to group the results
#' @param groupby the variable from dataframe that you want to color your plot by
#' @param frame a logical value indicating whether your groups should have a geometric shape drawn around them
#' @param frame.type the shape of the frame, can be 'norm', 't', or 'euclid'. Default is 'norm' which creates an ellipse around the groups
#' @keywords PCA plot graph
#' @export
#' @examples
#' pca_plot(d = network_data[,10:130], dataframe = network_data, groupby = 'History3', frame = TRUE, frametype = 'norm')

pca_plot<- function(d, dataframe, groupby, frame = TRUE, frametype = 'norm'){
  temp<- pca_temp(d)
  plot<- ggplot2::autoplot(temp, data = dataframe, colour = groupby, frame = frame, frame.type = frametype)
  return(plot)
}
