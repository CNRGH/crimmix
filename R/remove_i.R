#' Remove observation i from data
#'
#' @param data List of matrices.
#' @param i observation to remove
#'
#' @return data set without observation i
#' @export
#'
remove_i <- function (data, i){
  lapply(data, function (dd){
    dd[-i, ]
  })
}
