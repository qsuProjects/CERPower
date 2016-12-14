expand_matrix = function(.matrix, .obs) {
  library(plyr)
  cat("Expanding a matrix\n")
  
  .n = nrow(.matrix)
  .expanded = matrix(c(NA), nrow = .n*.obs, ncol = ncol(.matrix) )
  
  adply(.matrix, 1, function(..subject, ..obs) {
    matrix(rep(..subject, ..obs), nrow = ..obs, byrow = TRUE)
  }, .obs)
  
}
