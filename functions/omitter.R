omitter <- function(data, desiredCols) { 
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
} # https://stackoverflow.com/questions/11254524/omit-rows-containing-specific-column-of-na
