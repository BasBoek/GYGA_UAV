# Function to add zero to single numbers
add_zero <- function(df, colname){
  NR_ID <- as.character(df[, colname])
  for(i in 1:length(NR_ID)){
    if(nchar(NR_ID[i]) == 1){      NR_ID[i] <- paste0("0", NR_ID[i])    }
  }
  result <- cbind(df, NR_ID)
  return(result)
}
