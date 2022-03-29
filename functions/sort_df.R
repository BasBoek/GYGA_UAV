# sort_df: Putting a series of columns in front of the table
sort_df <- function(data, colnames){
  for(i in length(colnames):1){
    if(i == length(colnames)){
      new_df <- data[,c(which(colnames(data)==colnames[i]), which(colnames(data)!=colnames[i]))]
    } else {      new_df <- new_df[,c(which(colnames(new_df)==colnames[i]), which(colnames(new_df)!=colnames[i]))]    }}
  return(new_df)
}
