orderbykey <- function(df) {
  df$Reference.Key. <- as.character(df$Reference.Key.)
  return(df[order(df$Reference.Key.),])
}

diffcols <- function(df1, df2, colname1, colname2) {
  diff <- df1[, c("Reference.Key.", colname1)]
  diff <- merge(diff, df2[, c("Reference.Key.", colname2)], by="Reference.Key.", all.x=T)
  return(diff)
}
