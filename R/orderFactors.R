orderFactors <- function(data, ..., values, labels=NULL) {
  .df <- data
  
  .vars <- getVarsList(substitute(list(...)))
  if (is.null(.vars)) {
    msg = "You must provide one or more column(s) name(s)"
    stop(msg)
  }
  
  for (v in .vars) {
    if (!is.null(labels)) {
      .df[,v] <- factor(.df[,v], levels=values, ordered = TRUE, labels=labels)
    } else {
      .df[,v] <- factor(.df[,v], levels=values, ordered = TRUE)
    }
  }
  .df
}
