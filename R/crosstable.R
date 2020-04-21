crossTable <- function(data, var1, var2, percent = "none", statistic = "none") {

  .cases <- as.character(substitute(var1))
  .by <- as.character(substitute(var2))
  # .percent <-ifelse(is.null(percent), "none", as.character(substitute(percent)))
  # .stats <-ifelse(is.null(statistic), "none", as.character(substitute(statistic)))
  .percent <- as.character(substitute(percent))
  .stats <- as.character(substitute(statistic))

  .Cases <- data[, .cases]
  .By <- data[, .by]

  if (!is.factor(.Cases)) {
    .Cases <- as.factor(.Cases)
  }

  if (!is.factor(.By)) {
    .By <- as.factor(.By)
  }

  # --- Test if we have ordered factors
  .ordered <- TRUE
  .ord.cases <- ("ordered" %in% class(.cases))
  .ord.expos <- ("ordered" %in% class(.by))
  if ((.ord.cases & .ord.expos) == FALSE) {

    .ordered <- FALSE
  }

  .Title <- sprintf("%s / %s", .by, .cases)
  .ColExpo <- c(levels(.By))
  .ColCases <- c(levels(.Cases))

  DF <- NULL

  # Contingency table
  # ===========================================================================
  FR = table(.By, .Cases)
  I1E1 = FR[1,1] <- as.numeric(FR[1,1])
  I1E0 = FR[1,2] <- as.numeric(FR[1,2])
  I0E0 = FR[2,2] <- as.numeric(FR[2,2])
  I0E1 = FR[2,1] <- as.numeric(FR[2,1])

  # Totaux
  # ===========================================================================
  TE1 <- I1E1 + I0E1
  TE0 <- I1E0 + I0E0
  TI1 <- I1E1 + I1E0
  TI0 <- I0E1 + I0E0

  if (.percent == "none") {
    col1 <- c(.ColExpo, "Total")
    col2 <- c(I1E1, I0E1, TE1)
    col3 <- c(I1E0, I0E0, TE0)
    col4 <- c(TI1, TI0, TI1+TI0)


    DF <- data.frame(col1, col2, col3, col4, stringsAsFactors = FALSE)
    colnames(DF) <- c(.Title, .ColCases, "Total")
    #print(DF)
  }

  # ---- Compute vars
  # --------------------------------------------------------------------
  TI1I0 <- TI1+TI0
  PI1E1 <- I1E1/TE1 * 100
  PI0E1 <- I0E1/TE1 * 100
  PI1E0 <- I1E0/TE0 * 100
  PI0E0 <- I0E0/TE0 * 100
  PT11  <- TI1/TI1I0 * 100
  PT10  <- TI0/TI1I0 * 100
  PE1T1 <- I1E1/TI1 * 100
  PE0T1 <- I1E0/TI1 * 100
  PE1T0 <- I0E1/TI0 * 100
  PE0T0 <- I0E0/TI0 * 100
  TT <- TI1+TI0

  if (.percent == "col") {
    col1 <- c(.ColExpo[1], "%", .ColExpo[2], "%", "Total", "%")
    col2 <- c(I1E1, S2(PI1E1), I0E1, S2(PI0E1), TE1, S2(PI1E1+PI0E1))
    col3 <- c(I1E0, S2(PI1E0), I0E0, S2(PI0E0), TE0, S2(PI1E0+PI0E0))
    col4 <- c(TI1, S2(PT11), TI0, S2(PT10), TI1+TI0, S2(PT11+PT10))

    DF <- data.frame(col1, col2, col3, col4, stringsAsFactors = FALSE)
    colnames(DF) <- c(.Title, .ColCases, "Total")

  }

  if (.percent == "row") {
    col1 <- c(.ColExpo[1], .ColExpo[2], "Total")
    col2 <- c(I1E1, I0E1, TE1)
    col3 <- c(S2(PE1T1), S2(PE1T0), S2((TE1/TT)*100))
    col4 <- c(I1E0, I0E0, TE0)
    col5 <- c(S2(PE0T1), S2(PE0T0), S2((TE0/TT)*100))
    col6 <- c(TI1, TI0, TT)
    col7 <- rep(100.00, 3)
    DF <- data.frame(col1, col2, col3, col4, col5, col6, col7, stringsAsFactors=FALSE)
    colnames(DF) <- c(.Title, .ColCases[1], "%", .ColCases[2], "%", "Total", "%")
  }

  if (.percent == "both") {
    col1 <- c(.ColExpo[1], "%", .ColExpo[2], "%", "Total", "%")
    col2 <- c(I1E1, S2(PI1E1), I0E1, S2(PI0E1), TE1, S2(PI1E1+PI0E1))
    col3 <- c(S2(PE1T1), "-", S2(PE1T0), "-",  S2((TE1/TT)*100), "-")
    col4 <- c(I1E0, S2(PI1E0), I0E0, S2(PI0E0), TE0, S2(PI1E0+PI0E0))
    col5 <- c(S2(PE0T1), "-", S2(PE0T0), "-", S2((TE0/TT)*100), "-")
    col6 <- c(TI1, "-", TI0, "-", TT, "100.00")
    col7 <- c("100.00", "-", "100.00", "-", "100.00", "-")

    DF <- data.frame(col1, col2, col3, col4, col5, col6, col7, stringsAsFactors=FALSE)

    colnames(DF) <- c(.Title, .ColCases[1], "%", .ColCases[2], "%", "Total", "%")

  }

  if (.stats == "none") {
    return(DF)
  }


  if (.stats == "chi2") {
    .nbcol <- ncol(DF)
    .cmp <- rep("-", .nbcol-4)
    S = computeKHI2(I1E1, I0E1, I1E0, I0E0);
    DF <- rbind(DF, rep("-", .nbcol))
    DF <- rbind(DF, c("Pearson CHI2", round(S[1], 4), "Pr", round(S[2], 3), .cmp))
    return(DF)
  }

  if (.stats == "fisher") {
    .nbcol <- ncol(DF)
    .cmp <- rep("-", .nbcol-2)
    S = computeFisher(I1E1, I0E1, I1E0, I0E0);
    DF <- rbind(DF, rep("-", .nbcol))
    DF <- rbind(DF, c("Fisher's exact", round(S[1], 3), .cmp))
    return(DF)
  }

}
