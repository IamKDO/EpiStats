cc <- CC <- function( data,
                      cases,
                      exposure,
                      exact = FALSE,
                      full = FALSE,
                      title = "CC"
) UseMethod("CC", data)

CC.data.frame <- function(  data,
                            cases,
                            exposure,
                            exact = F,
                            full = FALSE,
                            title = "CC"
)
{
  .cases <- as.character(substitute(cases))
  .expos <- as.character(substitute(exposure))

  .Cases <- data[, .cases]
  .Exposure <- data[, .expos]
  .T1 = title
  .T2 = title

  Rownames1 <- c("Exposed", "Unexposed", "Total", "Proportion exposed")
  Colnames1 <- c("Cases", "Controls", "Total")
  Colnames2 <- c("Point estimate", "95%CI-ll", "95%CI-ul")
  PLabel <- c("chi2(1)", "Pr>chi2")


  # Contingency table
  # ===========================================================================
  FR = table(.Cases, .Exposure)
  I1E1 = FR[2,2] <- as.numeric(FR[2,2])
  I1E0 = FR[2,1] <- as.numeric(FR[2,1])
  I0E0 = FR[1,1] <- as.numeric(FR[1,1])
  I0E1 = FR[1,2] <- as.numeric(FR[1,2])

  STAT = computeKHI2(FR[1,1], FR[2,1], FR[1,2], FR[2,2]);

  FISHER <- NA
  if (exact == TRUE) {
    FISHER = computeFisher(FR[1,1], FR[2,1], FR[1,2], FR[2,2]);
    PLabel = c(PLabel, "Fisher p-value")
  }

  # Compute Total
  # ---------------------------------------------------------------------------
  TE = I0E1 + I1E1        ; # Total exposed
  TU = I0E0 + I1E0        ; # Total unexposed
  TCA = I1E1 + I1E0       ; # Total cases
  TNC = I0E1 + I0E0       ; # Total non-cases



  # Compute Proportions
  # ---------------------------------------------------------------------------
  PCAEX = I1E1/TCA
  PCTEX = I0E1/TNC
  PTOEX = TE/(TCA+TNC)

  COL_C <- as.character(c(I1E1, I1E0, TCA, sprintf("%2.2f",PCAEX)))
  COL_N <- as.character(c(I0E1, I0E0, TNC, sprintf("%2.2f",PCTEX)))
  COL_T <- as.character(c(TE, TU, TE+TU, sprintf("%2.2f",PTOEX)))
  #COL_R <- c(PCAEX, PCTEX, PTOEX)
  R1 <- data.frame(COL_C, COL_N, COL_T, row.names = Rownames1)
  colnames(R1) <- Colnames1


  # Estimate
  # ===========================================================================

  # ODD Ratio
  # ---------------------------------------------------------------------------

  R = or(FR);
  OREST = ODD <- R[1]
  ORCIL = R[2]
  ORCIH = R[3]


  rr = list(point_estimate = R, CI95.ll = ORCIL, CI95.ul = ORCIH)

  # Attr.frac.pop. <- Attr.frac.ex. x proportion.of.cases.exposed.
  # Prev.frac.pop. <- Prev.frac.ex. x proportion.of.controls.exposed
  if (R[1] >= 1.0) {
    R = CC_STATS(FR);
    # AFEST = as.numeric(R$AFest[1])

    AFEST = as.numeric((ODD-1) / ODD)
    AFCIL = as.numeric(R$AFest[2])
    AFCIH = as.numeric(round(R$AFest[3], 8))

    PAEST <- AFEST * PCAEX
    Rownames2 <- c("Odds ratio", "Attr. frac. ex.", "Attr. frac. pop", PLabel)

    RES <- c(point_estimate = AFEST, CI95.ll = AFCIL, CI95.ul = AFCIH)

    stats = list(odds_ratio=rr,
                 Attr.frac.ex = RES,
                 Attr.frac.pop = PAEST,
                 chi2 = as.numeric(STAT[1]),
                 p.chi2 = as.numeric(round(STAT[2], 8))    )
  }
  else {
    AFEST = as.numeric(1 - R[1])
    AFCIL = 1 - R[3]
    AFCIH = 1 - R[2]

    # Pe = TE / (TE + TU);
    # PAEST = Pe * (1 - R[1])
    PAEST <- AFEST * PCTEX
    Rownames2 <- c("Odds ratio", "Prev. frac. ex.", "Prev. frac. pop", PLabel)

    RES <- c(point_estimate = AFEST, CI95.ll = AFCIL, CI95.ul = AFCIH)

    stats = list(odds_ratio=rr,
                 Prev.frac.ex = RES,
                 Prev.frac.pop = PAEST,
                 chi2 = as.numeric(STAT[1]),
                 p.chi2 = as.numeric(round(STAT[2], 8))
    )

  }
  str_PCHI2 <- sprintf("%3.3f", as.numeric(round(STAT[2], 8)))
  str_FISH <- sprintf("%3.3f", round(FISHER, 5))

  .COL_PES <- c(S2(OREST), S2(AFEST), S2(PAEST), S2(as.numeric(STAT[1])), str_PCHI2)
  .COL_CIL <- c(ORCIL, AFCIL, NA, NA, NA)
  .COL_CIH <- c(ORCIH, AFCIH, NA, NA, NA)

  if (exact == TRUE) {
    .COL_PES <- c(.COL_PES, str_FISH)
    .COL_CIL <- c(.COL_CIL, NA)
    .COL_CIH <- c(.COL_CIH, NA)
  }

  R2 <- data.frame(.COL_PES, S2(.COL_CIL), S2(.COL_CIH), row.names = Rownames2)
  colnames(R2) <- Colnames2

  if (full == TRUE) {
    ret <- list(t1 = .T1, df1 = R1, df2 = R2,
                df1.align = c("rrr"),
                df2.align = c("rrr"),
                st = stats)
  } else {
    ret <- list(df1 = R1, df2 = R2)
  }

  ret

}
