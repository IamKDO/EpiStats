cctable <- CCTable <- function(x,
                    cases,
                    exposure=c(),
                    exact=FALSE,
                    sort = "pvalue",
                    full = FALSE
) UseMethod("CCTable", x)

CCTable.data.frame <- function(x,
                               cases,
                               exposure = c(),
                               exact = FALSE, 
                               sort = "pvalue",
                               full = FALSE
)
{
  .Cases <- x[, cases]
  
  if (length(exposure) < 1) {
    stop("Exposure list is empty.");
  }
  
  PLabel = ifelse(exact == TRUE, "p(Fisher)", "p(Chi2)")
  

  .TotalCases <- c()
  .TotalExposed <- c()
  .CasesExposed <- c()
  .CtrlExposed <- c()
  .TotalUnexposed <- c()
  .CasesUnexposed <- c()
  .TotalCtrl <- c()
  .OddsRatio <- c()
  .PCAExposed <- c()
  .PCTExposed <- c()
  # .ARExposed <- c();
  # .ARUnexposed <- c()
  # .RiskRatio <- c()
  .CILow <- c()
  .CIHight <- c()
  .Pvalue <- c()

  for (N in exposure) {
    # print(sprintf("%s", N))
    
    FR = table(.Cases, x[,N])
    # print(FR)
    KE  <- FR[2,2]    # Cases exposed
    KU  <- FR[2,1]    # Cases unexposed
    CTE <- FR[1,2]    # Controls exposed
    CTU <- FR[1,1]    # Controls unexposed
    TE  <- KE + CTE   # Total exposed
    TU  <- KU + CTU   # Total unexposed
    TCA <- KE + KU
    TCT = CTE + CTU
    P1 = (FR[2,2]/TCA)*100;
    P0 = (FR[1,2]/TCT)*100;
    
    if (exact == TRUE) {
      .Stat <- computeFisher(FR[1,1], FR[2,1], FR[1,2], FR[2,2]);
    } else {
      .Stat <- computeKHI2(FR[1,1], FR[2,1], FR[1,2], FR[2,2])[2];
    }
    
    #.Pvalue <- c(.Pvalue, .Stat)
    
    .Colnames = c("Tot.Cases", "Exposed", "%", "Tot.Ctrls", "Exposed", "%",
                         "OR", "CI ll", "CI ul", PLabel);
    
    .TotalCases   <- c(.TotalCases, TCA)
    .CasesExposed <- c(.CasesExposed, KE)
    .PCAExposed   <- c(.PCAExposed, P1)
    .TotalCtrl    <- c(.TotalCtrl, TCT)
    .CtrlExposed  <- c(.CtrlExposed, CTE)
    .PCTExposed   <- c(.PCTExposed, P0)
    
    ODR = (FR[2,2]/FR[2,1]) / (FR[1,2]/FR[1,1]);
#    RR = ODR
    
    .OddsRatio = c(.OddsRatio, ODR)
    
    .CI = computeExactORCI(CTU, CTE, KU, KE);
    .CILow   = c(.CILow, .CI[1])
    .CIHight = c(.CIHight, .CI[2])
    if (exact == TRUE) {
      .PVal = as.numeric(round(.CI[3], 5))
      .Pvalue = c(.Pvalue, .PVal)
    } else {
      .PVal = as.numeric(round(.Stat, 5))
      .Pvalue = c(.Pvalue, .PVal)
    }
  }

  DF <- data.frame(.TotalCases,
                   .CasesExposed,
                   S2(.PCAExposed),
                   .TotalCtrl,
                   .CtrlExposed,
                   S2(.PCTExposed),
                   S2(.OddsRatio),
                   S2(.CILow),
                   S2(.CIHight),
                   round(.Pvalue,3),
                   row.names = exposure
  )
  
  colnames(DF) <- .Colnames
  
  if (sort == "or") {
    DF <- DF[order(-fmt(DF[,7])),]
  } else if (sort == "pe") {
    DF <- DF[order(-fmt(DF[,3])),]
  } else if (sort == "pvalue") {
    DF <- DF[order(DF[,10]),]
  }
  
  if (full == TRUE) {
    ret <- list(df = DF, align="ccrccrrrrr")
  } else {
    ret <- list(df = DF)
  }
  ret
  
}
