
cstable <- CSTable <- function(x,
                    cases,
                    exposure=c(),
                    exact=FALSE,
                    sort = "pvalue",
                    full = FALSE
) UseMethod("CSTable", x)

CSTable.data.frame <- function(x,
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
  
  Colnames = c("Tot.Exp.", "Exp.Cases", "AR%", "Tot.Unex.", "Unex.Cases", "AR%",
                       "RR", "CI ll", "CI ul", PLabel);
  
  .TotalExposed <- c()
  .CasesExposed <- c()
  .TotalUnexposed <- c()
  .CasesUnexposed <- c()
  .ARExposed <- c();
  .ARUnexposed <- c()
  .RiskRatio <- c()
  .RiskCILow <- c()
  .RiskCIHight <- c()
  .Pvalue <- c()
  
  
  for (N in exposure) {
    # print(sprintf("%s", N))
    
    FR = table(.Cases, x[,N])
    I1E1 = FR[2,2]
    I1E0 = FR[2,1]
    I0E0 = FR[1,1]
    I0E1 = FR[1,2]

    if (exact == TRUE) {
      .Stat <- computeFisher(FR[1,1], FR[2,1], FR[1,2], FR[2,2]);
    } else {
      .Stat <- computeKHI2(FR[1,1], FR[2,1], FR[1,2], FR[2,2])[2];
    }
    
    .Stat <- as.numeric(sprintf("%3.6f", .Stat))
    
    .Pvalue <- c(.Pvalue, .Stat)
    
    # Compute Total
    # ---------------------------------------------------------------------------
    TE = I0E1 + I1E1        ; # Total exposed
    TU = I0E0 + I1E0        ; # Total unexposed
    TCA = I1E1 + I1E0       ; # Total cases
    TNC = I0E1 + I0E0       ; # Total non-cases
    
    P1 = (I1E1/TE)*100;
    P0 = (I1E0/TU)*100;
    
    
    .TotalExposed <- c(.TotalExposed, TE)
    .CasesExposed <- c(.CasesExposed, I1E1)
    .TotalUnexposed <- c(.TotalUnexposed, TU)
    .CasesUnexposed <- c(.CasesUnexposed, I1E0)
    .ARExposed <- c(.ARExposed, P1);
    .ARUnexposed <- c(.ARUnexposed, P0)
    
    RR = P1/P0;
    .RiskRatio = c(.RiskRatio, RR);

    CI = computeRiskCI(RR, FR[2,2], TE, FR[2,1], TU);
    .RiskCILow   = c(.RiskCILow, CI[1]);
    .RiskCIHight = c(.RiskCIHight, CI[2]);

  }

  DF <- data.frame(.TotalExposed,
                   .CasesExposed,
                   S2(.ARExposed),
                   .TotalUnexposed,
                   .CasesUnexposed,
                   S2(.ARUnexposed),
                   S2(.RiskRatio),
                   S2(.RiskCILow),
                   S2(.RiskCIHight),
                   round(.Pvalue,3),
                   row.names = exposure
  );
  
  colnames(DF) <- Colnames

  if (sort == "rr") {
    DF <- DF[order(-fmt(DF[,7])),];
  } else if (sort == "ar") {
    DF <- DF[order(-fmt(DF[,3])),];
  } else if (sort == "pvalue") {
    DF <- DF[order(DF[,10]),];
  }
  
  if (full == TRUE) {
    ret <- list(df = DF, digits=c(0,0,1,0,0,1,2,2,2,3), align="ccrccrrrrr")
  } else {
    ret <- list(df = DF)
  }

  ret
  
}
