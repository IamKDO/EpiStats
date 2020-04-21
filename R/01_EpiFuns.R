S2 <- function(x) {
  sprintf("%3.2f", x)
}

fmt <- function(V) {
  as.numeric(as.character(V))
}

# =============================================================================
# Utilities for Epiconcept
# =============================================================================
toNumeric <- function(x, loop) {
  .T <- x
  for (i in 1:loop) {
    for(j in 1:2) {
      for(k in 1:2) {
        if (loop > 1) {
          .T[j,k, i] <- as.numeric(x[j,k, i])
        } else {
          .T[j,k] <- as.numeric(x[j,k])
        }
      }
    }
  }
  .T
}

rr <- function(Tb)
{

  TE = Tb[2,1]+Tb[2,2];
  TU = Tb[1,1]+Tb[1,2];
  CE = Tb[2,2];
  CU = Tb[1,2];
  TO = TE + TU;

  RE = CE / TE;
  RU = CU / TU;

  RR  = RE/RU;
  CI = computeRiskCI(RR, CE, TE, CU, TU);
  RRCIL = CI[1];
  RRCIH = CI[2];

  return(c(RR, RRCIL, RRCIH))
}

# rr2 <- function(Tb)
# {
#   TE = Tb[2,1]+Tb[2,2];
#   TU = Tb[1,1]+Tb[1,2];
#   CE = Tb[2,2];
#   CU = Tb[1,2];
#   TO = TE + TU;

#   X   The number of disease occurence among exposed cohort.
#   Y	  The number of disease occurence among non-exposed cohort.
#   m1  The number of individuals in exposed cohort group.
#   m2  The number of individuals in non-exposed cohort group.
#   conf.level  Probability for confidence intervals. Default is 0.95.

#   R <- riskratio(CE, CU, TE, TU, conf.level=0.95)
#   return(c(R$estimate, R$conf.int[1], R$conf.int[2]))
# }

# ========================================================================
# CASES CONTROLS STUDY
# ========================================================================
CC_AR <- function(C)
{
  T = epi.2by2(dat=C, method="case.control");
  S <- summary(T);
  return(c(S$AFest[1], S$AFest[2], S$AFest[3]));
}

CC_PAR <- function(C)
{
  .T = epi.2by2(dat=C, method="case.control", outcome="as.columns");
  S <- summary(.T);
  return(S$PAFest.strata.wald$est);
}

CC_STATS <- function(C)
{
  .T = epi.2by2(dat=C, method="case.control", outcome="as.columns");
  S <- summary(.T);
  return(S);
}

CS_STATS <- function(C)
{
  .T = epi.2by2(dat=C, method="cohort.count", outcome="as.columns");
  S <- summary(.T);
  return(S);
}

computeORCI <- function(OR, a,b,c,d)
{
  LNO = log(OR)
  R1 = sqrt((1/a)+(1/b)+(1/c)+(1/d))
  CIL = exp(LNO - 1.96 * R1)
  CIH = exp(LNO + 1.96 * R1)
  return(c(CIL, CIH))
}

computeExactORCI <- function(a,b,c,d)
{
  x <- matrix(c(a, b, c, d), 2, byrow = TRUE);
  R <- fisher.test(x);
  CIL <- R$conf.int[1];
  CIH <- R$conf.int[2];
  return(c(CIL, CIH, R$p.value));
}

computeDiffRiskCI <- function(RE, RU, NE, NU)
{
  A = RE - RU;
  B = (RE * (1-RE))/NE;
  C = (RU * (1-RU))/NU;
  D = 1.96*sqrt(B + C);
  R1 = A + D;
  R2 = A - D;

  return(c(R2, R1));
}

GetStrateVector <- function(A) {
  CE = A[2,2]    ; # Cases exposed
  CU = A[1,2]    ; # Cases unexposed
  HE = A[2,1]    ; # Healthy exposed
  HU = A[1,1]    ; # Healthy unexposed

  TE = CE + HE   ; # Total exposed
  TU = CU + HU   ; # Total unexposed
  TS = TE + TU   ; # Total strate
  H  = HU + HE   ; # Total healthy
  C  = CU + CE   ; # Total cases

  c(CE, CU, HE, HU, TE, TU, TS, H, C)
}

MANTEL_RR <- function(M) {
  colnames(M) <- c("CE", "CU", "HE", "HU", "TE", "TU", "TS", "H", "C")
  df <- data.frame(M)

  R1 = sum((df$CE * df$TU) / df$TS)
  R2 = sum((df$CU * df$TE) / df$TS)
  rrmh = R1 / R2

  R <- vector()
  for(I in 1:nrow(df)) {
    d2 = df[I, "TS"]^2
    V = ((df[I,"C"]*df[I,"TE"]*df[I,"TU"]) - (df[I,"CE"]*df[I,"CU"]*df[I,"TS"])) / d2
    R <- c(R, V)
  }

  NUMER = sum(R)
  DENOM = R1 * R2
  RES = sqrt(NUMER / DENOM)

  L = log(rrmh) - (1.96 * RES);
  H = log(rrmh) + (1.96 * RES);

  CIL = exp(L);
  CIH = exp(H);

  c(rrmh, CIL, CIH)
}

CMHrr <- function(A, B)
{
  # Stratum 1 =====================
  Ce1 = A[2,2]    ; # Cases exposed
  Cu1 = A[1,2]    ; # Cases unexposed
  He1 = A[2,1]    ; # Healthy exposed
  Hu1 = A[1,1]    ; # Healthy unexposed

  Te1 = Ce1 + He1 ; # Total exposed
  Tu1 = Cu1 + Hu1 ; # Total unexposed
  T1  = Te1 + Tu1 ; # Total strate 1
  H1  = Hu1 + He1 ; # Total healthy
  C1  = Cu1 + Ce1 ; # Total cases

  # Stratum 2 =====================
  Ce2 = B[2,2]    ; # Cases exposed
  Cu2 = B[1,2]    ; # Cases unexposed
  He2 = B[2,1]    ; # Healthy exposed
  Hu2 = B[1,1]    ; # Healthy unexposed

  Te2 = Ce2 + He2 ; # Total exposed
  Tu2 = Cu2 + Hu2 ; # Total unexposed
  T2  = Te2 + Tu2 ; # Total strate 2
  H2  = Hu2 + He2 ; # Total healthy
  C2  = Cu2 + Ce2 ; # Total cases



  R1 = ((Ce1 * Tu1) / T1) + ((Ce2 * Tu2) / T2);
  R2 = ((Cu1 * Te1) / T1) + ((Cu2 * Te2) / T2);
  rrmh = R1 / R2;

  R3 = ((C1*Te1*Tu1) - (Ce1*Cu1*T1)) / T1^2;
  R4 = ((C2*Te2*Tu2) - (Ce2*Cu2*T2)) / T2^2;
  R5 = R3 + R4;
  R6 = R5 / (R1 * R2);
  R7 = sqrt(R6);

  L = log(rrmh) - (1.96 * R7);
  H = log(rrmh) + (1.96 * R7);

  CIL = exp(L);
  CIH = exp(H);

  return(c(rrmh, CIL, CIH));
}

MH_HomogeneityTest <- function(mht)
{
  T = epi.2by2(dat=mht);
  print(summary(T))
  S <- summary(T);
  return(c(S$RR.homog[1], S$RR.homog[3]));
}

# computeKHI2 <- function(A, B, C, D)
# {
#   t <- chisq.test(matrix(c(A,B,C,D),ncol=2), correct=FALSE);
#   return(c(t$statistic, t$p.value));
# }

getVarsList <- function(list) {
  input <- as.list(list)
  if (length(input) < 2) {
    return(NULL)
  }
  .L <- length(input)
  input <- input[2:.L]
  .mode <- mode(input[[1]])
  # ---- .mode == "name" : we have one or more variable
  # ---------------------------------------------------
  if (.mode == "name") {
    ret <- c()
    for (v in input) {
      ret <- c(ret, as.character(v))
    }
#    cat("1 ===> ", ret, "\n")
    return(ret)
  }
  if (.mode == "call") {
    ret <- eval(as.vector(input[[1]]))
#    cat("2 ===> ", ret, "\n")
    return(ret)
  }

  if (.mode == "character" | .mode == "numeric") {
    ret <- c()
    for (v in input) {
      ret <- c(ret, v)
    }
#    cat("3 ===> ", ret, "\n")
    return(ret)
  }
}


