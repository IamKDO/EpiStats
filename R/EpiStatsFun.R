# ===========================================================================
# Internal functions
# ===========================================================================

computeKHI2 <- function(A, B, C, D)
{
  t <- chisq.test(matrix(c(A,B,C,D),ncol=2), correct=FALSE);
  return(c(t$statistic, t$p.value));
}

computeFisher <- function(A, B, C, D)
{
  t <- fisher.test(matrix(c(A,B,C,D),ncol=2));
  t$p.value
}

computeRiskCI <- function(risk, X1, N1, X2, N2)
{
  A = ((N1-X1)/X1)/N1;
  B = ((N2-X2)/X2)/N2;
  R1 = log(risk) + (1.96*sqrt(A + B));
  R2 = log(risk) - (1.96*sqrt(A + B));
  E1 = exp(R1);
  E2 = exp(R2);
  
  return(c(E2, E1));
}

# CASES CONTROLS STUDY
# ========================================================================
CC.table <-function(x, cases, exposure, strate, level) {
  A  = table(x[strate == level, exposure], x[strate == level, cases])
  CE = A[2,2]    ; # Cases exposed
  CU = A[1,2]    ; # Cases unexposed
  HE = A[2,1]    ; # Healthy exposed
  HU = A[1,1]    ; # Healthy unexposed
  M <- matrix(c(CE,HE,CU,HU),ncol = 2, byrow = TRUE)
  as.table(M)
}

# Comppute ODDS ratio
# -----------------------------------------------------------------------------
or <- function(.T)
{
  O <- (.T[1,1]/.T[1,2]) / (.T[2,1]/.T[2,2]);
  x <- matrix(.T, 2, byrow = TRUE);
  R <- fisher.test(x);
  CIL <- R$conf.int[1];
  CIH <- R$conf.int[2];
  return(c(O, CIL, CIH));
}
# Comppute ODDS ratio
# -----------------------------------------------------------------------------
ODD.ratio <- function(.T) {
#  O <- (.T[1,1]/.T[1,2]) / (.T[2,1]/.T[2,2]);
  O <- (.T[1,1]/.T[2,1]) / (.T[1,2]/.T[2,2]);
  x <- matrix(.T, 2, byrow = TRUE);
  R <- fisher.test(x);
  CIL <- R$conf.int[1];
  CIH <- R$conf.int[2];
  return(c(O, CIL, CIH));
}
