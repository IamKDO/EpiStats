regLog <- function(data, outcome, exposure) {
  .df <- data
  .cases <- as.character(substitute(outcome))
  .exp <- substitute(exposure)

  .cmd <- paste("glm(", .cases, " ~ ", .exp, ", data = .df, family = binomial(logit))")
  .ret <- eval(parse(text=.cmd))

  print(summary(.ret))

  print(autoplot(.ret, which = 1:6, label.size = 3) + theme_bw())

  .ret
}

