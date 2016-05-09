plot.glmmadmb <- function(x, ...)
{
  plot(x$fitted, x$residuals, xlab="Fitted values", ylab="Residuals", ...)
}
