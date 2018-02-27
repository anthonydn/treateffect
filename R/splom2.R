splom2 <- function (R, cutoff = 0) {
  x = PerformanceAnalytics::checkData(R, method = "matrix")
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
      cex.cor) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r2 <- summary(lm(x ~ y, data = R))$r.square
    ci <- r2bootci(na.omit(data.frame(x,y)), x~y)
    txt2 <- format(c(r2, 0.123456789), digits = digits)[1]
    txt <- paste(txt2, " [", as.character(round(ci, 2))[4], " - ", as.character(round(ci, 2))[5], "]", sep = '')

    if (r2>cutoff) text(0.5, 0.5, txt, cex = 1)}
  
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd(x))}

  hist.panel = function(x) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    abline(v = hotspots::hotspots(x)$positive.cut, col = "red")
    rug(x)}
  
  pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel)
}

r2bootci <- function(data, formula) {
  results <- boot::boot(data = data, statistic = rsq, R = 1000, formula = formula)
  boot::boot.ci(results, type="bca")$bca
  }

rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 

