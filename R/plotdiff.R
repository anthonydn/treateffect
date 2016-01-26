#######plotting differences
plotdiff <- function(x, treatcol = NULL, panel.eq = NULL, dodge = 0,
  x_axis = "time", pscale = FALSE) {

d <- attributes(x$data)$design
if (is.null(d$times)) x_axis <- "treatment"
x2 <- x$analyses
nn <- names(d$extract)

if (length(nn) == 2) {
  x2$center <- x2[,nn[1]]
  x2$lo <- x2[,nn[1]]-x2[,nn[2]]
  x2$hi <- x2[,nn[1]]+x2[,nn[2]]}
if (length(d$extract)==3) {
  x2$center <- x2[,nn[1]]
  x2$lo <- x2[,nn[2]]
  x2$hi <- x2[,nn[3]]}
#does not currently work with "fracdiff" type functions because it is looking for 0-100
#needs to be rethought though to deal with 0s and infinities
if (pscale) {
  x2$center <- pscale(x2$center)
  x2$lo <- pscale(x2$lo)
  x2$hi <- pscale(x2$hi)}
if (is.null(treatcol)) treatcol <- treatcol.default(length(d$comparisons)+1)[-1]
pd <- position_dodge(dodge)

gg <- ggplot(data = x2) +
  xlab(d$times) +
  ylab(paste(d$response,", ", d$extract[1], " from ", d$control, " ± ",
	  d$extract[2], sep = '')) +
  scale_color_manual(values = treatcol)

if (x_axis == "time") {
  if (is.null(panel.eq)) panel.eq <- panel.default(d$panel)
  gg <- gg +
    geom_line(aes_string(d$times, "center", col = "comparison"),
	    position = pd, lwd = 1) +
    geom_linerange(aes_string(d$times, ymax = "hi", ymin = "lo",
	    col = "comparison"), position = pd, lwd = 1) +
    geom_hline(yintercept = c(0), col = "gray")
    }

if (x_axis == "treatment") {
  if (is.null(panel.eq) & is.null(d$times))
    panel.eq <- panel.default(c(d$panel))
  if (is.null(panel.eq) &! is.null(d$times))
    panel.eq <- panel.default(c(d$panel, d$times))
  gg <- gg + geom_pointrange(aes_string("comparison", "center",
    ymax = "hi", ymin = "lo", col = "comparison"), lwd = 1) +
  geom_hline(yintercept = c(0), col = "gray")
	}

if (!is.null(panel.eq)) gg <- gg + facet_grid(panel.eq,scales="free_y")
gg
}
