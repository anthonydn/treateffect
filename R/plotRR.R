#######plotting differences
plotRR <- function(x, panel.eq = NULL, scales = "free_y", horizontal = TRUE) {

d <- x$design
x2 <- x$comparisons
x2$response <- fct_rev(x2$response)

gg <- ggplot(data = x2) +
  ylab("Log response ratio")

  if (is.null(panel.eq) & is.null(d$times))
    panel.eq <- panel_default(c("comparison", d$panel))
  if (is.null(panel.eq) &! is.null(d$times))
    panel.eq <- panel_default(c("comparison", d$panel, d$times))
  gg <- gg + geom_pointrange(aes_string("response", "effect_size",
    ymax = "upr", ymin = "lwr"), lwd = 1) +
    geom_hline(yintercept = 1, lwd = 1, col = "red") +
    scale_y_log10(minor_breaks = RR_grid(c(x2$lwr,x2$upr))) +
    theme(panel.grid.minor = element_line(colour = "gray"),
      panel.grid.major = element_line(colour = "gray"),
      axis.title.y = element_blank())

if (!is.null(panel.eq)) gg <- gg + facet_grid(panel.eq, scales = scales)
h <- dim(x2)[1] #/ prod(unlist(lapply(d$panel, function(x) length(levels(x2[[x]])))))
if (horizontal) gg <- gg + coord_flip()
gg
}

RR_grid <- function(e) {
  ra <- floor(range(log10(e)))
  gr <- as.numeric(outer(1:9,10^c(ra[1]:ra[2])))
  i <- which(gr > min(e) & gr < max(e))
  gr[c(i-1,i,i+1)]}
