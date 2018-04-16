#make colors repeat
#implement that sweet Tufte graphics package

#######plotting differences
plotdiff <- function(x, treatcol = NULL, panel.eq = NULL, dodge = 0,
  x_axis = "time", pscale = FALSE, scales = "free_y") {

d <- x$design
if (is.null(d$times)) x_axis <- "treatment"
x2 <- x$treatment_comparisons

if (x_axis == "time" & is.null(treatcol)) treatcol <- treatcol_default(length(d$comparisons)+1)[-1]
if (x_axis == "treatment" & is.null(treatcol)) treatcol <- rep("black", (length(d$comparisons)))
#if treatcol does not = factor level, make it recycle
pd <- position_dodge(dodge)

gg <- ggplot(data = x2) +
  xlab(d$times) +
  ylab(paste(paste0(d$response, collapse = ", "), ": ",
    d$comp_function_name, sep = '')) +
  scale_color_manual(values = treatcol)

if (length(unique(treatcol)) < 2) gg <- gg + theme(legend.position="none")

if (x_axis == "time") {
  if (is.null(panel.eq)) panel.eq <- panel_default(c("response", d$panel))
  gg <- gg +
    geom_line(aes_string(d$times, "effect_size", col = "comparison"),
	    position = pd, lwd = 1) +
    geom_linerange(aes_string(d$times, ymax = "upr", ymin = "lwr",
	    col = "comparison"), position = pd, lwd = 1) +
    geom_hline(yintercept = c(0), col = "gray")
    }

if (x_axis == "treatment") {
  if (is.null(panel.eq) & is.null(d$times))
    panel.eq <- panel_default(c("response", d$panel))
  if (is.null(panel.eq) &! is.null(d$times))
    panel.eq <- panel_default(c("response", d$panel, d$times))
  gg <- gg + geom_pointrange(aes_string("comparison", "effect_size",
    ymax = "upr", ymin = "lwr", col = "comparison"), lwd = 1) +
  geom_hline(yintercept = c(0), col = "gray")
	}

if (!is.null(panel.eq)) gg <- gg + facet_grid(panel.eq, scales = scales)
gg
}
