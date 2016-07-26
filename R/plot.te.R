#main plot function
plot.te <- function(x, treatcol = NULL, panel.eq=NULL, dodge = 0,
  x_axis = "time", points = TRUE, cen = "mean", bars = "se", scales = "free_y", subset = NULL) {

#can't turn bars off
# maybe spend the ... on parameters for geom point

x2 <- x$treatment_summaries
subset <- eval(substitute(subset), x2)
if (!is.null(subset)) x2 <- x2[subset,]

d <- x$design
if (is.null(d$times)) x_axis <- "treatment"
if (is.null(treatcol)) treatcol <- treatcol.default(length(d$levels))
pd <- position_dodge(dodge)

x2$cen <- x2[[cen]]
if (length(bars) == 1) x2$lo <- x2[[cen]] - x2[[bars]]
if (length(bars) == 1) x2$hi <- x2[[cen]] + x2[[bars]]
if (length(bars) == 2) x2$lo <- x2[[bars[1]]]
if (length(bars) == 2) x2$hi <- x2[[bars[2]]]

gg <- ggplot(data = x2) +
  xlab(d$times) +
  ylab(paste(paste0(d$response, collapse = ", "), ": ", cen, " \u00B1 ", bars)) +
  scale_color_manual(values = treatcol)

if (x_axis == "time") {
  if (is.null(panel.eq)) panel.eq <- panel.default(c("variable", d$panel))
  gg <- gg + geom_line(aes_string(d$times, "cen", col = d$treatment),
      position = pd, lwd = 1) +
    geom_linerange(aes_string(d$times, "cen", ymax = "hi",
	  ymin = "lo", col = d$treatment), position = pd, lwd = 1)
  if (points) gg <- gg + geom_point(aes_string(d$times, "response",
    col = d$treatment), data = x$data, position = pd)
	#add if flag and ribbons
  }

if (x_axis == "treatment") {
  if (is.null(panel.eq) & is.null(d$times))
    panel.eq <- panel.default(c("variable", d$panel))
  if (is.null(panel.eq)&!is.null(d$times))
    panel.eq <- panel.default(c("variable", d$panel,"times"))
  if (bars == "box") gg <- gg +
    geom_boxplot(aes_string(d$treatment,"response", col = d$treatment),
	  data = x$data) else
    gg <- gg + geom_linerange(aes_string(d$treatment, "cen", ymax = "hi",
    ymin = "lo", col = d$treatment), position = pd, lwd = 1)
  if (points) gg <- gg + geom_point(aes_string(d$treatment, "response",
    col = d$treatment), data = x$data, position = pd)
  x2$gs1 <- as.numeric(x2[[d$treatment]]) - 0.3
  x2$gs2 <- as.numeric(x2[[d$treatment]]) + 0.3
  gg <- gg + geom_segment(aes_string("gs1", "cen", yend = "cen", xend = "gs2",
	col = d$treatment), data = x2, lwd = 1)
	}

if (!is.null(panel.eq)) gg <- gg + facet_grid(panel.eq, scales = scales)
gg
}


#utility functions
treatcol.default <- function(n)
	rep(c("black","red","orange","forestgreen",
	"blue","purple","cyan","green"),5)[1:n]

#could be tweaked to better handle multiple response variables (knowing to stack them vertically)
panel.default <- function(x) {
   if(is.null(x)) peq <- NULL
   if(length(x)==1) peq <- as.formula(paste(x,"~."))
   if(length(x)==2) peq <- as.formula(paste(x[1],"~",x[2]))
   if(length(x)>2) peq <- as.formula(paste(paste0(x[-2],collapse="+"), "~", x[2]))
 peq #peq and panel.eq refer to 'panel equation' should probably be panel formula
 }

#does not currently work with "fracdiff" type functions because it is looking for 0-100
pscale <- function(x) {
  out <- x
  out[x>100&!is.na(x)] <- (log10(x[x>100&!is.na(x)])-1)*100
  out[x<(-100)&!is.na(x)] <- -(log10(-x[x<(-100)&!is.na(x)])-1)*100
out}










