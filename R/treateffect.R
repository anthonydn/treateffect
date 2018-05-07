#things to add:

#expand pool of comp_functions, especially for pooled variance analyses
  #dunnet
  #bootperc
  #bootdiff
#try to cut down on dependencies by grabbing key functions
#if you have the same treatment names in multiple x variables it can cause errors for comparisons. add a warning
#work on a logical system for internal object names like d_f, data etc
#help page for comp_functions
#help page for plot.te and plotdiff
#help page for simulate_te_data
#add a bunch of the datasets the students looked at(go through carefully and make sure they all make sense)
#it's still worth adding examples of my data which tend to be messier and more complex than included data.
  #would also be cool to put them in the package if they fit
#go through that script with the Harrell example and see if it could somehow me included

## TREATEFFECT FUNCTION
treateffect <- function(data, formula, times = NULL, block = NULL,
  pool_variance = NULL, average_subsamples = FALSE,
  summary_functions = c("mean", "se", "CI68"), comp_groups = allcomps,
  control = NULL, comp_function = welchCI, conf.int = 0.95) {

##create a data frame with a standardized form
lpf <- lattice:::latticeParseFormula(formula, data, multiple = TRUE)
re <- unlist(strsplit(lpf$left.name, ' [+] '))
tr <- unlist(strsplit(lpf$right.name, ' [+] '))
N = length(lpf$left) / length(tr) / length(re)

ddl <- list(response = rep(re, length(tr), ea = N),
  treatment = rep(tr, ea = N * length(re)))
if (!is.null(lpf$condition)) ddl <- c(ddl, lapply(lpf$condition, as.factor))
if (!is.null(times)) ddl[[times]] <- rep(data[[times]], length(re) * length(tr))
if (!is.null(block)) ddl[[block]] <- rep(data[[block]], length(re) * length(tr))
ddl$y <- lpf$left
ddl$x <- if(!is.factor(lpf$right)) factor(lpf$right) else lpf$right
d_f <- data.frame(ddl)
#averaging across subsamples
if (average_subsamples) d_f <- d_f %>%
  group_by_at(setdiff(names(d_f), c("y"))) %>%
  summarise(subsamples = n(), y = mean(y, na.rm = T)) %>%
  ungroup

#create design list object defining the experimental design
d <- list(response = re, treatment = tr, times = times,
  block = block)
if (!is.null(lpf$condition)) d$panel <- names(lpf$condition)
d$levels <- lapply(tr, function(j) levels(factor(ddl$x[ddl$treatment == j])))
names(d$levels) <- tr
if (is.null(control)) d$control <- lapply(d$levels, `[[`, 1) else
  d$control <- control
names(d$control) <- tr
d$summary_functions <- summary_functions
d$comp_function <- comp_function
d$comp_function_name <- deparse(substitute(comp_function))
if (any(unlist(lapply(d$levels, function(x) length(x) == 1)))) {
  warning("no comparisons done because only one level in treatment")
  comparisons <- NULL
} else if (!is.null(comp_groups) & !is.null(d$levels)) {
  d$comparisons <-
    lapply(tr, function(x) define_comparisons(comp_groups(d$levels[[x]],
      d$control[[x]]))) %>% unlist(., recursive = FALSE) #a control does not always need to be specified.
  d$FUN <- function(x, ...) tryCatch(d$comp_function(x, ...), error =
    function(x) data.frame(effect_size = NA, lwr = NA, upr = NA))
  d$conf.int <- conf.int
  d$pool_variance <- pool_variance
} else comparisons <- NULL

#warning about the importance of ordering treatment variable. only if mcc and other conditions
#if (!is.null(comp_groups) & !is.null(d$levels) &
#    class(d_f[["x"]])[1] != "ordered")
#  warning("treatment is not an ordered factor.")

summaries <- treatment_summaries(d_f, d)
if (!is.null(d$comparisons))
comparisons <- treatment_comparisons(d_f, d) else comparisons <- NULL

output <- list(source_data = data, design = d, data = d_f,
  summaries = summaries, comparisons = comparisons)
class(output) <- "te"
output
}

#print te objects
print.te <- function(x) {
  cat(paste("Response variable(s):", x$design$response), "\n\n")
  cat("Treatment summaries:", "\n")
  try(print(x$summaries), silent = TRUE)
  if(!is.null(x$design$comparisons)) {
    cat("\n\n", "Treatment comparisons: ", x$design$comp_function_name, "\n", sep = '')
    try(print(x$comparisons), silent = TRUE)}}


# DEFINING WHICH COMPARISONS TO DO

define_comparisons <- function(groups) {
  groupA = groups$groupA
  groupB = groups$groupB
  nn <- mapply(function(x, y) as.list(c(groupA = x, groupB = y)),
    groupA, groupB, SIMPLIFY=FALSE)
  names(nn) <- lapply(1:length(groupA), function(x)
    paste(groupB[x], "-", groupA[x])) %>% as.character
  nn}

mcc <- function(levels, control)
  list(groupA = rep(control, length(levels)-1),
       groupB = setdiff(levels, control))

allcomps <- function(levels, control) #a control does not always need to be specified. awk to have to send control here and not use
  list(groupA = combn(levels, 2)[1,],
       groupB = combn(levels, 2)[2,])


# TREATMENT SUMMARIES

treatment_summaries <- function(d_f, d)
  d_f %>%
  mutate(response = ordered(response, d$response)) %>%
  group_by_at(c("response", "treatment", d$panel, d$times, "x")) %>%
  dplyr::filter(!is.na(y)) %>%
  summarise_at("y", c(n = "length", d$summary_functions)) %>%  #this possibly doesn't like illegal grouping variable names (ie backticked)
  arrange(response) #should probably implement for grouping vars and treatments (multiple x) too

#univariate stats
se <- function(x, na.rm = FALSE) sd(x, na.rm = na.rm) / sqrt(length(x))
CI <- function(x, p = 0.95, na.rm = FALSE) sd(x, na.rm = na.rm) /
  sqrt(length(x)) * qt(1 - (1 - p) / 2, length(x) - 1)
cv <- function(x, na.rm = FALSE) {sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)}
CI68 <- function(x) CI(x, p = 0.68)


# TREATMENT COMPARISONS

pergroup <- function(dg, d) lapply(d$comparisons, filter_treat, dg = dg) %>%
  lapply(d$FUN, d) %>% bind_rows

filter_treat <- function(compx, dg) dg %>%
  dplyr::filter(x %in% unlist(compx)) %>% droplevels

treatment_comparisons <- function(d_f, d) {

#comparison matrix
g <- c("response", d$panel, d$times)
cmat <- unique(d_f[g])
h <- names(d$comparisons)
cmat <- data.frame(cmat[rep(row.names(cmat),
  ea = length(h)), ,drop = FALSE], comparison = h)

#split
fnames <- setdiff(g, d$pool_variance)
ffac <- tidyr::unite(dplyr::select(d_f, rev(fnames)))[[1]]
ffac <- factor(ffac, levels = unique(ffac))
d_f_split_list <- split(d_f, ffac)

#apply and combine
treatpool <- any(d$treatment %in% d$pool.variance)
if (!treatpool) diffs <- lapply(d_f_split_list, pergroup, d) %>% bind_rows
if (any(d$treatment %in% d$pool_variance)) diffs <-
  lapply(d_f_split_list, d$comp_function, d) %>% bind_rows

#combine comparison matrix with comparison function output
bind_cols(cmat, diffs) %>% tbl_df
}
