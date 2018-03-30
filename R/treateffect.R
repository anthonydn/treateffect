#things to add:

#implement pool, block (paired t-tests for one)
#multiple comp_functions (this is going to be tricky will require walking through eg)
#make custom comp functions easier
#add in basic nlme functionality, possibly some rms too
#keep kicking tires w multiple data sets
#reporting functions
#try to cut down on dependencies by grabbing key functions
#broom?
#graphical alternative to sig stars - maybe text CIs. I think a new density-based visualization is warranted.
#for teaching, learning, tweaking, spit out code for graphs, lme etc
#custom geom for the gradient CI viz
#"main effects" or pooling over crossed treatments should be an option. maybe specified in comp_groups as crossed or something
#bootfrac could be done in a less janky way by figuring out how to use the boot function and potentially something like the BCa method.
#errors can be thrown when there is a bunch of lines of NAs due to a bad import
#if you have the same treatment names in multiple x variables it can cause errors for comparisons
#work on a logical system for internal object names like d_f, data etc


## TREATEFFECT FUNCTION
treateffect <- function(data, formula, times = NULL, pool_variance = NULL,
  block = NULL, replicate_id = NULL, average_subsamples = FALSE,
  summary_functions = c("mean", "se", "CI68"), comp_groups = mcc,
  control = NULL, comp_function = welchCI, conf.int = 0.95) {

##create a data frame with a standardized form
lpf <- lattice:::latticeParseFormula(formula, data, multiple = TRUE)
re <- unlist(strsplit(lpf$left.name, ' [+] '))
tr <- unlist(strsplit(lpf$right.name, ' [+] '))
N = length(lpf$left) / length(tr) / length(re)

ddl <- list(y_variable = rep(re, length(tr), ea = N),
  x_variable = rep(tr, ea = N * length(re)))
if (!is.null(lpf$condition)) ddl <- c(ddl, lapply(lpf$condition, as.factor))
if (!is.null(times)) ddl[[times]] <- rep(data[[times]], length(re) * length(tr))
if (!is.null(block)) ddl[[block]] <- rep(data[[block]], length(re) * length(tr))
if (!is.null(replicate_id)) ddl[[replicate_id]] <- rep(data[[replicate_id]],
  length(re) * length(tr))
ddl$y <- lpf$left
ddl$x <- if(!is.factor(lpf$right)) factor(lpf$right) else lpf$right
d_f <- data.frame(ddl)
#averaging across subsamples
if (!is.null(replicate_id)) average_subsamples = TRUE
if (average_subsamples) d_f <- d_f %>%
  group_by_at(setdiff(names(d_f), c("y"))) %>%
  summarise(subsamples = n(), y = mean(y, na.rm = T))

#create design list object defining the experimental design
d <- list(response = re, treatment = tr, times = times, replicate_id = replicate_id,
  block = block)
if (!is.null(lpf$condition)) d$panel <- names(lpf$condition)
d$levels <- lapply(tr, function(j) levels(factor(ddl$x[ddl$x_variable == j])))
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
  d$FUN <- function(x, ...) try(d$comp_function(x, ...), silent = TRUE)
  d$conf.int <- conf.int
  d$pool_variance <- pool_variance
} else comparisons <- NULL

#warning about the importance of ordering treatment variable
#if (!is.null(comp_groups) & !is.null(d$levels) &
#    class(d_f[["x"]])[1] != "ordered")
#  warning("treatment is not an ordered factor.")

summaries <- treatment_summaries(d_f, d)
if (!is.null(d$comparisons))
compare <- treatment_comparisons(d_f, d) else compare <- NULL

output <- list(source_data = data, design = d, data = d_f,
  treatment_summaries = summaries, treatment_comparisons = compare)
class(output) <- "te"
output
}

#print te objects
print.te <- function(x) {
  cat(paste("Response variable(s):", x$design$response), "\n\n")
  cat("Treatment summaries:", "\n")
  try(print(x$treatment_summaries), silent = TRUE)
  if(!is.null(x$design$comparisons)) {
    cat("\n\n", "Treatment comparisons:", "\n", sep = '')
    try(print(x$treatment_comparisons), silent = TRUE)}}


## DEFINING WHICH COMPARISONS TO DO

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

allcomps <- function(levels, control)#a control does not always need to be specified. awk to have to send control here and not use
  list(groupA = combn(levels, 2)[1,],
       groupB = combn(levels, 2)[2,])


#TREATMENT SUMMARIES

treatment_summaries <- function(data, design)
  data %>%
  group_by_at(c("y_variable", "x_variable", design$panel, design$times, "x")) %>%
  filter(!is.na(y)) %>%
  summarise_at("y", c(n = "length", design$summary_functions)) #this possibly doesn't like illegal grouping variable names (ie backticked)

#univariate stats
se <- function(x, na.rm = FALSE) sd(x, na.rm = na.rm) / sqrt(length(x))
CI <- function(x, p = 0.95, na.rm = FALSE) sd(x, na.rm = na.rm) /
  sqrt(length(x)) * qt(1 - (1 - p) / 2, length(x) - 1)
cv <- function(x, na.rm = FALSE) {sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)}
CI68 <- function(x) CI(x, p = 0.68)


#TREATMENT COMPARISONS

pergroup <- function(dg, d) lapply(d$comparisons, filter_treat, dg = dg) %>%
  lapply(d$FUN, d) %>% bind_rows

filter_treat <- function(compx, dg) dg %>%
  filter(x %in% unlist(compx)) %>% droplevels

treatment_comparisons <- function(d_f, d) {

  #comparison matrix
  g <- c("y_variable", d$panel, d$times)
  a <- lapply(g, function(x) unique(d_f[[x]]))
  names(a) <- g
  a$comparison = names(d$comparisons)
  cmat <- expand.grid(a[length(a):1])[length(a):1]

  #split
  fnames <- setdiff(c("y_variable", d$panel, d$times), d$pool.variance)
  d_f_split_list <- split(d_f, interaction(select(d_f, rev(fnames))))

  #apply and combine
  treatpool <- any(d$treatment %in% d$pool.variance)
  if (!treatpool) diffs <-
    lapply(d_f_split_list, pergroup, d) %>% bind_rows
  if (any(d$treatment %in% d$pool.variance)) diffs <- NULL #do this next

  bind_cols(cmat, diffs) %>% tbl_df
  }


##SPECIFIC FUNCTIONS FOR DOING COMPARISONS (t-test vs. bootstrap vs. bayes vs. lme etc)
#should output dataframe with diff, diffmin, diffmax. if it pools, it can output multiple rows

welchCI <- function(dfcf, d) {
  tt <- t.test(y ~ x, dfcf, conf.level = d$conf.int)
  data.frame(diff = tt$estimate[2] - tt$estimate[1],
    diffmin = -tt$conf.int[2], diffmax = -tt$conf.int[1],
    row.names = NULL)}

pairedttest <- function(dfcf, d) {
  dfcf <- arrange(dfcf, dfcf[[d$block]])
  tt <- t.test(y ~ x, dfcf, conf.level = d$conf.int, paired = TRUE)
  data.frame(diff = -tt$estimate, diffmin = -tt$conf.int[2],
    diffmax = -tt$conf.int[1], row.names = NULL)}


##NOT WORKING

#http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/

bootfrac <- list(FUN = function(data, conf.int) {
  l <- levels(data$x)
  groupA <- data$y[data$x == l[1]] %>% na.omit
  groupB <- data$y[data$x == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(fracdiff = mean(groupB)/mean(groupA),
    lo = quantile(trt/con, (1-conf.int)/2 ,na.rm=T),
    hi = quantile(trt/con, 1-(1-conf.int)/2 ,na.rm=T))},
  default_extract = alist(fracdiff = fracdiff, fracdiff_lo = lo, fracdiff_hi = hi))

bootperc <- list(FUN = function(data, conf.int) {
  l <- levels(data$x)
  groupA <- data$y[data$x == l[1]] %>% na.omit
  groupB <- data$y[data$x == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(percdiff = (mean(groupB) - mean(groupA)) / mean(groupA) * 100,
       lo = quantile((trt - con)/con * 100, (1-conf.int)/2 ,na.rm=T),
       hi = quantile((trt - con)/con * 100, 1-(1-conf.int)/2 ,na.rm=T))},
  default_extract = alist(percdiff = percdiff, percdiff_lo = lo, percdiff_hi = hi))

bootdiff <- list(FUN = function(data, conf.int) {
  l <- levels(data$x)
  groupA <- data$y[data$x == l[1]] %>% na.omit
  groupB <- data$y[data$x == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(bootdiff = mean(groupB) - mean(groupA),
    lo = quantile(trt - con, (1 - conf.int) / 2 , na.rm = T),
    hi = quantile(trt - con, 1 - (1 - conf.int) / 2 , na.rm = T))},
  default_extract = alist(meandiff = bootdiff, bootdiff_lo = lo, bootdiff_hi = hi))
