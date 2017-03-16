#things to add:
#multiple comp_functions (this is going to be tricky will require walking through eg)
#make custom comp functions easier
#add in basic nlme functionality, possibly some rms too
#implement pool, block (paired t-tests for one)
#keep kicking tires w multiple data sets
#reporting functions
#try to cut down on dependencies by grabbing key functions
#broom?
#graphical alternative to sig stars - maybe text CIs
#for teaching, learning, tweaking, spit out code for graphs, lme etc
#custom geom for the gradient CI viz
#"main effects" or pooling over crossed treatments should be an option. maybe specified in comp_groups as crossed or something
#bootfrac could be done in a less janky way by figuring out how to use the boot function and potentially something like the BCa method.

##NOT SURE WHAT'S NEXT. POSSIBLY WORK ON POOLING/BLOCKING/SUBSAMPLING

#FORMULA PARSING INTO DATA FRAME AND DESIGN LIST
teParseFormula <- function(data, formula, times = NULL, pool = NULL, block = NULL,
  subsample = NULL, subset = NULL, summary_functions = c("mean", "se", "CI68"),
  comp_groups = mcc, control = NULL, comp_function = welchCI,
  extract = NULL, conf.int = 0.95) {

#create a data frame with a standardized form
lpf <- lattice:::latticeParseFormula(formula, data, multiple = TRUE)
re <- unlist(strsplit(lpf$left.name, ' [+] '))
tr <- unlist(strsplit(lpf$right.name, ' [+] '))
N = length(lpf$left) / length(tr) / length(re)

ddl <- list()
ddl$y_variable <- rep(re, length(tr), ea = N)
ddl$x_variable <- rep(tr, ea = N * length(re))
if (!is.null(lpf$condition)) ddl <- c(ddl, lapply(lpf$condition, as.factor))
if (!is.null(times)) ddl[[times]] <- rep(data[[times]], length(re) * length(tr))
if (!is.null(pool)) ddl[[pool]] <- rep(data[[pool]], length(re) * length(tr))
if (!is.null(block)) ddl[[block]] <- rep(data[[block]], length(re) * length(tr))
ddl$treatment <- if(!is.factor(lpf$right)) factor(lpf$right) else lpf$right
ddl$response <- lpf$left
ddf <- lapply(ddl, data.frame) %>% bind_cols
names(ddf) <- names(ddl)
#averaging across subsamples (not sure this works as intended in all cases)
if (!is.null(subsample)) ddf <- ddf %>%
  group_by_(.dots = setdiff(names(ddf), c("response"))) %>%
  summarise(n = n(), response = mean(response, na.rm = T))
#subsetting
subset <- eval(subset, ddf)
if (!is.null(subset)) ddf <- ddf[subset,]

#create design list object defining the experimental design
d <- list(response = re, treatment = tr, times = times, subsample = subsample,
  pool = pool, block = block)
if (!is.null(lpf$condition)) d$panel <- names(lpf$condition)
d$levels <- lapply(tr, function(x) levels(factor(data[[x]])))
names(d$levels) <- tr
if (is.null(control)) d$control <- lapply(d$levels, `[[`, 1) else
  d$control <- control
d$summary_functions <- summary_functions
d$comp_function <- comp_function
d$extract <- extract
if (!is.null(comp_groups) & !is.null(d$levels)) {
  d$comparisons <-
    lapply(tr, function(x) define_comparisons(comp_groups(d$levels[[x]],
      d$control[[x]]))) %>% unlist(., recursive = FALSE)
  d$FUN <- function(x, ...) try(d$comp_function$FUN(x, ...), silent = TRUE)
  if (is.null(d$extract)) d$extract <- d$comp_function$default_extract else
    d$extract <- d$extract
  d$conf.int <- conf.int
} else comparisons <- NULL

#warning about the importance of ordering treatment variable
if (!is.null(comp_groups) & !is.null(d$levels) &
    class(ddf[[lpf$right.name]])[1] != "ordered")
  warning("treatment is not an ordered factor.")

list(data = tbl_df(ddf), design = d)
}

#TREATMENT SUMMARIES
treatment_summaries <- function(data, design) data %>%
  group_by_(.dots = lapply(c("y_variable", "x_variable", design$panel,
    design$times, "treatment"), as.symbol)) %>%
  filter(!is.na(response)) %>%
  select(-suppressWarnings(one_of(design$block, "n"))) %>% #the problem is that summarise_all does ALL non-grouping columns, which includes block and "n" in the case of a subsample. (may need these at some point we will see)
  summarise_all(c("length", design$summary_functions)) %>% #this doesn't like illegal grouping variable names (ie backticked)
  rename(n = length)

#TREATMENT COMPARISONS
treatment_comparisons <- function(data, design) data %>%
  group_by_(.dots = lapply(c("y_variable", design$panel, design$times), as.symbol)) %>% #note: if you have the same treatment names in multiple x variables it will fuck up
  do(mod = pergroup(., comparisons = design$comparisons, FUN = design$FUN,
    extract = design$extract, conf.int = design$conf.int)) %>%
  expand_tbl_matrix(., comparisons = design$comparisons) %>%
  uldf

## TREATEFFECT FUNCTION (ties the room together does it not?)
treateffect <- function(data, formula, times = NULL, pool = NULL, block = NULL,
  subsample = NULL, subset = NULL, summary_functions = c("mean", "se", "CI68"),
  comp_groups = mcc, control = NULL, comp_function = welchCI,
  extract = NULL, conf.int = 0.95) {

ddf <- teParseFormula(data, formula, times, pool, block, subsample, substitute(subset),
  summary_functions, comp_groups, control, comp_function, extract, conf.int)
summaries <- treatment_summaries(ddf$data, ddf$design)
if (!is.null(ddf$design$comparisons)) comparisons <-
  treatment_comparisons(ddf$data, ddf$design) else compare <- NULL

output <- list(source_data = data, design = ddf$design, data = ddf$data,
  treatment_summaries = summaries, comparisons = comparisons)
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
    try(print(x$comparisons), silent = TRUE)}}




##PERGROUP FUNCTIONS

#comparisons,
pergroup <- function(data, comparisons, FUN, extract, conf.int) {
  m <- lapply(comparisons, filter_treatments, data = data) %>% lapply(FUN, conf.int = conf.int)
  do.call(sapply, c(list(m), extract_stats, extract)) %>%
  t}

filter_treatments <- function(data, comparisons) data %>%
  filter_(lazyeval::interp(~(t == comparisons$groupA | t == comparisons$groupB),
    t = as.name("treatment"))) %>% droplevels

extract_stats <- function(x, ...) #broom may help here
  eval(substitute(alist(...))) %>% lapply(function(y)
    tryCatch(with(x, eval(y)), error = function(e) NA))

expand_tbl_matrix <- function(x, comparisons) {
  #the name mod comes from the do function
  comps <- do.call(rbind, x$mod)
  comps <- suppressWarnings(data.frame(comparison = ordered(row.names(comps),
    levels = names(comparisons)), comps))
  groups <- select(x, -mod) %>% data.frame
  groupsrep <- groups[rep(row.names(groups), ea = dim(x$mod[[1]])[1]),TRUE, drop = FALSE] #causes failure w small data frames, possibly due to one group but maybe not
  cbind(groupsrep, comps)
}
#i have a gut feeling expand_tbl_matrix and uldf could be more elegant
uldf <- function(data) {
  n <- names(data)
  ListCols <- sapply(data, is.list)
  out <- cbind(data[!ListCols], t(apply(data[ListCols], 1, unlist)))
  names(out) <- n
  tbl_df(out)
  }

##UNIVARIATE STATS

se <- function(x, na.rm = FALSE) sd(x, na.rm = na.rm) / sqrt(length(x))

CI <- function(x, p = 0.95, na.rm = FALSE) sd(x, na.rm = na.rm) /
  sqrt(length(x)) * qt(1 - (1 - p) / 2, length(x) - 1)

cv <- function(x, na.rm = FALSE) {sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)}

CI68 <- function(x) CI(x, p = 0.68)

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

allcomps <- function(levels, control)
  list(groupA = combn(levels, 2)[1,],
       groupB = combn(levels, 2)[2,])

##SPECIFIC FUNCTIONS FOR DOING COMPARISONS (t-test vs. bootstrap vs. bayes vs. lme etc)

#http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/

welchCI = list(
  FUN = function(data, conf.int) t.test(response ~ treatment, data, conf.level = conf.int),
  default_extract = alist(meandiff = estimate[2] - estimate[1],
    CIlo = -conf.int[2], CIhi = -conf.int[1]))

bootfrac <- list(FUN = function(data, conf.int) {
  l <- levels(data$treatment)
  groupA <- data$response[data$treatment == l[1]] %>% na.omit
  groupB <- data$response[data$treatment == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(fracdiff = mean(groupB)/mean(groupA),
    lo = quantile(trt/con, (1-conf.int)/2 ,na.rm=T),
    hi = quantile(trt/con, 1-(1-conf.int)/2 ,na.rm=T))},
  default_extract = alist(fracdiff = fracdiff, fracdiff_lo = lo, fracdiff_hi = hi))

bootperc <- list(FUN = function(data, conf.int) {
  l <- levels(data$treatment)
  groupA <- data$response[data$treatment == l[1]] %>% na.omit
  groupB <- data$response[data$treatment == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(percdiff = (mean(groupB) - mean(groupA)) / mean(groupA) * 100,
       lo = quantile((trt - con)/con * 100, (1-conf.int)/2 ,na.rm=T),
       hi = quantile((trt - con)/con * 100, 1-(1-conf.int)/2 ,na.rm=T))},
  default_extract = alist(percdiff = percdiff, percdiff_lo = lo, percdiff_hi = hi))

bootdiff <- list(FUN = function(data, conf.int) {
  l <- levels(data$treatment)
  groupA <- data$response[data$treatment == l[1]] %>% na.omit
  groupB <- data$response[data$treatment == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(bootdiff = mean(groupB) - mean(groupA),
    lo = quantile(trt - con, (1 - conf.int) / 2 , na.rm = T),
    hi = quantile(trt - con, 1 - (1 - conf.int) / 2 , na.rm = T))},
  default_extract = alist(meandiff = bootdiff, bootdiff_lo = lo, bootdiff_hi = hi))
