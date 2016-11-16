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


## TREATEFFECT FUNCTION

treateffect <- function(data, formula, control = NULL,
  times = NULL, subsample = NULL, pool = NULL, block = NULL,
  summary_functions = c("mean", "se", "CI68"), comp_groups = mcc, comp_function = welchCI,
  extract = NULL, conf.int = 0.95, subset = NULL) {

#formula parsing and data frame construction
lpf <- lattice:::latticeParseFormula(formula, data, multiple = TRUE)
if ("groups" %in% names(lpf)) r <- levels(lpf$groups) else r <- lpf$left.name
ddl <- list()
if (length(r) == 1) ddl$variable <- rep(r, dim(data)[1])
if (length(r)>1) ddl$variable <- lpf$groups
if (!is.null(lpf$condition)) ddl <- c(ddl, lapply(lpf$condition, as.factor))
if (!is.null(times)) ddl[[times]] <- rep(data[[times]], length(r))
if (!is.null(pool)) ddl[[pool]] <- rep(data[[pool]], length(r))
if (!is.null(block)) ddl[[block]] <- rep(data[[block]], length(r))
ddl[[lpf$right.name]] <- lpf$right
ddl$response <- lpf$left
ddf <- lapply(ddl, data.frame) %>% bind_cols
names(ddf) <- names(ddl)

#averaging across subsamples
#(this is tricky and requires some more conceptual thinking about how to do it right)
if (!is.null(subsample)) ddf <- ddf %>%
  group_by_(.dots = setdiff(names(ddf), c("response"))) %>%
  summarise(n = n(), response = mean(response, na.rm = T))

#parsing the subset argument
subset <- eval(substitute(subset), ddf)
if (!is.null(subset)) ddf <- ddf[subset,]

#just a warning about the importance of ordering treatment variable
if (class(ddf[[lpf$right.name]])[1] != "ordered")
  warning("treatment is not an ordered factor.") #this was not triggering when I sent a character vector in as the treatment variable. needs to warn or else you get a cryptic error form the define_comparisons function

#create design list object to keep track of experimental design
d <- list(response = r, treatment = lpf$right.name,
  times = times, subsample = subsample, pool = pool, block = block)
if (!is.null(lpf$condition)) d$panel <- names(lpf$condition)
d$levels <- levels(lpf$right)
if (is.null(control)) d$control <- levels(lpf$right)[1] else
 d$control <- control
d$basic_formula <- formula(paste("response ~ ", d$treatment))
d$conf.int <- conf.int
d$comp_function <- comp_function
d$extract <- extract
#if (length(r) > 1) design$variable <- "variable"

#summaries
treatment_summaries <- ddf %>%
  group_by_(.dots = lapply(c("variable", d$panel, d$times, d$treatment),
    as.symbol)) %>%
  filter(!is.na(response)) %>%
  summarise_all(c("length", summary_functions)) %>%
  rename(n = length) #i believe this is problematic. run it on the "subsampled" labile N data

#comparisons
if (!is.null(comp_groups)) {
  d$comparisons <- comp_groups(d) %>% define_comparisons
  d$FUN <- function(x, ...) try(d$comp_function$FUN(x, ...), silent = TRUE)
  if (is.null(d$extract)) d$extract <- d$comp_function$default_extract else
    d$extract <- d$extract
###########  comparisons <- ddf %>%
    group_by_(.dots = lapply(c("variable", d$panel, d$times), as.symbol)) %>%
    do(mod = pergroup(., d = d)) %>%
    expand_tbl_matrix(., d = d) %>%
    uldf
} else comparisons <- NULL
#option for multiple comparison functions

#output
output <- list(source_data = data, design = d, data = tbl_df(ddf),
  treatment_summaries = treatment_summaries, comparisons = comparisons)
class(output) <- "te"
output
}

print.te <- function(x) {
  cat(paste("Response variable(s):", x$design$response), "\n\n")
  cat("Treatment summaries:", "\n")
  try(print(x$treatment_summaries), silent = TRUE)
  cat("\n\n", "Treatment comparisons:", "\n", sep = '')
  try(print(x$comparisons), silent = TRUE)
}

##PERGROUP FUNCTIONS

pergroup <- function(x, d) {
  m <- lapply(d$comparisons, filter_treatments, data = x, d = d) %>%
    lapply(d$FUN, d = d)
  do.call(sapply, c(list(m), extract_stats, d$extract)) %>%
  t}

filter_treatments <- function(data, d = d, comparisons) data %>%
  filter_(lazyeval::interp(~(t == comparisons$groupA | t == comparisons$groupB),
    t = as.name(d$treatment))) %>% droplevels

extract_stats <- function(x, ...) #broom may help here
  eval(substitute(alist(...))) %>% lapply(function(y)
    tryCatch(with(x, eval(y)), error = function(e) NA))

expand_tbl_matrix <- function(x, d = d) {
  #the name mod comes from the do function
  comps <- do.call(rbind,x$mod)
  comps <- suppressWarnings(data.frame(comparison = ordered(row.names(comps),
    levels = names(d$comparisons)), comps))
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

mcc <- function(design)
  list(groupA = rep(design$control, length(design$levels)-1),
  groupB = setdiff(design$levels,design$control))

allcomps <- function(design)
  list(groupA = combn(design$levels, 2)[1,],
       groupB = combn(design$levels, 2)[2,])

##SPECIFIC FUNCTIONS FOR DOING THE ACTUAL COMPARISONS (t-test vs. bootstrap vs. bayes vs. lme etc)

#http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/

welchCI = list(
  FUN = function(data, d = d) t.test(d$basic_formula, 
    data, conf.level = d$conf.int),
  default_extract = alist(meandiff = estimate[2] - estimate[1],
    CIlo = -conf.int[2], CIhi = -conf.int[1]))

bootfrac <- list(FUN = function(data, d = d) {
  l <- levels(data[[d$treatment]])
  groupA <- data$response[data[[d$treatment]] == l[1]] %>% na.omit
  groupB <- data$response[data[[d$treatment]] == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(fracdiff = mean(groupB)/mean(groupA),
    lo = quantile(trt/con, (1-d$conf.int)/2 ,na.rm=T),
    hi = quantile(trt/con, 1-(1-d$conf.int)/2 ,na.rm=T))},
  default_extract = alist(fracdiff = fracdiff, fracdiff_lo = lo, fracdiff_hi = hi))

bootdiff <- list(FUN = function(data, d = d) {
  l <- levels(data[[d$treatment]])
  groupA <- data$response[data[[d$treatment]] == l[1]] %>% na.omit
  groupB <- data$response[data[[d$treatment]] == l[2]] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(groupA, B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(groupB, B=2000, reps=TRUE), "reps")
  list(bootdiff = mean(groupB) - mean(groupA),
    lo = quantile(trt - con, (1 - d$conf.int) / 2 , na.rm = T),
    hi = quantile(trt - con, 1 - (1 - d$conf.int) / 2 , na.rm = T))},
  default_extract = alist(meandiff = bootdiff, bootdiff_lo = lo, bootdiff_hi = hi))
