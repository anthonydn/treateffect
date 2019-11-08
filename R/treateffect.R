## TREATEFFECT FUNCTION
treateffect <- function(data, formula = NULL,
  response = NULL, treatment = NULL, groups = NULL, times = NULL, block = NULL,
  pool_variance = NULL, average_subsamples = FALSE,
  summary_functions = c("mean", "se", "CI68"), comp_groups = allpairwise,
  control = NULL, comp_function = NULL, conf.int = 0.95,
  CI_derivation = "ML", effect_size_type = "difference") {

#create a formula if not specified
if (is.null(formula)) {f <- paste(paste(response, collapse='+'),
  "~", paste(treatment, collapse='+'))
if (!is.null(groups)) f <- paste(f, "|", paste(groups, collapse='+'))
formula <- as.formula(f)}

#parse te formula
f <- teParseFormula(formula)
formula = f$formula
if (is.null(times)) times <- f$times
if (is.null(block)) block <- f$block
if (is.null(pool_variance)) pool_variance <- f$pool_variance

##create an analysis data frame (d_f) with a standardized form
lpf <- lattice:::latticeParseFormula(formula, data, multiple = TRUE)
re <- unlist(strsplit(lpf$left.name, ' [+] ')) #this causes an error because lattice inserts random spaces in a response variable name for long formulas
tr <- unlist(strsplit(lpf$right.name, ' [+] '))
N = length(lpf$left) / length(tr) / length(re)
if (!is.null(lpf$condition) &
  any(names(lpf$condition) %in% c("response", "treatment", "y", "x")))
  names(lpf$condition) <- paste0(names(lpf$condition), "_group")

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
  summarise(subsamples = n(), y = mean(y, na.rm = T)) %>% #this causes a problem because it resports the response variables
  ungroup

#create design list object
d <- list(response = re, treatment = tr, times = times,
  block = block)
if (!is.null(lpf$condition)) d$panel <- names(lpf$condition)
d$levels <- lapply(tr, function(j) levels(factor(ddl$x[ddl$treatment == j])))
names(d$levels) <- tr
if (is.null(control)) d$control <- lapply(d$levels, `[[`, 1) else
  d$control <- as.list(control)
names(d$control) <- tr
d$summary_functions <- lapply(summary_functions, get)
names(d$summary_functions) <- summary_functions
F <- comp_function_selector(comp_function, CI_derivation,
  effect_size_type, block, pool_variance,
  comp_function_name = deparse(substitute(comp_function)))
d$comp_function <- F[[1]]
d$comp_function_name <- F[[2]]
d$contrasts <- deparse(substitute(comp_groups))
if (any(unlist(lapply(d$levels, function(x) length(x) == 1)))) {
  warning("no comparisons done because only one level in treatment")
  comparisons <- NULL
} else if (!is.null(comp_groups) & !is.null(d$levels)) {
  d$comparisons <-
    lapply(tr, function(x) define_comparisons(comp_groups(d$levels[[x]],
      d$control[[x]]))) %>% unlist(., recursive = FALSE)
  #a control does not always need to be specified.
  d$FUN <- function(x, ...) tryCatch(d$comp_function(x, ...), error =
    function(x) data.frame(effect_size = NA, lwr = NA, upr = NA))
  d$conf.int <- conf.int
  d$pool_variance <- pool_variance
} else comparisons <- NULL

#calculations of summaries and comparisons
summaries <- treatment_summaries(d_f, d)

if (!is.null(d$comparisons))
comparisons <- treatment_comparisons(d_f, d) else
comparisons <- NULL

output <- list(source_data = data, design = d, data = d_f,
  summaries = summaries, comparisons = comparisons)
class(output) <- "te"
output
}

#parse te formulas
teParseFormula <- function(formula) {

  s1 <- as.character(formula)[3]
  pool_variance <- str_extract_all(s1, '\\w+__pool') %>%
    unlist %>% str_replace_all('__pool', '')
  s2 <- str_replace_all(s1, '__pool', '')
  if(length(pool_variance) == 0) pool_variance <- NULL
  times <- str_extract_all(s2, '\\w+__time') %>%
    unlist %>% str_replace_all('__time', '')
  s3 <- str_replace_all(s2, ' [\\+\\|] \\w+__time', '')
  if(length(times) == 0) times <- NULL
  block <- str_extract_all(s2, '\\w+__block') %>%
    unlist %>% str_replace_all('__block', '')
  s4 <- str_replace_all(s3, ' [\\+\\|] \\w+__block', '')
  if(length(block) == 0) block <- NULL

  f <- as.formula(paste(as.character(formula)[2], "~", s4))

  list(formula = f, pool_variance = pool_variance, times = times, block = block)
}

#print te objects
print.te <- function(x) {
  cat(paste("Response variable(s):", x$design$response), "\n\n")
  cat("Treatment summaries:", "\n")
  try(print(x$summaries), silent = TRUE)
  if(!is.null(x$design$comparisons)) {
    cat("\n\n", "Treatment comparisons: ", x$design$comp_function_name, "\n", sep = '')
    if (!is.null(x$design$pool_variance)) cat("Variance pooled by: ",
      paste(x$design$pool_variance, collapse = ", "), "\n", sep = '') else
      cat("No variances pooled\n")
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

allpairwise <- function(levels, control)
  list(groupA = combn(levels, 2)[1,],
       groupB = combn(levels, 2)[2,])

allothers <- function(levels, control)
  list(groupA = rep("all others", length(levels)), groupB = levels)

# TREATMENT SUMMARIES

treatment_summaries <- function(d_f, d)
  d_f %>%
  mutate(response = ordered(response, d$response)) %>%
  group_by_at(c("response", "treatment", d$panel, d$times, "x")) %>%
  dplyr::filter(!is.na(y)) %>%
  summarise_at("y", c(n = "length", d$summary_functions)) %>%
  arrange(response)

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
  ea = length(h)), ,drop = FALSE], comparison = ordered(h, levels = h))

#split
fnames <- setdiff(g, d$pool_variance)
ffac <- tidyr::unite(d_f[rev(fnames)], "ffac")[[1]]
ffac <- factor(ffac, levels = unique(ffac))
d_f_split_list <- split(d_f, ffac)

#apply and combine
if (!any(d$treatment %in% d$pool_variance))
  diffs <- lapply(d_f_split_list, pergroup, d) %>% bind_rows else
  diffs <- lapply(d_f_split_list, d$comp_function, d) %>% bind_rows

#sample sizes
h <- lapply(d_f_split_list, function(j) summary(j$x))
iA <- match(unlist(lapply(d$comparisons, function(x) x$groupA)), names(h[[1]]))
iB <- match(unlist(lapply(d$comparisons, function(x) x$groupB)), names(h[[1]]))
cmat$left_n <- unlist(lapply(h, function(x) x[iB]))
cmat$right_n <- unlist(lapply(h, function(x) x[iA]))
if (!is.null(d$pool_variance) & (d$contrasts %in% c("allpairwise", "mcc")))
  cmat <- filter(cmat, left_n != 0, right_n != 0)

#combine comparison matrix with comparison function output
bind_cols(cmat, diffs) %>% tbl_df
}
