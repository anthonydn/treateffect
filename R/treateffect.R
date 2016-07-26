#things to add:
#multiple comp_functions and flexible univariate stats
#think of universal system and solution for NAs. this has been less of a problem recently
#work on the user friendliness of the options like which comp function to use and which comparisons to make
#try to cut down on dependencies by grabbing key functions
#implement pool, block
#next step. work w roman/keats data. chance to do all_comp
#keep kicking tires w multiple data sets
#reporting functions

theme_te <- function() theme_set(theme_bw() %+replace%
  theme(axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  legend.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  plot.background = element_blank() ))

## 1. EXPERIMENTAL DESIGN

treateffect <- function(data, formula, control = NULL,
  times = NULL, subsample = NULL, pool = NULL, block = NULL,
  summary_function = mean.se, comparisons = mcc, comp_function = welchCI,
  extract = NULL, conf.int = 0.95, subset = NULL) {

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

if (!is.null(subsample)) ddf <- ddf %>%
  group_by_(.dots = setdiff(names(ddf), c("response"))) %>%
  summarise(n = n(), response = mean(response, na.rm = T))
#subsample is tricky and requires some more conceptual thinking about how to do it right

subset <- eval(substitute(subset), ddf)
if (!is.null(subset)) ddf <- ddf[subset,]

if (class(ddf[[lpf$right.name]])[1] != "ordered") warning("treatment is not an ordered factor.")

design <- list(response = r, treatment = lpf$right.name,
  times = times, subsample = subsample, pool = pool, block = block)
if (!is.null(lpf$condition)) design$panel <- names(lpf$condition)
design$levels <- levels(lpf$right)
if (is.null(control)) design$control <- levels(lpf$right)[1] else
 design$control <- control
design$basic_formula <- formula(paste("response ~ ",design$treatment))
design$conf.int <- conf.int
design$comp_function <- comp_function
design$extract <- extract
#if (length(r) > 1) design$variable <- "variable"

treatment_summaries <- summary_function(ddf, design) #option for multiple summary functions

comps <- comparisons(ddf, design) #option for multiple comparison functions
design <- comps$design

output <- list(source_data = data, design = design, data = tbl_df(ddf),
  treatment_summaries = treatment_summaries, comparisons = comps$comparisons)
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


## 2. COMPARISON FUNCTIONS
mcc <- function(data, d) {
d$comparisons <- mcc_groups(d) %>% define_comparisons
d$FUN <- function(x, ...) try(d$comp_function$FUN(x, ...), silent = TRUE)
if (is.null(d$extract)) d$extract <- d$comp_function$default_extract else
  d$extract <- d$extract

comparisons <- data %>%
  group_by_(.dots=lapply(c("variable", d$panel, d$times), as.symbol)) %>%
  do(mod = pergroup(., d = d)) %>%
  expand_tbl_matrix(., d = d) %>%
  uldf
list(comparisons = comparisons, design = d)
}
#this is where I can add all_comparisons etc

##PERGROUP FUNCTIONS

pergroup <- function(x, d) {
  m <- lapply(d$comparisons, filter_treatments, data = x, d = d) %>%
  lapply(d$FUN, d = d)
  do.call(sapply, c(list(m), extract_stats, d$extract)) %>%
  t}

filter_treatments <- function(data, d = d, comparisons) data %>%
  filter_(lazyeval::interp(~(t == comparisons$groupA | t == comparisons$groupB),
    t=as.name(d$treatment))) %>% droplevels

extract_stats <- function(x, ...)
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

mean.se <- function(data, d)  {
se <- function(x) sd(x,na.rm=T)/sqrt(length(x))
data %>%
  group_by_(.dots = lapply(c("variable", d$panel, d$times, d$treatment), as.symbol)) %>%
  summarise_(n = ~n(),
    mean = lazyeval::interp(~mean(var, na.rm = T), var = as.name("response")),
    se = lazyeval::interp(~se(var), var = as.name("response")))
}

## DEFINING WHICH COMPARISONS TO DO

define_comparisons <- function(groups) {
  groupA = groups$groupA
  groupB = groups$groupB
  nn <- mapply(function(x, y) as.list(c(groupA = x, groupB = y)),
    groupA, groupB, SIMPLIFY=FALSE)
  names(nn) <- lapply(1:length(groupA), function(x)
    paste(groupB[x], "-", groupA[x])) %>% as.character
  nn
}

mcc_groups <- function(design)
  list(groupA = rep(design$control, length(design$levels)-1),
  groupB = setdiff(design$levels,design$control))


##SPECIFIC FUNCTIONS FOR DOING THE ACTUAL COMPARISONS (t-test vs. bootstrap vs. bayes vs. lme etc)

welchCI = list(FUN = function(data, d = d) tt <- t.test(d$basic_formula, data, conf.level = d$conf.int),
  default_extract = alist(meandiff = estimate[2] - estimate[1],
    CIlo = -conf.int[2], CIhi = -conf.int[1]))

bootfrac <- list(FUN = function(data, d = d) {
  data <- data.frame(data)
  control <- data[data[[d$treatment]] == d$control, "response"] %>% na.omit
  treatment <- data[data[[d$treatment]] != d$control, "response"] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(control,B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(treatment,B=2000, reps=TRUE), "reps")
  list(fracdiff = mean(treatment)/mean(control),
    lo = quantile(trt/con, (1-d$conf.int)/2 ,na.rm=T),
    hi = quantile(trt/con, 1-(1-d$conf.int)/2 ,na.rm=T))
  },
  default_extract = alist(fracdiff = fracdiff, fracdiff_lo = lo, fracdiff_hi = hi))

bootdiff <- list(FUN = function(data, d = d) {
  data <- data.frame(data)
  control <- data[data[[d$treatment]] == d$control, "response"] %>% na.omit
  treatment <- data[data[[d$treatment]] != d$control, "response"] %>% na.omit
  con <- attr(Hmisc::smean.cl.boot(control,B=2000, reps=TRUE), "reps")
  trt <- attr(Hmisc::smean.cl.boot(treatment,B=2000, reps=TRUE), "reps")
  list(bootdiff = mean(treatment)-mean(control),
    lo = quantile(trt-con, (1-d$conf.int)/2 ,na.rm=T),
    hi = quantile(trt-con, 1-(1-d$conf.int)/2 ,na.rm=T))
  },
  default_extract = alist(meandiff = bootdiff, bootdiff_lo = lo, bootdiff_hi = hi))



