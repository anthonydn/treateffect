#things to add:
#multiple response variables
#multiple comp_functions and flexible univariate stats
#transformations in formula
#pass conf.int
#think of universal system and solution for NAs
#work on the user friendliness of the options like which comp function to use and whcih comparisons to make
#try to cut down on dependencies by grabbing key functions

##UNIVARIATE STATS

mean.se.boot <- function(data)  { #probably break these up a bit
d <- attributes(data)$design

se <- function(x) sd(x,na.rm=T)/sqrt(length(x))
dat.sum <- data %>%
  group_by_(.dots = lapply(c(d$panel, d$times, d$treatment), as.symbol)) %>%
  summarise_(n = ~n(),
  mean = interp(~mean(var, na.rm = T), var = as.name(d$response)),
  se = interp(~se(var), var = as.name(d$response)))

bci <- data %>%
  group_by_(.dots = lapply(c(d$panel, d$times, d$treatment), as.symbol)) %>%
  do(mod = smean.cl.boot(.[,d$response][[1]], conf.int = 0.68)) %>%
  cbind(.,
  summarise(.,
    bootselo = nth(mod, 2),
    bootsehi = nth(mod, 3))) %>%
  select(-mod)

univariate <- merge(dat.sum,bci)
output <- list(data = data, univariate = univariate)
class(output) <- "te"
output
}

##COMPARISON FUNCTIONS
mcc <- function(data, comp_function = welchCI, extract = NULL, conf.int = 0.95) {

d <- attributes(data)$design
d$conf.int <- conf.int
d$comparisons <- mcc_groups(d) %>% define_comparisons
d$FUN <- function(x, ...) try(comp_function$FUN(x, ...), silent = TRUE)
if (is.null(extract)) d$extract <- comp_function$default_extract else
  d$extract <- extract
attributes(data)$design <- d

analyses <- data %>%
  group_by_(.dots=lapply(c(d$panel,
    d$times),as.symbol)) %>%
  do(mod = pergroup(., d = d)) %>%
  expand_tbl_matrix %>%
  uldf

output <- list(data = data, analyses = analyses)
class(output) <- "te"
output
}

#this is where I can add all_comparisons etc

#print function
print.te <- function(x) {
  try(print(x$univariate), silent = TRUE)
  try(print(x$analyses), silent = TRUE)
  }

##EXPERIMENTAL DESIGN

##this function should warn you if you do not have an ordered factor for the treatment variable
exp_design <- function(data, formula, control = NULL, block = NULL,
  times = NULL) {
#implement model.frame and transformations here.
#also needs to deal with multiple response variables
lhs <- all.vars(formula[[2]])
rhs <- all.vars(formula[[3]])
design <- list(response = lhs, treatment = rhs[1], block = block, times = times)
if (length(rhs)>1) design$panel <- rhs[2:length(rhs)]
design$levels <- levels(data[,design$treatment])
if (is.null(control)) design$control <- levels(data[,design$treatment])[1] else
 design$control <- control
design$basic_formula <- formula(paste(design$response,"~",design$treatment))

attributes(data)$design <- design
class(data) <- c(class(data), "design")
data
}

##GENERAL PERGROUP FUNCTIONS

pergroup <- function(x, d) {
  m <- lapply(d$comparisons, filter_treatments, data = x, d = d) %>%
  lapply(d$FUN, d = d)
  do.call(sapply, c(list(m), extract_stats, d$extract)) %>%
  t}

filter_treatments <- function(data, d = d, comparisons) data %>%
  filter_(interp(~(t == comparisons$groupA | t == comparisons$groupB),
    t=as.name(d$treatment))) %>% droplevels

extract_stats <- function(x, ...)
  eval(substitute(alist(...))) %>% lapply(function(y)
    tryCatch(with(x, eval(y)), error = function(e) NA))

##CLEAN OUPUT

expand_tbl_matrix <- function(x) {
  #the name mod comes from the do function below. it could be made generic
  comps <- do.call(rbind,x$mod)
  comps <- suppressWarnings(data.frame(comparison = ordered(row.names(comps),
    levels = names(attributes(x)$design$comparisons)), comps))
  groups <- select(x, -mod) %>% data.frame
  groupsrep <- groups[rep(row.names(groups), ea = dim(x$mod[[1]])[1]),]
  cbind(groupsrep, comps)
}

uldf <- function(data) {
  n <- names(data)
  ListCols <- sapply(data, is.list)
  out <- cbind(data[!ListCols], t(apply(data[ListCols], 1, unlist)))
  names(out) <- n
  out
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
welchCI = list(FUN = function(data, d = d) tt <- t.test(d$basic_formula, data),
  default_extract = alist(meandiff = estimate[2] - estimate[1],
    CIlo = -conf.int[2], CIhi = -conf.int[1]))

bootfrac <- list(FUN = function(data, d = d) {
  data <- data.frame(data)
  control <- data[data[,d$treatment] == d$control, d$response] %>% na.omit
  treatment <- data[data[,d$treatment] != d$control, d$response] %>% na.omit
  con <- attr(smean.cl.boot(control,B=2000, reps=TRUE), "reps")
  trt <- attr(smean.cl.boot(treatment,B=2000, reps=TRUE), "reps")
  list(fracdiff = mean(treatment)/mean(control),
    lo = quantile(trt/con, .025 ,na.rm=T),
    hi = quantile(trt/con, .975 ,na.rm=T))
  },
  default_extract = alist(fracdiff = fracdiff, fracdiff_lo = lo, fracdiff_hi = hi))

bootdiff <- list(FUN = function(data, d = d) {
  data <- data.frame(data)
  control <- data[data[,d$treatment] == d$control, d$response] %>% na.omit
  treatment <- data[data[,d$treatment] != d$control, d$response] %>% na.omit
  con <- attr(smean.cl.boot(control,B=2000, reps=TRUE), "reps")
  trt <- attr(smean.cl.boot(treatment,B=2000, reps=TRUE), "reps")
  list(bootdiff = mean(treatment)-mean(control),
    lo = quantile(trt-con, .025 ,na.rm=T),
    hi = quantile(trt-con, .975 ,na.rm=T))
  },
  default_extract = alist(meandiff = bootdiff, bootdiff_lo = lo, bootdiff_hi = hi))



