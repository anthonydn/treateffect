##SPECIFIC FUNCTIONS FOR DOING COMPARISONS

comp_function_selector <- function(comp_function, CI_derivation,
  effect_size_type, block, pool_variance, comp_function_name) {
	x = NULL
  if(CI_derivation == "ML" & effect_size_type == "difference" &
    is.null(block) & is.null(pool_variance))
	  x = list(welchCI, "welchCI")
  if(CI_derivation == "ML" & effect_size_type == "difference" &
    !is.null(block) & is.null(pool_variance))
	  x = list(pairedttest, "pairedttest")
  if(CI_derivation == "ML" & effect_size_type == "difference" &
    is.null(block) & !is.null(pool_variance))
	  x = list(mult_compCI, "mult_compCI")
  if(CI_derivation == "ML" & effect_size_type == "difference" &
    !is.null(block) & !is.null(pool_variance))
	  x = list(mixed_effectsCI, "mixed_effectsCI")

  if(CI_derivation == "ML" & effect_size_type == "ratio" &
    is.null(block) & is.null(pool_variance))
	  x = list(mratioRR, "mratioRR")
  if(CI_derivation == "ML" & effect_size_type == "ratio" &
    !is.null(block) & is.null(pool_variance))
	  x = list(pairedRR, "pairedRR")

  if(CI_derivation == "bootstrap" & effect_size_type == "difference" &
    is.null(block) & is.null(pool_variance))
	  x = list(bootdiff_bca, "bootdiff_bca")
  if(CI_derivation == "bootstrap" & effect_size_type == "difference" &
    !is.null(block) & is.null(pool_variance))
	  x = list(bootdiff_bca_paired, "bootdiff_bca_paired")
  if(CI_derivation == "bootstrap" & effect_size_type == "ratio" &
    is.null(block) & is.null(pool_variance))
	  x = list(bootRR_bca, "bootRR_bca")
  if(CI_derivation == "bootstrap" & effect_size_type == "ratio" &
    !is.null(block) & is.null(pool_variance))
	  x = list(bootRR_bca_paired, "bootRR_bca_paired")

  if(CI_derivation == "bayesian" & effect_size_type == "difference" &
    is.null(block) & is.null(pool_variance))
	  x = list(BESTHDI, "BESTHDI")

  if(!is.null(comp_function))
    x = list(comp_function, comp_function_name)

  if(!is.null(x)) x else stop("Specified comparison approach not implemented")}


#Maximum likelihood (ML) differences including for pooled variances
welchCI <- function(dfcf, d) {
  tt <- tryCatch(t.test(y ~ x, dfcf, conf.level = d$conf.int),
    error = function(x) list(conf.int = c(NA, NA)))
  data.frame(effect_size = diff(aggregate(y ~ x, dfcf, mean)$y),
    lwr = -tt$conf.int[2], upr = -tt$conf.int[1],
    row.names = NULL)}

twosamplettest <- function(dfcf, d) {
  tt <- t.test(y ~ x, dfcf, conf.level = d$conf.int, var.equal = TRUE)
  data.frame(effect_size = tt$estimate[2] - tt$estimate[1],
    lwr = -tt$conf.int[2], upr = -tt$conf.int[1],
    row.names = NULL)}

pairedttest <- function(dfcf, d) {
  dfcf <- arrange(dfcf, dfcf[[d$block]])
  tt <- t.test(y ~ x, dfcf, conf.level = d$conf.int, paired = TRUE, na.action = na.omit)
  data.frame(effect_size = -tt$estimate, lwr = -tt$conf.int[2],
    upr = -tt$conf.int[1], row.names = NULL)}

mult_compCI  <- function(dfcf, d) {
require("emmeans")
gvars <- setdiff(d$pool_variance, d$treatment)
if (length(gvars) > 0 & !is.null(d$times))
  gvars[gvars == d$times] <- paste("factor(", d$times, ")")
fixed <- as.formula(paste("y ~", paste(c("x", gvars), collapse = "*")))
if (d$contrasts == "allpairwise") co <- "revpairwise" else
if (d$contrasts == "mcc") co <- "trt.vs.ctrl"
if (d$contrasts == "allothers") co <- "del.eff"
if (length(gvars) > 0) contrasts <- as.formula(paste(co, "~x|",
  paste(gvars, collapse = "+"))) else
  contrasts <- as.formula(paste(co, "~x"))
lm(fixed, data = dfcf) %>%
  emmeans(contrasts) %>%
  confint %>%
  `$`(contrasts) %>%
  as_tibble %>%
  transmute(effect_size = estimate, lwr = lower.CL, upr = upper.CL)}

mixed_effectsCI  <- function(dfcf, d) {
require("nlme")
require("emmeans")
gvars <- setdiff(d$pool_variance, d$treatment)
if (length(gvars) > 0 & !is.null(d$times))
  gvars[gvars == d$times] <- paste("factor(", d$times, ")")
fixed <- as.formula(paste("y ~", paste(c("x", gvars), collapse = "*")))
random <- as.formula(paste("~ 1 |", d$block))
if (d$contrasts == "allpairwise") co <- "revpairwise" else
if (d$contrasts == "mcc") co <- "trt.vs.ctrl"
if (d$contrasts == "allothers") co <- "del.eff"
if (length(gvars) > 0) contrasts <- as.formula(paste(co, "~x|",
  paste(gvars, collapse = "+"))) else
  contrasts <- as.formula(paste(co, "~x"))
lme(fixed, random, data = dfcf) %>%
  emmeans(contrasts) %>%
  confint %>%
  `$`(contrasts) %>%
  as_tibble %>%
  transmute(effect_size = estimate, lwr = lower.CL, upr = upper.CL)}


# ML ratios
mratioRR <- function(dfcf, d) {
  require("mratios")
  tt <- ttestratio(y ~ x, data = dfcf, var.equal = FALSE,
    conf.int = d$conf.int, base = 1)
  data.frame(effect_size = tt$estimate[3],
    lwr = tt$conf.int[1],
    upr = tt$conf.int[2]) }

pairedRR <- function(dfcf, d) {
  dfcf <- arrange(dfcf, dfcf[[d$block]])
  tt <- t.test(-log(y) ~ x, dfcf, conf.level = d$conf.int, paired = TRUE, na.action = na.omit)
  data.frame(effect_size = exp(tt$estimate), lwr = exp(tt$conf.int[1]),
    upr = exp(tt$conf.int[2]), row.names = NULL)}

## Ratios for zero-inflated unpaired log-normal data
#note: it doesn't do NA
zhouRR <- function(dfcf,d) {
  l <- unique(dfcf$x) %>% as.character
  dfcf_l <- split(dfcf$y, dfcf$x)
  W1 <- dfcf_l[[1]]
  W2 <- dfcf_l[[2]]
  n1 <- length(W1)
  n2 <- length(W2)
  yi1 <- log(W1[W1 != 0])
  yi2 <- log(W2[W2 != 0])
  n11 <- length(yi1)
  n12 <- length(yi2)
  d1 = sum(W1 == 0)/length(W1)
  d2 = sum(W2 == 0)/length(W2)
  u1 = mean(mean(yi1))
  u2 = mean(mean(yi2))

  s1 <- sqrt((1/(n11-1))*sum((yi1 - u1)^2))
  s2 <- sqrt((1/(n12-1))*sum((yi2 - u2)^2))

  M1 <- (1-d1)*exp(u1 + (s1^2/2))
  M2 <- (1-d2)*exp(u2 + (s2^2/2))

  se2 <- d1/(n1*(1-d1)) + s1^2/(n1*(1-d1)) + s1^4/(2*n1*(1-d1)) +
    d2/(n2*(1-d2)) + s2^2/(n2*(1-d2)) + s2^4/(2*n2*(1-d2))

  data.frame(effect_size = M2/M1, #mean(W2)/mean(W1),
    lwr = exp(log(M2) - log(M1) - qnorm(.975) * sqrt(se2)),
    upr = exp(log(M2) - log(M1) + qnorm(.975) * sqrt(se2)))}


# Bootstrapping (differences and ratios)
bootdiff_bca <- function(dfcf, d) {
  dfcf <- na.omit(dfcf)
  l <- levels(dfcf$x)
  diff <- function(d, i) {
    E = d[i,] ; mean(E$y[E$x == l[2]]) - mean(E$y[E$x == l[1]])}
  b <- boot::boot(dfcf, diff, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}

bootdiff_bca_paired <- function(dfcf, d) {
  dfcf <- na.omit(dfcf)
  l <- levels(dfcf$x) %>% as.character
  dfcfs <- tidyr::spread(dfcf[,c(d$block, "y", "x")], x, y)
  ratio <- function(d, i) {E = d[i,] ; mean(E[[l[2]]]) - mean(E[[l[1]]])}
  b <- boot::boot(dfcfs, ratio, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}

bootRR_bca <- function(dfcf, d) {
  dfcf <- na.omit(dfcf)
  l <- levels(dfcf$x)
  ratio <- function(d, i) {
    E = d[i,] ; mean(E$y[E$x == l[2]]) / mean(E$y[E$x == l[1]])}
  b <- boot::boot(dfcf, ratio, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}

bootRR_bca_paired <- function(dfcf, d) {
  dfcf <- na.omit(dfcf)
  l <- levels(dfcf$x) %>% as.character
  dfcfs <- tidyr::spread(dfcf[,c(d$block, "y", "x")], x, y)
  ratio <- function(d, i) {E = d[i,] ; mean(E[[l[2]]]) / mean(E[[l[1]]])}
  b <- boot::boot(dfcfs, ratio, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}


#Bayesian t-test
#not sure you can change 95% to something else in this one
BESTHDI <- function(dfcf, d){
  require(BEST)
  l <- unique(dfcf$x) %>% as.character
  dfcf_l <- split(dfcf$y, dfcf$x)
  m <- summary(BESTmcmc(dfcf_l[[l[1]]],dfcf_l[[l[2]]]))[3,c(1,5,6)]
  m <- as_tibble(data.frame(as.list(m)))
  names(m) <- c("effect_size", "lwr", "upr")
  m}
