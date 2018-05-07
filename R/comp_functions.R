##SPECIFIC FUNCTIONS FOR DOING COMPARISONS
#(t-test vs. bootstrap vs. bayes vs. lme etc)
#should output dataframe with effect_size, lwr, upr
#if it pools, it can output multiple rows

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
  tt <- t.test(y ~ x, dfcf, conf.level = d$conf.int, paired = TRUE)
  data.frame(effect_size = -tt$estimate, lwr = -tt$conf.int[2],
    upr = -tt$conf.int[1], row.names = NULL)}

bootdiff_bca <- function(dfcf, d) {
  l <- unique(dfcf$x)
  diff <- function(d, i) {
    E = d[i,] ; mean(E$y[E$x == l[2]]) - mean(E$y[E$x == l[1]])}
  b <- boot::boot(dfcf, diff, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}

bootdiff_bca_paired <- function(dfcf, d) {
  l <- unique(dfcf$x) %>% as.character
  dfcfs <- tidyr::spread(dfcf[,c(d$block, "y", "x")], x, y)
  ratio <- function(d, i) {E = d[i,] ; mean(E[[l[2]]]) - mean(E[[l[1]]])}
  b <- boot::boot(dfcfs, ratio, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}

bootRR_bca <- function(dfcf, d) {
  l <- unique(dfcf$x)
  ratio <- function(d, i) {
    E = d[i,] ; mean(E$y[E$x == l[2]]) / mean(E$y[E$x == l[1]])}
  b <- boot::boot(dfcf, ratio, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}

bootRR_bca_paired <- function(dfcf, d) {
  l <- unique(dfcf$x) %>% as.character
  dfcfs <- tidyr::spread(dfcf[,c(d$block, "y", "x")], x, y)
  ratio <- function(d, i) {E = d[i,] ; mean(E[[l[2]]]) / mean(E[[l[1]]])}
  b <- boot::boot(dfcfs, ratio, R = 1000) %>%
    boot::boot.ci(conf = d$conf.int, type = "bca")
  data.frame(effect_size = b$t0, lwr = b$bca[4], upr = b$bca[5])}

# for this to work, pool_variance has to be set to equal the treatment variable
# and comp_groups set to allcomps
TukeyCI  <- function(dfcf, d) {
  require(multcomp)
  lm(y ~ x, data = dfcf) %>%
  glht(linfct = mcp(x = "Tukey")) %>%
  confint(level = d$conf.int) %>%
  `$`(confint) %>%
  tbl_df}

#not sure you can change 95% to something else in this one
BESTHDI <- function(dfcf, d){
  l <- unique(dfcf$x) %>% as.character
  dfcf_l <- split(dfcf$y, dfcf$x)
  m <- summary(BESTmcmc(dfcf_l[[l[1]]],dfcf_l[[l[2]]]))[3,c(1,5,6)]
  m <- tbl_df(data.frame(as.list(m)))
  names(m) <- c("effect_size", "lwr", "upr")
  m}

