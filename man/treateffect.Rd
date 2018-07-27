\name{treateffect}
\alias{treateffect-package}
\alias{treateffect}
\title{Quantify treatment effect sizes in experiments}
\description{
Estimate summary statistics and effect size comparisons for each independent set of replicates in a structured experiment. Choice of functions for summary statistics and comparisons are flexible. Inputs use a formula interface and human-understandable variable types such as responses, treatments, time, and block.
}
\usage{
treateffect(data, formula = NULL,
  response = NULL, treatment = NULL, groups = NULL, times = NULL, block = NULL,
  pool_variance = NULL, average_subsamples = FALSE,
  summary_functions = c("mean", "se", "CI68"), comp_groups = allpairwise,
  control = NULL, comp_function = NULL, conf.int = 0.95,
  CI_derivation = "ML", effect_size_type = "difference")
}
\arguments{
 \item{data}{
a data frame
}
 \item{formula}{
The formula is specified as \code{y1 + y2 + ... ~ x1 + x2 + ... | g1 + g2 + ...} where \code{y1...ym} are response variables \code{1...m} (should be numeric), \code{x1...xi} are treatment variables \code{1...i} (should be categorical). \code{g1...gm} are grouping variables \code{1...g} (factors or character vectors) by which data will be split before analysis, also appearing on separate panels during graphing. This formula interface is similar to the lattice formula interface and in fact uses the lattice formula parser. As in lattice, factor interactions can be specified with the : operator as long as the relevant input vectors are factors (not character). The * operator is not implemented but can be specified as a + b + a:b. On-the-fly transformations and other manipulations can also be used.

Using a double underscore, treatments or groups can be tagged for variance pooling as \code{xi__pool} or \code{gi__pool}. Groups can additionally be tagged as time or blocking variables as \code{gi__time} or \code{gi__block}.
}
 \item{response}{
If a formula is not specified, the response variables y1 + y2 + ... can be specified as a character vector.
}
 \item{treatment}{
If a formula is not specified, the treatment variables x1 + x2 + ... can be specified as a character vector.
}
 \item{groups}{
If a formula is not specified, the response variables g1 + g2 + ... can be specified as a character vector.
}
 \item{times}{
an optional character string indicating the name of a time variable (which can be numeric or of any of the typical date and time formats). Only one time variable can be specified.
}
 \item{block}{
an optional character string indicating the names of variables indicating blocking or pairing structure in the data. This is used by \code{comp_function} functions such as \code{pairedttest}.
}
 \item{pool_variance}{
An optional character vector indicating over which variables to pool the variance. This information can be used by comp_function functions that pool variance such as \code{multcompCI}.
}
 \item{average_subsamples}{
a logical argument indicating whether any remaining replicates should be averaged.
}
 \item{summary_functions}{
a character vector including the name of one or multiple functions to be used to summarize each treatment at each time point within each group. Functions \code{se} and \code{CI68} are included. Others such as \code{SD} or custom functions can be used too.
}
 \item{comp_groups}{
a function that will be used to set up specific comparisons to be made. The three currently available functions are \code{mcc}, \code{allpairwise}, and \code{allothers} specifying multiple comparisons with the control (the control group defined by the control argument), all pairwise comparisons, or comparisons of each treatment vs. all others combined, respectively. The \code{allothers} option currently works only with comp_functions that pool the variance. Viewing the \code{mcc} and \code{allpairwise} functions will show how a set of custom comparisons could be defined. If none is specified \code{comp_function_selector} will select one based on which types of variables have been supplied, e.g., if a blocking variable is defined, it will select a function that appropriately pairs responses.
}
 \item{control}{
an optional character vector the length of the number of specified treatement variables) indicating which treatment level(s) is/are the "control" for each treatment variable if comparisons with the control are desired. If nothing is specified, the first level in each treatment vector is used. Create an ordered variable or specify a control to change which level will be used as the control.
}
 \item{comp_function}{
a function that will be used to make comparisons. See \code{comp_function_selector} for current options (note: this help file does not yet exist!). If none is specified \code{comp_function_selector} will select one for you.
}
 \item{conf.int}{
The alpha level for confidence interval calculations.
CI_derivation = "ML", effect_size_type = "difference"
}
 \item{CI_derivation}{
character vector that will be used by \code{comp_function_selector} to specify how confidence intervals will be calculated. Current options include "ML" for maximum likelihood methods, "bootstrap" for bootstrap methods, and "bayesian" for Bayesian approaches. ML has the most methods implemented and is the default. It allows pooling of variance, ratios or differences, and paired or unpaired designs. However, complex error nestings and such are not currently implmented. Bootstrap approaches can currently do paired/unpaired and ratios or differences, though the current options implemented are not recommended for small samples. The only Bayesian option is simple unpaired t-test for all specified comparisons.
}
 \item{effect_size_type}{
chacter vector specifying how an effect size be reported. The two options are: "difference" (the default) and  "ratio" Keep in mind that ratios often do not make sense for data that contain 0 or negative values.
}
}

\details{
Using a formula interface, treateffect facilitates efficient calculation of the size of treatment effect sizes when comparing multiple categorical treatments. Standard variable types such as response variables (of which there can be multiple for efficient analysis), treatment categories (multiple variables also allowed), a time variable, a blocking variable, variables over which to pool the data, and panel variables over which to divide the data before analysis (e.g. two sites) can be specified. Functions used to calculate the size of the effect are flexible. The specific summary statistics (e.g., mean, SE, SD), comparison function (e.g., confidence interval from Welch t-test or a bootstrapped confidence interval) and the comparisons to perform (e.g., all pairwise comparisons or multiple comparisons with a control) can be specified.
}

\value{
Returns an object of class "\code{te}". A \code{print} method shows the results in tabular form and the \code{plot} and \code{plotdiff}, can be used to plot the results.
 An object of class "\code{te}" is a list containing some or all of the following components:
\item{source_data}{the unmodified data supplied via the \code{data} argument.}
\item{design}{a list of all of the "design" features used to shape the analysis and plotting such as the identities of the treatment, response, time, and block variables}
\item{data}{the data frame used for analysis after restructuring to accommodate for example multiple response variables}
\item{treatment_summaries}{one of two main output data frames showing the output of the \code{summary_functions} functions that were specified for each variable for each treatment group at each time point in each group.}
\item{treatment_comparisons}{the second of two main output data frames showing the output of the analyses by the specified \code{comp_function} for each variable for each comparison specified by the \code{comp_groups} function at each time point in each group.}
}
\author{Anthony Darrouzet-Nardi}

\seealso{
\code{\link{plot.te}}, \code{\link{plotdiff}}, \code{\link{define_comparisons}}
}
\examples{
theme_te() #a more spartan aesthetic

### EXAMPLES WITH SIMULATED DATA

# very basic case: 1 response, 1 treatment with 2 levels
ex1 <- tedatasim(n = 10, response = 1, levels = 2, time = 1,
  groups = 1, subsample = 1, block = FALSE)
ex1.te <- treateffect(ex1, resp_var1 ~ pred_var1)
ex1.te
plot(ex1.te)
plotdiff(ex1.te)

# 2 responses
ex2 <- tedatasim(n = 10, response = 2, levels = 2, time = 1,
  groups = 1, subsample = 1, block = FALSE)
ex2.te <- treateffect(ex2, resp_var1 + resp_var2 ~ pred_var1)
ex2.te
plot(ex2.te)
plotdiff(ex2.te)

# 2 predictors
ex3 <- tedatasim(n = 3, response = 1, levels = c(2,2), time = 1,
  groups = 1, subsample = 1, block = FALSE)
ex3.te <- treateffect(ex3, resp_var1 ~ pred_var1 + pred_var2)
ex3.te
plot(ex3.te)
plotdiff(ex3.te)

# univariate data
ex4.te <- tedatasim(n = 10, response = 1, levels = 1, time = 1,
    groups = 1, subsample = 1, block = FALSE) \%>\%
  treateffect(resp_var1 ~ pred_var1)
ex4.te
plot(ex4.te)

## Examples with time variables

# One treatment over time
ex5 <- tedatasim(n = 5, response = 1, levels = 1, time = 3,
  groups = 1, subsample = 1, block = FALSE)
ex5.te <- treateffect(ex5, resp_var1 ~ pred_var1 | time_var__time)
ex5.te
plot(ex5.te)

# Two treatments over time
ex6 <- tedatasim(n = 5, response = 1, levels = 2, time = 10,
  groups = 1, subsample = 1, block = FALSE)
ex6.te <- treateffect(ex6, resp_var1 ~ pred_var1 | time_var__time)
ex6.te
plot(ex6.te, dodge = 0.3)
plotdiff(ex6.te, dodge = 0.3)

# 3 groups, + time
ex7 <- tedatasim(n = 10, response = 1, levels = 2, time = 30,
  groups = 3, subsample = 1, block = FALSE)
ex7.te <- treateffect(ex7, resp_var1 ~ pred_var1 | group_var1 + time_var__time)
ex7.te
plot(ex7.te, panel_formula = group_var1 ~ .)


## EXAMPLES WITH VARIOUS INTERNAL DATA SETS
#warpbreaks - a commonly used R data set
wb.te <- treateffect(warpbreaks, breaks ~ wool:tension)
wb.te
plot(wb.te)
plotdiff(wb.te)

#barley (using comparisons with all other treatments)
barley.te <- lattice::barley \%>\%
treateffect(yield ~ variety__pool | year, comp_groups = allothers)
barley.te
plot(barley.te)
plotdiff(barley.te)

#city example used in many bootstrapping examples
city_tidy <- boot::city \%>\%
  gather(year, population) \%>\%
  mutate(year = factor(year, labels = c(1920, 1930)), city = rep(1:10, 2))
treateffect(city_tidy, population ~ year + city__block,
  comp_function = bootRR_bca_paired)

#pool variance among treatments with Tukey test
amod <- aov(breaks ~ tension, data = warpbreaks)
confint(multcomp::glht(amod, linfct = mcp(tension = "Tukey")))$confint
treateffect(warpbreaks, breaks ~ tension__pool)

#Dunnett multiple comparisons with control example
amod <- aov(breaks ~ tension, data = warpbreaks)
confint(multcomp::glht(amod, linfct = mcp(tension = "Dunnett")))$confint
treateffect(warpbreaks, breaks ~ tension__pool, comp_groups = mcc)
#specify a different control
treateffect(warpbreaks, breaks ~ tension__pool, comp_groups = mcc,
  control = "H")

#Bayesian BESTmcmc "supersedes the t-test"
y1 <- c(5.77, 5.33, 4.59, 4.33, 3.66, 4.48)
y2 <- c(3.88, 3.55, 3.29, 2.59, 2.33, 3.59)
summary(BESTmcmc(y1,y2))
data.frame(r = c(y1,y2), t = rep(c("y1","y2"), ea = 6)) \%>\%
treateffect(r ~ t, comp_function = BESTHDI)

#starwars dataset from tidyr with lots missing and categories with 1 number
treateffect(starwars, height + mass + birth_year ~ gender) \%>\% plot

#Lots of ways to skin a cat
treateffect(NSE_7mo, NSE_7mo ~ treatment | year) \%>\% plot
treateffect(NSE_7mo, NSE_7mo ~ treatment) \%>\% plotdiff
treateffect(NSE_7mo, NSE_7mo ~ treatment | year) \%>\% plotdiff
treateffect(NSE_7mo, NSE_7mo ~ treatment | year + block__block) \%>\% plotdiff
treateffect(NSE_7mo, NSE_7mo ~ treatment__pool | year + block__block) \%>\% plotdiff
treateffect(NSE_7mo, NSE_7mo ~ treatment__pool | year__pool +
  block__block) \%>\% plotdiff
treateffect(NSE_7mo, NSE_7mo ~ treatment__pool | year__pool) \%>\% plotdiff
treateffect(NSE_7mo, NSE_7mo ~ treatment | year,
  comp_function = bootdiff_bca) \%>\% plotdiff
treateffect(NSE_7mo, NSE_7mo ~ treatment | year__time + block__block,
  comp_groups = mcc) \%>\% plotdiff(dodge = 0.5)
#this next one gives wrong answers and I should find out why bc someone will try it.
treateffect(NSE_7mo, NSE_7mo ~ treatment | year__pool + block__block) \%>\% plotdiff
}
