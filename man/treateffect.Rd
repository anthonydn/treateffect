\name{treateffect}
\alias{treateffect-package}
\alias{treateffect}
\title{Quantify treatment effect sizes in structured experiments}
\description{
Estimate summary statistics and effect size comparisons for each independent set of replicates in a structured experiment. Choice of functions for summary statistics and comparisons are flexible. Inputs use a formula interface and human-understandable variable types such as responses, treatments, time, and block.
}
\usage{
treateffect(data, formula, times = NULL, pool_variance = NULL,
  block = NULL, replicate_id = NULL, average_subsamples = FALSE,
  summary_functions = c("mean", "se", "CI68"), comp_groups = mcc,
  control = NULL, comp_function = welchCI, conf.int = 0.95)
}
\arguments{
 \item{data}{
a data frame
}
 \item{formula}{
The formula is specified as \code{y1 + y2 + ... ~ x1 + x2 + ... | g1 + g2 + ...} where \code{y1...ym} are response variables \code{1...m} (should be numeric), \code{x1...xi} are treatment variables \code{1...i} (should be categorical). \code{g1...gm} are grouping variables \code{1...g} (factors or character vectors) by which data will be split before analysis, also appearing on separate panels during graphing. This formula interface is similar to the lattice formula interface and in fact uses the lattice formula parser. On-the-fly transformations and other manipulations can be used.
}
 \item{times}{
an optional character string indicating the name of a time variable (which can be numeric or of any of the typical date and time formats). Only one time variable can be specified.
}
 \item{pool_variance}{
An optional character vector indicating over which variables to pool the variance. This information can be used by  comp_function functions that pool variance (none are yet implemented, but a comp_function based on lme for example would use this).
}
 \item{block}{
an optional character string indicating the names of variables indicating blocking or pairing structure in the sample. This is used by comp_function functions such as pairedttest.
}
 \item{replicate_id}{
an optional character string indicating the names of variables over which subsamples should be averaged. The \code{mean()} function is used for averaging.
}
 \item{average_subsamples}{
a logical argument indicating whether any remaining replicates should be averaged. This can be specified in the case that subsamples exist but no replicate_id variable is specified.
}
 \item{summary_functions}{
a character vector including the name of one or multiple functions to be used to summarize each treatment at each time point within each group. Functions \code{se} and \code{CI68} are included. Others such as \code{SD} or custom functions can be used too.
}
 \item{comp_groups}{
a function that will be used to set up specific comparisons to be made. The two currently available functions are \code{mcc} and \code{allcomps}, specifying multiple comparisons with the control (the control group defined by the control argument) or all pairwise comparisons, respectively. Viewing the \code{mcc} and \code{allcomps} functions will show how a set of custom comparisons could be defined.
}
 \item{control}{
an optional character string indicating which treatment level(s) is/are the "control" for each treatment variable if comparisons with the control are desired. Comparisons of each treatment using the \code{mcc()} function are the default. If nothing is specified, the first level in the treatment vector (preferably an ordered factor) is used.
}
 \item{comp_function}{
a function that will be used to make comparisons. Currently available options include \code{welchCI} for confidence intervals generated with a Welch t-test around the subtractive difference treatment - control; \code{bootdiff} for a bootstrapped confidence interval around the subtractive difference treatment - control; \code{bootfrac} for the treatment/control quotient with a bootstrapped confidence interval; and \code{bootperc} for a "percent difference," (treatment - control) / control * 100 with a bootstrap confidence interval. Custom functions can be defined and more are planned to be added. Ultimately, more complex analyses including approaches such as mixed models could be added.
}
 \item{conf.int}{
The alpha level for confidence interval calculations.
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
    groups = 1, subsample = 1, block = FALSE) %>%
  treateffect(resp_var1 ~ pred_var1)
ex4.te
plot(ex4.te)

## Examples with time variables

# One treatment over time
ex5 <- tedatasim(n = 5, response = 1, levels = 1, time = 3,
  groups = 1, subsample = 1, block = FALSE)
ex5.te <- treateffect(ex5, resp_var1 ~ pred_var1, times = "time_var")
ex5.te
plot(ex5.te)

# Two treatments over time
ex6 <- tedatasim(n = 5, response = 1, levels = 2, time = 10,
  groups = 1, subsample = 1, block = FALSE)
ex6.te <- treateffect(ex6, resp_var1 ~ pred_var1, times = "time_var")
ex6.te
plot(ex6.te, dodge = 0.3)
plotdiff(ex6.te, dodge = 0.3)

# 3 groups, + time
ex7 <- tedatasim(n = 10, response = 1, levels = 2, time = 30,
  groups = 3, subsample = 1, block = FALSE)
ex7.te <- treateffect(ex7, resp_var1 ~ pred_var1 | group_var1, times = "time_var")
ex7.te
plot(ex7.te, panel_formula = group_var1 ~ .)


## EXAMPLES WITH CLASSIC DATA SETS
wb.te <- treateffect(warpbreaks, breaks ~
  interaction(wool, tension), comp_groups = allcomps)
wb.te
plot(wb.te)
plotdiff(wb.te)

barley.te <- treateffect(lattice::barley, yield ~ variety | year)
barley.te
plot(barley.te)
plotdiff(barley.te)

city_tidy <- boot::city \%\>%
  gather(year, population) \%\>%
  mutate(year = factor(year, labels = c(1920, 1930)), city = rep(1:10, 2))
treateffect(city_tidy, population ~ year, block = "city",
  comp_function = bootfrac_bca_paired)
}
