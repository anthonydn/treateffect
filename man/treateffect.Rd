\name{treateffect}
\alias{treateffect-package}
\alias{treateffect}
\title{Quantify treatment effect sizes in structured experiments}
\description{
Organize data and perform a "many t-tests" type analysis in which effect sizes are independently calculated for each specficic contrast in lieu of a complex omnibus model. Inputs use a formula interface and human-understandable data variable types such as responses, treatments, time, and block.
}
\usage{
treateffect(data, formula, times = NULL, pool.variance = NULL, block = NULL,
  replicate_id = NULL, average_subsamples = FALSE, summary_functions = c("mean", "se", "CI68"),
  comp_groups = mcc, control = NULL, comp_function = welchCI,
  extract = NULL, conf.int = 0.95)
}
\arguments{
 \item{data}{
a data frame or tibble
}
 \item{formula}{
The formula is specified as \code{y1 + y2 + ... ~ x1 + x2 + ... | g1 + g2 + ...} where \code{y1...ym} are response variables \code{1...m} (should be numeric), \code{x1...xi} are treatment variables \code{1...i} (should be categorical). \code{g1...gm} are grouping variables \code{1...g} (factors or character vectors) by which data will be split before analysis, also appearing on separate panels during graphing. This formula interface is similar to the lattice formula interface and in fact uses the lattice formula parser. On-the-fly transformations and other manipulations can be used.
}
 \item{control}{
an optional character string indicating which treatment level(s) is/are the "control" for each treatment variable if comparisons with the control are desired. Comparisons of each treatment using the \code{mcc()} function are the default. If nothing is specified, the first level in the treatment vector (preferably an ordered factor) is used.
}
 \item{times}{
an optional character string indicating the name of a time variable (which can be numeric or of any of the typical date and time formats). Only one time variable can be specified.
}
 \item{pool.variance}{
Not implemented but will allow for confidence intervals to be calculated with more complex models that pool variance.
}
 \item{block}{
an optional character string indicating the names of variables indicating blocking. Should be specified if subsample is specified as it indicates how averaging should be done. Otherwise this is not implemented but will be for use with analyses such as paired t-tests
}
 \item{replicate_id}{
an optional character string indicating the names of variables over which subsamples should be averaged. The \code{mean()} function is used for averaging.
}
 \item{average_subsamples}{
a logical argument indicating whether any remaining replicates should be averaged. This can be specified in the case that subsamples exist but no replicate id variable is specified.
}
 \item{summary_functions}{
a character vector including the name of one or multiple functions to be used to summarize each treatment at each time point within each group. Functions \code{se} and \code{CI68} are included. Others such as \code{SD} or custom functions can be used too.
}
 \item{comp_groups}{
a function that will be used to set up specific comparisons to be made. The two currently available functions are \code{mcc} and \code{allcomps}, specifying multiple comparisons with the control (the control group defined by the control argument) or all pairwise comparisons, respectively. Viewing the \code{mcc} and \code{allcomps} functions will show how a set of custom comparisons could be easily defined.
}
 \item{comp_function}{
a function that will be used to make the specific comparisons. Currently available options include \code{welchCI} for confidence intervals generated with a Welch t-test around the subtractive difference treatment - control; \code{bootdiff} for a bootstrapped confidence interval around the subtractive difference treatment - control; \code{bootfrac} for the treatment/control quotient with a bootstrapped confidence interval; and \code{bootperc} for a "percent difference," (treatment - control) / control * 100 with a bootstrap confidence interval. Custom functions can be defined and more are planned to be added, especially for paired tests. Ultimately, more complex analyses including approaches such as mixed models could be added.
}
 \item{extract}{
a list describing which components of the output of \code{comp_function} to extract. Each \code{comp_function} includes a default list.
}
 \item{conf.int}{
The alpha level for confidence interval calculations.
}
}

\details{
Using a formula interface, treateffect facilitates efficient calculation of the size of treatment effect sizes when comparing multiple categorical treatments. Standard variable types such as response variables (of which there can be multiple for efficient analysis), treatment categories, a time variable, a blocking variable, variables over which to pool the data, and panel variables over which to divide the data before analysis (e.g. two sites) can be specified. Functions used to calculate the size of the effect are flexible. The specific summary statistics (e.g., mean, SE, SD), comparison function (e.g., confidence interval from Welch t-test or a bootstrapped confidence interval) and the comparisons to perform (e.g., all pairwise comparisons or multiple comparisons with a control) can all be specified.
}

\value{
Returns an object of class "\code{te}". A \code{print} method shows the results in tabular form and the \code{plot} and \code{plotdiff}, can be used to plot the results.
 An object of class "\code{te}" is a list containing some or all of the following components:
\item{source_data}{the unmodified data supplied via the \code{data} argument.}
\item{design}{a list of all of the "design" features used to shape the analysis and plotting such as the identities of the treatment, response, time, and block variables}
\item{data}{the data frame used for analysis after restructuring to accommodate for example multiple response variables}
\item{treatment_summaries}{one of two main output data frames showing the output of the \code{summary_functions} functions that were specified for each variable for each treatment group at each time point in each group.}
\item{comparisons}{the second of two main output data frames showing the output of the analyses by the specified \code{comp_function} for each variable for each comparison specified by the \code{comp_groups} function at each time point in each group.}
}
\author{Anthony Darrouzet-Nardi}

\seealso{
\code{\link{plot.te}}, \code{\link{plotdiff}}, \code{\link{define_comparisons}}
}
\examples{
set.seed(100)
tst1 <- data.frame(response = c(rnorm(5), rnorm(5) + 2),
  treatment = ordered(rep(c("control", "treatment"), ea = 5)),
  block = factor(rep(1:5,2)))
tst1.te <- treateffect(tst1, response~treatment, comp_function = bootperc)
tst1.te
plot(tst1.te)
plotdiff(tst1.te)
}
