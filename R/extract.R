pmean <- function(a,b,c, units, digits = 2) paste0(round(a, digits = digits),
  units, " [", round(b, digits = digits),", ", round(c, digits = digits), "]")

fracdifftable <- function(x, factor = 1, units = '', digits = 2)
  mutate(x, fracdiff = pmean(fracdiff*factor, fracdiff_lo*factor,
    fracdiff_hi*factor, units, digits = digits)) %>%
  select(-fracdiff_lo, -fracdiff_hi)

summex <- function(x, ..., cen = "mean", pm = "se", digits = 2, units) {
  ex <- x$treatment_summaries %>%
    filter_(.dots = lazyeval::lazy_dots(...))
  paste0(round(ex[cen], digits = digits),
    " \u00b1 ", round(ex[pm], digits = digits), units)
}

#npool.te %>% summex(variable == "labile", tussock == "T", snow == "A", doy == 148, units = " ug N g-1 soil")

compex_fracdiff <- function(x, ..., cen = "fracdiff", lo = "fracdiff_lo",
  hi = "fracdiff_hi", digits = 2, units = "", factor = 1) {
  ex <- x$comparisons %>%
    filter_(.dots = lazyeval::lazy_dots(...))
  paste0(round(ex[cen] * factor, digits = digits),
    units, " [", round(ex[lo] * factor, digits = digits),", ",
    round(ex[hi] * factor, digits = digits), "]")
}

#n.te %>% compex_fracdiff(variable == "labile", tussock == "T", comparison == "A.N - C.N", doy == 148, digits = 0)
