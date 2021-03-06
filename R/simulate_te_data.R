tedatasim <- function(n, na = 0, response = 1, time = 1, groups = 1, levels = 2,
  subsample = 1, block = FALSE) {

  resp_var = matrix(rnorm(n * response * prod(groups) * prod(levels)
    * time * subsample), ncol = response) %>% data.frame
  colnames(resp_var) <- paste0(rep("resp_var", response), 1:response)
  resp_var <- as.data.frame(lapply(resp_var, function(cc) cc[sample(c(TRUE, NA),
    prob = c(1-na, na), size = length(cc), replace = TRUE) ]))

  splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
  code_num <- c(levels, groups)
  group_vectors <- splitAt(rev(letters[1:sum(code_num)]), cumsum(code_num) + 1)
  units <- expand.grid(c(group_vectors, list(1:time)))
  units <- units[rep(1:nrow(units), each = n * subsample),]
  row.names(units) <- 1:nrow(units)
  colnames(units) <- c(
    paste0(rep("pred_var", length(levels)), 1:length(levels)),
    paste0(rep("group_var", length(groups)), 1:length(groups)),
    "time_var")
  if(time < 2) units$time_var <- NULL
  if(!any(groups > 1)) units$group_var1 <- NULL

  if(subsample > 1)
    subsample_var = data.frame(subsample_var = rep(1:subsample,
      times = prod(groups) * time * prod(levels) * n)) else subsample_var = NULL

  if(block)
    block_var = data.frame(block_var =
      rep(rev(letters[10:(n+9)]),
      times = time * prod(groups) * prod(levels), ea = subsample)
      ) else
      block_var = NULL

  k <- eval(lapply(as.list(c("units", "block_var",
    "subsample_var", "resp_var")), FUN = function(x) eval(parse(text = x))))
  k <- k[!unlist(lapply(k, is.null))]
  m <- bind_cols(k) %>% tbl_df

#add some treatment/group/block differences
l <- data.frame(letters, numbers = 1:26, stringsAsFactors = FALSE)
new <- m
new[] <- lapply(m, as.character)
new <- dplyr::select(new, -contains("resp"), -contains("time"), -contains("sample"))
new[] <- l$numbers[match(unlist(new), l$letters)]
s <- rowSums(new)
if ("time_var" %in% colnames(m)) s <- s + m$time_var
for (r in 1:response) {m[names(resp_var)[r]] <-
  m[names(resp_var)[r]] + s + r * 3}
m
}
