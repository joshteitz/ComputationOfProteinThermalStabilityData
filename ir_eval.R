library(conflicted)
library(tidyverse)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

# Compute precision at k.
prec_at_k <- function(rk, rel, k) {
  sum(rk[1:k] %in% rel) / k
}

# Compute average precision
avg_prec <- function(rk, rel) {
  k_vals <- which(rk %in% rel)
  prec_at_k_vals <- map_dbl(k_vals, ~ prec_at_k(rk, rel, .x))
  ap <- sum(prec_at_k_vals) / length(k_vals)
  return(ap)
}

# Compute improvement probability
improv_prob <- function(rk, rel, num) {
  
  # average precision of input ranking
  ap <- avg_prec(rk, rel)
  
  # Generate random rankings
  rand_rks <- replicate(num, sample(rk))
  
  # average precisions of random rankings
  # pb <- progress_bar$new(total = num)
  rand_aps <- map_dbl(1:num, ~ {avg_prec(rand_rks[,.x], rel)})
  
  # improvement probability
  ip <- sum(ap > rand_aps) / num
  
  return(ip)
}