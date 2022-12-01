library(conflicted)
library(tidyverse)
library(progress)
library(reticulate)
library(here)
library(Rfast)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

source_python(here("ir_mc.py"))

# IR-MC
ir_mc <- function(d = c("Eucl", "Pear", "rand"), query, docs, mdata) {
  
  # Ensure that `query` has a corresponding melting curve in `mdata`
  mc_query <- tibble(Protein = query) %>% inner_join(mdata, by = "Protein")
  if (nrow(mc_query) == 1) {
    mc_query <- mc_query %>% pull(Curve) %>% unlist(.)
  } else {
    return(NULL)
  }
  
  # Ensure that at least one protein in `docs` has a corresponding melting curve in `mdata`
  mc_docs <- tibble(Protein = docs) %>% inner_join(mdata, by = "Protein")
  if (nrow(mc_docs) > 0) {
    docs1 <- mc_docs$Protein
    mc_docs <- matrix(unlist(mc_docs$Curve), byrow = T, nrow = length(mc_docs$Curve))
    rownames(mc_docs) <- docs1
  } else {
    return(NULL)
  }
  
  if (d[1] == "Eucl") {
    # Compute distances
    all_dists <- dista(t(as.matrix(mc_query)), mc_docs) %>% as.vector(.)
    names(all_dists) <- docs1
    # Generate ranking
    rk <- sort(all_dists) %>% names(.)
    
  } else if (d[1] == "Pear") {
    # Compute distances
    all_dists <- apply(mc_docs, 1, pearson_dissim, x2 = mc_query)
    # Generate ranking
    rk <- sort(all_dists) %>% names(.)
    
  } else if (d[1] == "rand") {
    # Generate ranking
    rk <- sample(docs1)
    
  } else {
    stop("Unrecognizable dissimilarity measure! For argument `d`, please enter one of:
    \"Eucl\", \"Pear\", or \"rand\".")
  }
  
  return(rk)
}

# IR-MC-ML
ir_mc_ml <- function(query, docs, mdata, known_ints) {
  
  # Ensure that `query` has a corresponding melting curve in `mdata`
  mc_query <- tibble(Protein = query) %>% inner_join(mdata, by = "Protein")
  if (nrow(mc_query) == 1) {
    mc_query <- mc_query %>% pull(Curve) %>% unlist(.)
  } else {
    return(NULL)
  }
  
  # Ensure that at least one protein in `docs` has a corresponding melting curve in 'mdata
  mc_docs <- tibble(Protein = docs) %>% inner_join(mdata, by = "Protein")
  if (nrow(mc_docs) > 0) {
    docs1 <- mc_docs$Protein
    mc_docs <- matrix(unlist(mc_docs$Curve), byrow = T, nrow = length(mc_docs$Curve))
    rownames(mc_docs) <- docs1
  } else {
    return(NULL)
  }
  
  # Ensure that at least one protein in `known_ints` has a corresponding melting curve in `mdata`
  mc_ints <- tibble(Protein = known_ints) %>% inner_join(mdata, by = "Protein")
  if (nrow(mc_ints) > 0) {
    known_ints1 <- mc_ints$Protein
    mc_ints <- matrix(unlist(mc_ints$Curve), byrow = T, nrow = length(mc_ints$Curve))
    rownames(mc_ints) <- known_ints1
  } else {
    return(NULL)
  }
  
  # Form input for ITML
  mc_itml <- rbind(mc_query, mc_ints)
  
  # Learn a matrix M by applying ITML mc_itml
  M <- learn_itml(mc_itml)
  
  # Compute the Mahalanobis distance between the query's melting curve and each document's melting
  # curve.
  all_dists <- mahalanobis(x = mc_docs, center = mc_query, cov = M, inverted = T)
  
  # Rank the documents from smallest distance to largest
  rk <- all_dists %>% sort(.) %>% names(.)
  
  return(rk)
}

# Pearson dissimilarity
pearson_dissim <- function(x1, x2) {
  (1 - cor(x1, x2, method = "pearson")) / 2
}
