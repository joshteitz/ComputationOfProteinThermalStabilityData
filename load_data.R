library(conflicted)
library(tidyverse)
library(readxl)
library(here)
library(progress)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

# paths to datasets
S1_TO_S18 <- here("datasets", "tabless1_to_s27", "aan0346_Tables_S1_to_S18.xlsx")
S19_TO_S27 <- here("datasets", "tabless1_to_s27", "aan0346__Tables_S19_to_S27.xlsx")
BIOPLEX <- here("datasets", "BioPlex_BaitPreyPairs_noFilters_293T_10K_Dec_2019.tsv")
BIOPLEX_PPIS <- here("datasets", "BioPlex_293T_Network_10K_Dec_2019.tsv")
GENEID_TO_UNIPROT <- here("datasets", "geneid_to_uniprot.tsv")

# Load melting data for a cell line
load_melting_curves <- function(cell_line = "K562") {
  
  path_sheet = switch(
    cell_line,
    "K562 lysate" = c(S1_TO_S18, "Table S6"),
    "K562" = c(S1_TO_S18, "Table S7"),
    "A375" = c(S19_TO_S27, "Table S19"),
    "HCT116" = c(S19_TO_S27, "Table S20"),
    "HEK293T" = c(S19_TO_S27, "Table S21"),
    "HL60" = c(S19_TO_S27, "Table S22"),
    "MCF7" = c(S19_TO_S27, "Table S23"),
    "mouse liver" = c(S19_TO_S27, "Table S24")
  )
  
  if (is.null(path_sheet)) {
    stop("Unrecognizable cell line! Please enter one of: \"K562 lysate\", \"K562\", \"A375\", \"HCT116\",
         \"HEK293T\", \"HL60\", \"MCF7\", or \"mouse liver\".")
  }
  
  mdata <- read_excel(path = path_sheet[1], sheet = path_sheet[2], skip = 2) %>%
    select(Protein = Accession, matches("^T\\d+$")) %>%
    select(-T37) %>%
    nest(Curve = starts_with("T")) %>%
    mutate(Curve = map(Curve, ~ .x %>% as.numeric)) %>%
    distinct(Protein, .keep_all = T)
  
  return(mdata)
}

# Load protein-protein interactions
# `num_pubs`: minimum number of publications a PPI must have to be returned
load_ppis <- function(num_pubs = 2) {
  
  # load data
  ppidata <- read_excel(path = S1_TO_S18, sheet = "Table S2", skip = 2) %>%
    select(ProteinA = `Protein A`, ProteinB = `Protein B`, Num_pubs = `No. of Publications`) %>%
    filter(Num_pubs >= num_pubs) %>%
    select(-Num_pubs)
  
  # # make sure ProteinA < ProteinB and make sure there are no duplicate PPIs
  # ppidata <- ppidata %>% mutate(Proteins = map2(ProteinA, ProteinB, ~ sort(c(.x, .y)))) %>%
  #   mutate(ProteinA = map_chr(Proteins, ~ .x[1])) %>%
  #   mutate(ProteinB = map_chr(Proteins, ~ .x[2])) %>%
  #   select(-Proteins) %>%
  #   distinct(ProteinA, ProteinB, .keep_all = T)
  
  # Convert PPI data to `Protein` and `Interactors` columns
  ppidata <- bind_rows(
    ppidata,
    ppidata %>% select(ProteinA = ProteinB, ProteinB = ProteinA)
  ) %>%
    nest(Interactors = ProteinB) %>%
    mutate(Interactors = map(Interactors, ~ .x$ProteinB %>% unique)) %>%
    select(Protein = ProteinA, Interactors)
  
  return(ppidata)
}

load_BioPlex_PPIs <- function() {
  
  # Load mapping from gene ID to Uniprot
  geneid_to_uniprot <- read_tsv(GENEID_TO_UNIPROT, show_col_types = F)
  
  # Load BioPlex PPI data with `pInt > .99`
  ppidata <- read_tsv(BIOPLEX_PPIS, show_col_types = F) %>%
    filter(pInt > .99) %>%
    select(ProteinA = UniprotA, ProteinB = UniprotB)
  
  # Convert PPI data to `Protein` and `Interactors` columns
  ppidata <- bind_rows(
    ppidata,
    ppidata %>% select(ProteinA = ProteinB, ProteinB = ProteinA)
  ) %>%
    nest(Interactors = ProteinB) %>%
    mutate(Interactors = map(Interactors, ~ .x$ProteinB %>% unique)) %>%
    select(Protein = ProteinA, Interactors)
  
  return(ppidata)
}

# Load BioPlex pull-down assay results
load_pulldown_results <- function() {
  
  # Load mapping from gene ID to Uniprot
  geneid_to_uniprot <- read_tsv(GENEID_TO_UNIPROT, show_col_types = F)
  
  # Load BioPlex pull-down assay data and convert gene IDs to Uniprot
  pdata <- read_tsv(BIOPLEX, show_col_types = F) %>%
    inner_join(geneid_to_uniprot, by = c("bait_geneid" = "GeneId")) %>%
    inner_join(geneid_to_uniprot, by = c("gene_id" = "GeneId"), suffix = c("Bait", "Prey")) %>%
    select(Bait = UniprotBait, Prey = UniprotPrey) %>%
    filter(Bait != Prey) %>%
    nest(Prey = Prey) %>%
    mutate(Prey = map(Prey, ~ .x$Prey %>% unique))
  
  return(pdata)
}

# mdata <- load_melting_data(cell_line = "HEK293T")
# ppis <- load_ppi_data(num_pubs = 2)
# 
# pdata1 <- pdata %>% filter(Bait %in% c("P00813", "P55039", "Q9NZL4", "Q6ICG6")) %>% select(Bait, Prey)
# 
# pdata1 <- pdata1 %>%
#   nest(Prey = Prey) %>%
#   inner_join(mdata, by = c("Bait" = "Protein")) %>%
#   select(Bait = )
#   mutate(Prey = map(Prey, ~ .x %>% inner_join(mdata, by = c("Prey" = "Protein")))) %>%
#   mutate(Prey = map(Prey, ~ {
#     m <- matrix(unlist(.x$Sols), byrow = T, nrow = length(.x$Sols));
#     rownames(m) <- .x$Prey;
#     m
#   }))
