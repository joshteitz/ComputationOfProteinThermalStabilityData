library(here)

source(here("load_data.R"))
source(here("ir_mc.R"))
source(here("ir_eval.R"))

###################################################################################################
## Information retrieval of melting curves tends to rank interactors higher than non-interactors ##
###################################################################################################

# load BioPlex pull-down assay results
pdata <- load_pulldown_results()

# load melting curves
mdata <- load_melting_curves(cell_line = "HEK293T")
nrow(mdata) # number of melting curves

# load PPIs cited by at least two publications
ppis <- load_ppis(num_pubs = 2)

# Restrict PPIs to those where both proteins have a melting curve
ppis <- ppis %>% 
  semi_join(mdata, by = "Protein") %>%
  mutate(Interactors = map(Interactors, ~ .x[.x %in% mdata$Protein])) %>%
  mutate(Num_ints = map_int(Interactors, length)) %>%
  filter(Num_ints > 0) %>%
  select(-Num_ints)
map_int(ppis$Interactors, length) %>% sum(.) / 2 # number of PPIs

# Remove any pull-down assay without at least one bait-prey interaction.
# `Interactors` are all prey that interact with the bait.
pdata <- pdata %>%
  inner_join(ppis, by = c("Bait" = "Protein")) %>%
  mutate(Interactors = map2(Interactors, Prey, ~ .x[.x %in% .y])) %>%
  mutate(Num_ints = map_int(Interactors, length)) %>%
  filter(Num_ints >= 1) %>%
  select(-Num_ints)
nrow(pdata) # number of remaining pull-down assays

# Apply IR-MC to pull-down assays

# Euclildean distance 
pb <- progress_bar$new(total = nrow(pdata))
res_eucl <- pdata %>% mutate(Rk = map2(Bait, Prey, ~ {
  pb$tick();
  ir_mc(d = "Eucl", query = .x, docs = .y, mdata)
})) %>%
  select(-Prey)

# Pearson dissimilarity
pb <- progress_bar$new(total = nrow(pdata))
res_pear <- pdata %>% mutate(Rk = map2(Bait, Prey, ~ {
  pb$tick();
  ir_mc(d = "Pear", query = .x, docs = .y, mdata)
})) %>%
  select(-Prey)
  
# Random ranking
set.seed(28)
pb <- progress_bar$new(total = nrow(pdata))
res_rand <- pdata %>% mutate(Rk = map2(Bait, Prey, ~ {
  pb$tick();
  ir_mc(d = "rand", query = .x, docs = .y, mdata)
})) %>%
  select(-Prey)

# Compute average precision

# Euclidean distance
pb <- progress_bar$new(total = nrow(res_eucl))
res_eucl <- res_eucl %>%
  mutate(AP = map2_dbl(Rk, Interactors, ~ {pb$tick(); avg_prec(.x, .y)}))

# Pearson dissimilarity
pb <- progress_bar$new(total = nrow(res_pear))
res_pear <- res_pear %>%
  mutate(AP = map2_dbl(Rk, Interactors, ~ {pb$tick(); avg_prec(.x, .y)}))

# Random ranking
pb <- progress_bar$new(total = nrow(res_rand))
res_rand <- res_rand %>%
  mutate(AP = map2_dbl(Rk, Interactors, ~ {pb$tick(); avg_prec(.x, .y)}))

# Paired t-test comparing mAP of Euclidean distance to mAP of Pearson dissimilarity
t.test(x = res_eucl$AP, y = res_pear$AP, alternative = "two.sided", paired = T)

# Paired t-test comparing mAP of Euclidean distance to mAP of random ranking
t.test(x = res_eucl$AP, y = res_rand$AP, alternative = "two.sided", paired = T)

# Compute improvement probability
set.seed(28)
pb <- progress_bar$new(total = nrow(res_eucl))
res_eucl <- res_eucl %>%
  mutate(IP = map2_dbl(Rk, Interactors, ~ {pb$tick(); improv_prob(.x, .y, 1000)}))

# Number of IR-MC rankings with improvement probability greater than .8 and .95
res_eucl %>% filter(IP > .8) %>% nrow(.)
res_eucl %>% filter(IP > .95) %>% nrow(.)

# histogram of improvement probability
ggplot(res_eucl, aes(x = IP)) + 
  geom_bar() +
  scale_y_continuous(breaks = seq(0, 500, by = 50)) +
  scale_x_binned(
    breaks = seq(0, 1, by = .05),
    labels = c("0", "", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5",
               "", "0.6", "", "0.7", "", "0.8", "", "0.9", "", "1.0")
  ) +
  geom_hline(yintercept = nrow(res_eucl) / 20, linetype = "dashed", color = "black", size = 3) +
  ylab("# of IR-MC rankings") +
  xlab("Improvement probability") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 28, color = "black"),
    axis.text.x = element_text(size = 28, color = "black"),
    axis.title = element_text(size = 32),
    panel.grid.major.y = element_line(colour = "grey", size = 1),
    panel.grid.major.x = element_line(colour = "grey", size = 1)
  )

###################################################################################################
## Quantifying the number of experiments necessary to find one interactor and all interactors    ##
###################################################################################################

# Standardize each pull-down assay to contain exactly five interactors and 45 non-interactors

# `Non_ints` are all prey that are not interactors and have melting curves
pdata1 <- pdata %>%
  mutate(Non_ints = map2(Prey, Interactors, ~ .x[!(.x %in% .y)])) %>%
  mutate(Non_ints = map(Non_ints, ~ .x[.x %in% mdata$Protein]))

# Remove any pull-down assay that does not contain > 5 interactors and > 45 non-interactors
pdata1 <- pdata1 %>%
  mutate(Num_ints = map_int(Interactors, length)) %>%
  mutate(Num_non_ints = map_int(Non_ints, length)) %>%
  filter(Num_ints >= 5 & Num_non_ints >= 45) %>%
  select(-Num_ints, -Num_non_ints)

# Standardize each of the remaining pull-down assays.
set.seed(28)
pdata1 <- pdata1 %>%
  mutate(Prey = map(Non_ints, ~ sample(.x, 45))) %>%
  mutate(Interactors = map(Interactors, ~ sample(.x, 5))) %>%
  mutate(Prey = map2(Interactors, Prey, ~ c(.x, .y))) %>%
  select(-Non_ints)

# Apply IR-MC to each standardized pull-down assay
pb <- progress_bar$new(total = nrow(pdata1))
res_eucl1 <- pdata1 %>%
  mutate(Rk = map2(Bait, Prey, ~ {pb$tick(); ir_mc(d = "Eucl", query = .x, docs = .y, mdata)})) %>%
  select(Bait, Interactors, Rk)

# Randomly rank each standardized pull-down assay
set.seed(28)
pb <- progress_bar$new(total = nrow(pdata1))
res_rand1 <- pdata1 %>%
  mutate(Rk = map2(Bait, Prey, ~ {pb$tick(); ir_mc(d = "rand", query = .x, docs = .y, mdata)})) %>%
  select(Bait, Interactors, Rk)

# Number of experiments needed to find one bait-prey interaction and all five bait-prey
# interactions for each IR-MC ranking
res_eucl1 <- res_eucl1 %>%
  mutate(Rk_interactors = map2(Interactors, Rk, ~ which(.y %in% .x))) %>%
  mutate(Num_expts = map_int(Rk_interactors, ~ .x[1])) %>%
  mutate(Num_expts_all = map_int(Rk_interactors, ~ .x[5])) %>%
  select(-Rk_interactors)

# Number of experiments needed to find one bait-prey interaction and all five bait-prey
# interactions for each random ranking
res_rand1 <- res_rand1 %>%
  mutate(Rk_interactors = map2(Interactors, Rk, ~ which(.y %in% .x))) %>%
  mutate(Num_expts = map_int(Rk_interactors, ~ .x[1])) %>%
  mutate(Num_expts_all = map_int(Rk_interactors, ~ .x[5])) %>%
  select(-Rk_interactors)

# Scatter plots of number of experiments needed to find one and all five interactors
# for standardized pull-down assays

# plot data
plot_data <- bind_rows(
  tibble(Method = "IR-MC", Num_ints = 1L, Num_expts = res_eucl1$Num_expts),
  tibble(Method = "IR-MC", Num_ints = 5L, Num_expts = res_eucl1$Num_expts_all),
  tibble(Method = "Random", Num_ints = 1L, Num_expts = res_rand1$Num_expts),
  tibble(Method = "Random", Num_ints = 5L, Num_expts = res_rand1$Num_expts_all)
) %>%
  mutate(Category = paste(Method, Num_ints)) %>%
  select(Category, everything()) %>%
  mutate(Category = factor(Category, levels = c("Random 5", "IR-MC 5", "Random 1", "IR-MC 1")))

# ggplot default color scheme
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# scatter plot
ggplot(plot_data, aes(y = Category, x = Num_expts, color = Method)) +
  geom_jitter(size = 8, alpha = .4) +
  # geom_vline(xintercept = 4.97) +
  # geom_vline(xintercept = 8.40) +
  # geom_vline(xintercept = 32.2) +
  # geom_vline(xintercept = 42.5) +
  scale_color_manual(values = gg_color_hue(2) %>% rev(.)) +
  scale_y_discrete(labels = c("IR-MC 1" = "1", "Random 1" = "1", "IR-MC 5" = "5", "Random 5" = "5")) +
  xlab("Number of experiments required") +
  ylab("Number of bait-prey interactions desired") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 28, color = "black"),
    axis.title = element_text(size = 32),
    panel.grid.major.y = element_line(colour = "grey", size = 1),
    legend.title = element_text(size = 32, color = "black"),
    legend.text = element_text(size = 28, color = "black"),
    legend.position = "none"
  )

###################################################################################################
## IR-MC-ML tends to rank interactors higher than IR-MC                                          ##
###################################################################################################

# Remove pull-down assays with five or fewer interactors.
pdata <- pdata %>%
  mutate(Num_ints = map_int(Interactors, length)) %>%
  filter(Num_ints > 5) %>%
  select(-Num_ints)
nrow(pdata)

# Modify each pull-down assay by randomly sampling five prey proteins categorized as interactors
# and extracting them. Store these extracted proteins in a column called `Known_ints`.
set.seed(28)
pdata <- pdata %>%
  mutate(Known_ints = map(Interactors, ~ sample(.x, 5))) %>%
  mutate(Interactors = map2(Interactors, Known_ints, ~ .x[!(.x %in% .y)])) %>%
  mutate(Prey = map2(Prey, Known_ints, ~ .x[!(.x %in% .y)]))

# Apply IR-MC to each pull-down assay
pb <- progress_bar$new(total = nrow(pdata))
res_ir_mc <- pdata %>%
  mutate(Rk = map2(Bait, Prey, ~ {pb$tick(); ir_mc(d = "Eucl", query = .x, docs = .y, mdata)})) %>%
  select(Bait, Interactors, Rk)

# Apply IR-MC-ML to each pull-down assay
pb <- progress_bar$new(total = nrow(pdata))
res_ir_mc_ml <- pdata %>%
  mutate(Rk = pmap(list(Bait, Prey, Known_ints), ~ {
    pb$tick();
    ir_mc_ml(query = ..1, docs = ..2, mdata, known_ints = ..3)
  })) %>%
  select(Bait, Interactors, Rk)

# Compute AP of each IR-MC ranking
pb <- progress_bar$new(total = nrow(pdata))
res_ir_mc <- res_ir_mc %>%
  mutate(AP = map2_dbl(Rk, Interactors, ~ {pb$tick(); avg_prec(.x, .y)}))

# Compute IP of each IR-MC ranking
pb <- progress_bar$new(total = nrow(pdata))
res_ir_mc <- res_ir_mc %>%
  mutate(IP = map2_dbl(Rk, Interactors, ~ {pb$tick(); improv_prob(.x, .y, 1000)}))

# Compute AP of each IR-MC-ML ranking
pb <- progress_bar$new(total = nrow(pdata))
res_ir_mc_ml <- res_ir_mc_ml %>%
  mutate(AP = map2_dbl(Rk, Interactors, ~ {pb$tick(); avg_prec(.x, .y)}))

# Compute IP of each IR-MC-ML ranking
pb <- progress_bar$new(total = nrow(pdata))
res_ir_mc_ml <- res_ir_mc_ml %>%
  mutate(IP = map2_dbl(Rk, Interactors, ~ {pb$tick(); improv_prob(.x, .y, 1000)}))

# Standardize each pull-down assay to contain exactly five interactors and 45 non-interactors

# `Non_ints` are all prey that are not interactors and have melting curves
pdata_stand <- pdata %>%
  mutate(Non_ints = map2(Prey, Interactors, ~ .x[!(.x %in% .y)])) %>%
  mutate(Non_ints = map(Non_ints, ~ .x[.x %in% mdata$Protein]))

# Remove any pull-down assay that does not contain > 5 interactors and > 45 non-interactors
pdata_stand <- pdata_stand %>%
  mutate(Num_ints = map_int(Interactors, length)) %>%
  mutate(Num_non_ints = map_int(Non_ints, length)) %>%
  filter(Num_ints >= 5 & Num_non_ints >= 45) %>%
  select(-Num_ints, -Num_non_ints)

# Standardize each of the remaining pull-down assays.
set.seed(28)
pdata_stand <- pdata_stand %>%
  mutate(Prey = map(Non_ints, ~ sample(.x, 45))) %>%
  mutate(Interactors = map(Interactors, ~ sample(.x, 5))) %>%
  mutate(Prey = map2(Interactors, Prey, ~ c(.x, .y))) %>%
  select(-Non_ints)

# Apply IR-MC to each standardized pull-down assay
pb <- progress_bar$new(total = nrow(pdata_stand))
res_ir_mc_stand <- pdata_stand %>%
  mutate(Rk = map2(Bait, Prey, ~ {pb$tick(); ir_mc(d = "Eucl", query = .x, docs = .y, mdata)})) %>%
  select(Bait, Interactors, Rk)

# Apply IR-MC-ML to each pull-down assay
pb <- progress_bar$new(total = nrow(pdata_stand))
res_ir_mc_ml_stand <- pdata_stand %>%
  mutate(Rk = pmap(list(Bait, Prey, Known_ints), ~ {
    pb$tick();
    ir_mc_ml(query = ..1, docs = ..2, mdata, known_ints = ..3)
  })) %>%
  select(Bait, Interactors, Rk)

# Number of experiments needed to find one bait-prey interaction and all five bait-prey
# interactions for each IR-MC ranking
res_ir_mc_stand <- res_ir_mc_stand %>%
  mutate(Rk_interactors = map2(Interactors, Rk, ~ which(.y %in% .x))) %>%
  mutate(Num_expts = map_int(Rk_interactors, ~ .x[1])) %>%
  mutate(Num_expts_all = map_int(Rk_interactors, ~ .x[5])) %>%
  select(-Rk_interactors)

# Number of experiments needed to find one bait-prey interaction and all five bait-prey
# interactions for each IR-MC-ML ranking
res_ir_mc_ml_stand <- res_ir_mc_ml_stand %>%
  mutate(Rk_interactors = map2(Interactors, Rk, ~ which(.y %in% .x))) %>%
  mutate(Num_expts = map_int(Rk_interactors, ~ .x[1])) %>%
  mutate(Num_expts_all = map_int(Rk_interactors, ~ .x[5])) %>%
  select(-Rk_interactors)

## Euclidean distance ##

# Run IR-MC with Euclidean distance
rk_eucl <- ir_mc("Eucl", pdata, mdata)

# join interactors so that rankings can be evaluated
rk_eucl <- rk_eucl %>% inner_join(pdata, by = "Bait") %>% select(-Prey)
  
# Evaluate each ranking by average precision
pb <- progress_bar$new(total = nrow(rk_eucl))
rk_eucl <- rk_eucl %>% mutate(Avg_prec = map2_dbl(Rk, Interactors, ~{pb$tick(); avg_prec(.x, .y)}))

# Evaluate each ranking by improvement probability
set.seed(28)
pb <- progress_bar$new(total = nrow(rk_eucl))
rk_eucl <- rk_eucl %>% mutate(Improv_prob = map2_dbl(Rk, Interactors, ~ {pb$tick(); improv_prob(.x, .y, 1000)}))

## Pearson dissimilarity ##

# Run IR-MC with Pearson dissimilarity
rk_pear <- ir_mc("Pear", pdata, mdata)

# join interactors so that rankings can be evaluated
rk_pear <- rk_pear %>% inner_join(pdata, by = "Bait") %>% select(-Prey)

# Evaluate each ranking by average precision
pb <- progress_bar$new(total = nrow(rk_pear))
rk_pear <- rk_pear %>% mutate(Avg_prec = map2_dbl(Rk, Interactors, ~{pb$tick(); avg_prec(.x, .y)}))

# Evaluate each ranking by improvement probability
set.seed(28)
pb <- progress_bar$new(total = nrow(rk_pear))
rk_pear <- rk_pear %>% mutate(IP = map2_dbl(Rk, Interactors, ~ {pb$tick(); improv_prob(.x, .y, 1000)}))

## Random ranking ##
set.seed(28)
rk_rand <- ir_mc("Rand", pdata, mdata)

# join interactors so that rankings can be evaluated
rk_rand <- rk_rand %>% inner_join(pdata, by = "Bait") %>% select(-Prey)

# Evaluate each ranking by average precision
pb <- progress_bar$new(total = nrow(rk_rand))
rk_rand <- rk_rand %>% mutate(Avg_prec = map2_dbl(Rk, Interactors, ~{pb$tick(); avg_prec(.x, .y)}))

# Evaluate each ranking by improvement probability
set.seed(28)
pb <- progress_bar$new(total = nrow(rk_rand))
rk_rand <- rk_rand %>% mutate(Improv_prob = map2_dbl(Rk, Interactors, ~ {pb$tick(); improv_prob(.x, .y, 1000)}))

# paired t-tests
t.test(rk_eucl$Avg_prec, rk_rand$Avg_prec, alternative = "two.sided", paired = T)
t.test(rk_pear$Avg_prec, rk_rand$Avg_prec, alternative = "two.sided", paired = T)
t.test(rk_eucl$Avg_prec, rk_pear$Avg_prec, alternative = "two.sided", paired = T)

## Standardize the pull-down assays to each have 50 prey in total, exactly five of which are interactors.

# `Non_ints` are all prey that are not interactors and have melting curves
pdata1 <- pdata %>%
  mutate(Non_ints = map2(Prey, Interactors, ~ .x[!(.x %in% .y)])) %>%
  mutate(Non_ints = map(Non_ints, ~ .x[.x %in% mdata$Protein]))

# Remove any pull-down assay that does not contain > 5 interactors and > 45 non-interactors
pdata1 <- pdata1 %>%
  mutate(Num_ints = map_int(Interactors, length)) %>%
  mutate(Num_non_ints = map_int(Non_ints, length)) %>%
  filter(Num_ints >= 5 & Num_non_ints >= 45) %>%
  select(-Num_ints, -Num_non_ints)

# Standardize each of the remaining pull-down assays.
set.seed(28)
pdata1 <- pdata1 %>%
  mutate(Prey = map(Non_ints, ~ sample(.x, 45))) %>%
  mutate(Interactors = map(Interactors, ~ sample(.x, 5))) %>%
  mutate(Prey = map2(Interactors, Prey, ~ c(.x, .y))) %>%
  select(-Non_ints)

## Euclidean distance

# Run IR-MC with Euclidean distance
rk_eucl1 <- ir_mc("Eucl", pdata1, mdata)

# join interactors so that rankings can be evaluated
rk_eucl1 <- rk_eucl1 %>% inner_join(pdata1, by = "Bait") %>% select(-Prey)

# Ranks of interactors
rk_eucl1 <- rk_eucl1 %>% mutate(Rk_ints = map2(Rk, Interactors, ~ which(.x %in% .y)))

## Random ranking

# Run IR-MC with random distance
set.seed(28)
rk_rand1 <- ir_mc("Rand", pdata1, mdata)

# join interactors so that rankings can be evaluated
rk_rand1 <- rk_rand1 %>% inner_join(pdata1, by = "Bait") %>% select(-Prey)

# Ranks of interactors
rk_rand1 <- rk_rand1 %>% mutate(Rk_ints = map2(Rk, Interactors, ~ which(.x %in% .y)))

## Plot number of experiments to find one and all five interactors for standardized pull-down assays.

# Construct data frame
df <- bind_rows(
  tibble(
    Category = 0L,
    Num_expts = map_int(rk_eucl1$Rk_ints, ~ .x[1]),
  ),
  # tibble(
  #   Category = 5L,
  #   Num_expts = map_int(rk_rand1$Rk_ints, ~ .x[1]),
  #   Color = "2"
  # ),
  tibble(
    Category = 5,
    Num_expts = map_int(rk_eucl1$Rk_ints, ~ .x[5])
  )
  # tibble(
  #   Category = 15L,
  #   Num_expts = map_int(rk_rand1$Rk_ints, ~ .x[5]),
  #   Color = "2"
  # )
)

ggplot(df, aes(y = Category, x = Num_expts)) +
  geom_jitter(size = 8, alpha = .4, color = "#00BFC4") +
  # geom_segment(aes(x = -2.2, y = 4.973333, xend = 2.2, yend = 4.973333), size = 3, color = "#336699", linetype = "dashed") +
  # geom_segment(aes(x = 5-2.2, y = 8.4031, xend = 5+2.2, yend = 8.4031), size = 3, color = "#CC3300", linetype = "dashed") +
  # # geom_point(aes(x = 10, y = 32.18476), size = 8, color = "#336699", shape = 23, fill = "#336699") +
  # geom_segment(aes(x = 10-2.2, y = 32.18476, xend = 10+2.2, yend = 32.18476), size = 3, color = "#336699", linetype = "dashed") +
  # geom_point(aes(x = 15, y = 42.51048), size = 8, color = "#CC3300", shape = 23, fill = "#CC3300") +
  scale_y_reverse(breaks = c(0,5), labels = c(1, 5)) +
  xlab("# of experiments") +
  ylab("# of interactors") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 28, color = "black"),
    axis.title = element_text(size = 32),
    panel.grid.major.y = element_line(colour = "grey", size = 1),
    legend.title = element_text(size = 28, color = "black"),
    legend.text = element_text(size = 22, color = "black")
  )
