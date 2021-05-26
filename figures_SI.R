# Generate figures in the Supplement of Miller and Allesina, "Metapopulations with habitat modification"
# zachmiller@uchicago.edu

source("functions.R")

library("deSolve")
library("tidyverse")
library("scales")
library("RColorBrewer")
library("ggpubr")

# set up color palettes to use throughout
light <- brewer.pal(n = 8, name = "Paired")[2*(1:5) - 1]
dark <- brewer.pal(n = 8, name = "Paired")[2*(1:5)]


##### Figure S1 #####

# Comparison of diversity-robustness relationship for stable vs. unstable P matrices

n <- 5 # number of species
num_mats <- 500 # how many matrices to sample for each panel (stable vs. unstable)

results <- tibble(id = numeric(), k = numeric(), type = numeric(), max_m = numeric())

# keep count of stable (i) and unstable (j) matrices, and keep sampling until we find num_mats of each type 
i <- 0
j <- 0
count <- 0 # keep a counter to i.d. each matrix
while(i < num_mats | j < num_mats) {
  count <- count + 1
  success <- FALSE
  while(!success) { # keep the next sampled P if its type is still need
    
    # sample symmetric P matrices without regard for stability (but require feasibility)
    out <- sample_P(n, check.feas = TRUE, check.stab = FALSE, symmetric = TRUE)
    type <- out$stab
    P <- out$P
    
    if (type) { # if stable
      if (i < num_mats) success <- TRUE
    } else { # else unstable
      if (j < num_mats) success <- TRUE
    }
  }
  
  # loop over submatrices and check feasibility. If any sub-community is not feasible, discard. 
  # Else, keep and record m_max along assembly sequence.  
  for(k in 1:n) {
    success <- FALSE
    P_k <- P[1:k, 1:k]
    if (all(rowSums(solve(P_k)) > 0)) { # if feasible (i.e. individual feasibility condition holds)
      results <- results %>% add_row(id = count, k = k, type = type, max_m = 1/sum(solve(P_k)))
      success <- TRUE
    } else { 
      break # exit loop and don't use this matrix
    }
  }
  
  if (success) { # advance count if the assembly sequence was feasible
    if (type) i <- i + 1 else j <- j + 1
  }
}

# process output
results <- results %>% 
  group_by(id, type) %>% 
  mutate(length = n()) %>% # length of each assembly sequence
  filter(length == n) %>% # discard incomplete sequences (i.e. where some sub-communities were not feasible)
  mutate(normalized_m = max_m / max_m[1], # compute normalized measure of m_max
         type = ifelse(type, "Stable", "Unstable"), # label stability
         increasing = all(diff(max_m) >= 0)) # classify assembly sequences as non-decreasing or not

psi1 <- results %>% ggplot() + 
  aes(x = k, y = normalized_m, color = increasing) + 
  geom_line(aes(group = id), alpha = 0.2) + 
  geom_hline(yintercept = 1, linetype = "dashed") + # show baseline for reference
  geom_boxplot(aes(group = k), outlier.shape = NA, fill = NA, size = 0.75, width = 0.2, color = "black") + # add boxplots to summarize spread
  facet_wrap(.~type) +
  xlab("Richness") + 
  ylab(expression("Normalized m"[max])) +
  scale_color_manual(values = dark) + 
  theme_classic() +
  theme(legend.position = "none")

show(psi1)


##### Figure S2 #####

# 2-panel figure showing long-term dynamics for cyclic P across different n and m

c <- 1 # common colonization rate
max_n <- 10 # maximum number of species
n_ms <- 30 # how many ms to use

results <- tibble(n = numeric(), m = numeric(), result = character())
# loop over number of species and values of m, classifying dynamics at each combination
for(n in 2:max_n) {
  
  P <- build_cyclic_P(n, c)
  max_m <- c / n
  
  # loop over m values
  for(j in 1:n_ms) {
    
    m <- (j/(n_ms + 1)) * max_m # generate an even sequence of m values up to max_m
    type <- simulate_and_classify_dynamics(P, rep(m, n)) # classify long-term dynamics
    
    results <- results %>% add_row(n = n, m = m, result = type)
  }
}

# compute the analytic stability threshold from SI Eq. 57 (in terms of absolute m value)
cycle_threshold_absolute <- tibble(n = seq(2, max_n, length.out = 100), c = rep(c, 100)) %>% 
  mutate(m_c = c / (n * (cos((2 * pi)/n) + 2)))
# compute the analytic stability threshold from SI Eq. 57 (in terms of fraction of m_max)
cycle_threshold_relative <- cycle_threshold_absolute %>% mutate(m_c = n * m_c)

# Scatter plot showing outcomes across n vs. absolute m values
pa <- results %>% ggplot() + 
  aes(x = m, y = n) + 
  geom_point(aes(color = result), size = 2) + 
  geom_line(data = cycle_threshold_absolute, aes(x = m_c, y = n), size = 1) + # overlay analytic threshold
  scale_color_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "gray")) + 
  xlab("Local extinction rate, m") +
  ylab("Number of species, n") + 
  guides(color = FALSE) + 
  theme_classic()

# Tile plot showing outcomes across n vs. relative m values
pb <- results %>% group_by(n) %>%
  mutate(`Fraction of max m` = round(m * (n / c), 3)) %>%
  ungroup() %>%
  ggplot() + 
  aes(x = `Fraction of max m`, y = n) + 
  geom_tile(aes(fill = result), color = "white", size = 0.25) +
  geom_line(data = cycle_threshold_relative, aes(x = m_c, y = n), size = 1) + # overlay analytic threshold 
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "gray")) + 
  ylab("Number of species, n") + 
  theme_classic() + 
  theme(legend.title = element_blank())

# compose panels and plot
psi2 <- ggarrange(pa, pb, ncol = 2, labels = c("A", "B"), widths = c(1, 1.5), legend = "right")
show(psi2)


##### Figure S3 #####

# Bifurcation diagram for cyclic P

n <- 4 # number of species
c <- 1 # common colonization rate

P <- build_cyclic_P(n, c)

m_range <- seq(0.001, c / n - 0.001, by = 0.001) # consider a wide range of m values from near 0 to m_max

results <- tibble(m = numeric(), lower = numeric(), upper = numeric())
for(m in m_range) {
  
  # compute equilibrium frequencies
  ystar <- (m / c) * rep(1, n)
  xstar <- (1 / n - m/c) * rep(1, n)
  zstar <- c(xstar, ystar)
  
  # initialize dynamics near equilibrium
  z0 <- zstar * runif(2 * n, 0.9, 1.1)
  z0 <- z0 / sum(z0)
  
  # integrate for along time to ensure convergence to equilibrium / cycle
  times <- seq(0, 5000, by = 1)
  
  # numerically integrate dynamics
  out <- ode(z0, times, evaluate_base_model, parms = list(m = m, P = P, n = n))
  
  # pull out min and max value (of species 1) from the second half of dynamics (after transients have passed)
  lower <- min(out[mean(times):max(times), 2])
  upper <- max(out[mean(times):max(times), 2])
  
  results <- results %>% add_row(m = m, lower = lower, upper = upper)
}

psi3 <- results %>% 
  gather(key = type, value = Frequency, -m) %>% 
  ggplot() +
  aes(x = m, y = Frequency, group = type) + # plot min and max frequencies for different values of m
  geom_rect(xmin = min(results %>% filter(upper < THRESH) %>% pull(m)), # highlight range of extinctions with gray box
            xmax = max(results %>% filter(upper < THRESH) %>% pull(m)),
            ymin = -Inf,
            ymax = Inf,
            fill = "gray",
            alpha = 0.2) + 
  geom_abline(slope = -1, intercept = 0.25, linetype = "dashed", color = "blue") + # indicate equilibrium with blue dashed line
  geom_vline(xintercept = c / (n * (cos((2 * pi) / n) + 2)), color = "red") + # indicate stability threshold with red line
  geom_hline(yintercept = 0) + # show zero for reference
  geom_point(size = 1) +
  xlab("Extinction rate (m)") +
  theme_classic()

show(psi3)


##### Figure S4 #####

# Show long-term dynamics for many random nonsymmetric P, at different values of m

set.seed(46)

n <- 3 # number of species
n_mats <- 200 # how many matrices to sample
n_ms <- 50 # how many ms to use
eps <- 0.01 # small parameter for finding P stable near 0 or near m_max

results <- tibble(mat = numeric(), m = numeric(), j = numeric(), result = character(), start = numeric())
for(i in 1:(2 * n_mats)) {
  
  # find feasible P stable at LOW (first 1:n_mats) or HIGH (second 1:n_mats) m 
  if (i <= n_mats) scale_m <- eps else scale_m <- 1 - eps
  P <- sample_P(n, m = "proportion", eps = scale_m)$P
  
  # loop over m values and classify long-term dynamics
  for(j in 1:n_ms) {
    
    m <- (j/(n_ms + 1)) * max_m # generate an even sequence of m values from 0 to m_max
    type <- simulate_and_classify_dynamics(P, rep(m, n))
    
    results <- results %>% add_row(mat = i, m = m, j = j, result = type, start = scale_m)
  }
}

results %>% 
  group_by(mat) %>%
  mutate(non_eq = ifelse(result == "equilibrium", n_ms + 1, j), # find m value at which each matrix loses stability (for a nice ordering)
         first_non_eq = jitter(min(non_eq), amount = 0.001)) %>% # jitter to break ties
  ungroup() %>%
  group_by(start) %>% # order separately for high and low eps 
  mutate(mat = dense_rank(first_non_eq), # rank matrices by where they lose stability
         `Fraction of max m` = j / (n_ms + 1)) %>%
  ungroup() %>%
  ggplot() + 
  aes(x = `Fraction of max m`, y = mat, fill = result) + 
  geom_tile() +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#CC79A7", "gray")) + 
  ylab("Realizations of P") + 
  facet_wrap(.~start) + 
  theme_classic() + 
  theme(legend.title = element_blank()) + 
  theme(legend.position = "bottom")