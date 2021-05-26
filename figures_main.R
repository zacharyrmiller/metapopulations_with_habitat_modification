# Generate figures in the Main Text of Miller and Allesina, "Metapopulations with habitat modification"
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


##### Figure 1 #####

# Conceptual figure illustrating metapopulation model with patch memory effects
# Panels A and B are generated separately; the following code generates example time series shown in Fig. 1C

n <- 3 # number of species
m <- 1 # common local extinction rate
P <- matrix(c(2.2, 2.0, 3.2, 2.0, 0.6, 3.6, 3.2, 3.6, 2.6), # match matrix shown in Fig. 1B
            nrow = 3, ncol = 3)

# choose initial conditions and times to show transients and approach to equilibrium
init <- c(0.2, 0.05, 0.3, 0.15, 0.05, 0.25)
times <- seq(from = 0, to = 500, by = 0.1) 

# numerically integrate dynamics
out <- ode(y = init, times = times, func = evaluate_base_model, 
           parms = list(m = m, P = P, n = n), method = "ode45")

tidy_out <- tidy_ode_output(out, n)

p1 <- tidy_out %>% ggplot() + 
  aes(x = Time, y = Frequency, # plot time series for all species and patch types
      group = interaction(species, type),
      color = species, linetype = type) +
  geom_line(size = 1.5) +
  scale_y_continuous(trans = mysqrt) + 
  scale_linetype_manual(values = c("solid", "21"), guide = FALSE) +
  scale_color_manual(values = dark, guide = FALSE) + 
  theme_classic()

show(p1)


##### Figure 2 #####

# 3-panel figure showing some aspects of diversity and coexistence with species-specific memory effects

n <- 4 # number of species 

alphas <- c(-0.9, -1.9, -1.6, -0.4) # species-specific memory effects
beta <- 2 # baseline colonization rate
P <- matrix(beta, n, n) + diag(alphas)

## Panel A

# Time-series showing sequential invasion of species 1-4 

m <- 1 # common local extinction rate (for panels A and B)

init <- c(0.01, 0.99) # initialize first species at low frequency
times <- seq(from = 1, to = 100, by = 0.1) 

# numerically integrate dynamics
out <- ode(y = init, times = times, func = evaluate_base_model, 
           parms = list(m = m, P = P[1, 1], n = 1), method = "ode45")

for(i in 1:(n-1)) { # introduce species 2-4 one at a time
  
  current_state <- out[nrow(out), -1] # get current state
  init <- c(current_state[1:i], 0.01, current_state[(i+1):(i+i)], 0.01) # augment to introduce invader at low frequency
  init <- init / sum(init) # normalize to keep on the simplex
  next_out <- ode(y = init, times = times, func = evaluate_base_model, 
                  parms = list(m = m, P = P[1:(i+1), 1:(i+1)], n = i+1), method = "ode45")
  
  # augment matrix of time-series (pad each time to match new number of species)
  out <- rbind(cbind(out[, 1:(i+1)], rep(0, nrow(out)), # add column of zeroes for invader (previously absent)
                     out[, (i+2):(i+i+1)], rep(0, nrow(out))), # add column of zeroes for invader-state patches (previously absent)
               cbind(times + i * max(times), next_out[, -1])) # add newest dynamics
}

colnames(out) <- c("time", 1:(2*n))
tidy_out <- tidy_ode_output(out, n, vacant = FALSE)

pa <- tidy_out %>% ggplot() + 
  aes(x = Time, y = Frequency, # plot time series showing species only
      group = species, color = species) +
  geom_vline(xintercept = c(0, 100, 200, 300), size = 2, alpha = 0.2) + # indicate invasions with gray lines
  geom_line(size = 1.5) +
  scale_y_sqrt(limits = c(0.01, NA)) + 
  scale_color_manual(values = dark, guide = FALSE) + 
  theme_classic()

## Panel B

# Time-series showing change of m_max (increasing with each invasion)

R_alpha <- cumsum(1 / alphas) # compute sum of reciprocals after each invasion
m_max <- 1 / R_alpha + beta # compute m_max after each invasion
m_df <- data.frame(Time = out[, 1], # generate m_max through time
                   m_max = rep(m_max, each = nrow(out) / n))

pb <- m_df %>% ggplot() + # plot change in m_max vs. time
  aes(x = Time, y = m_max) + 
  geom_vline(xintercept = c(0, 100, 200, 300), size = 2, alpha = 0.2) + # indicate invasions with gray lines
  geom_hline(yintercept = beta, linetype = "dashed") + # indicate beta (asymptotic maximum for m_max) with dashed line
  geom_line(size = 1.5, color = "#69b3a2") + # use a distinct color 
  ylab(expression(m[max])) +
  scale_y_continuous(limits = c(0, NA)) + 
  theme_classic()

## Panel C

# Time-series showing community collapse when species 4 is removed (under m = 1.57)

m <- 1.57 # higher m (so coexistence depends on having all 4 species)

init <- c((1/alphas) * (1/(sum(1/alphas)) + beta - m) / (1 + beta * sum(1/alphas)), # initialize at 4-species equilibrium
          m/alphas * (1 /(1 + beta * (sum(1/alphas)))))
times <- seq(from = 0, to = 100, by = 0.01)

# numerically integrate dynamics 
out <- ode(y = init, times = times, func = evaluate_base_model, 
           parms = list(m = m, P = P, n = n), method = "ode45")

# at t = 100, set species 4 to zero frequency
new_init <- init
new_init[2*n] <- new_init[n] + new_init[2*n]
new_init[n] <- 0

# numerically integrate dynamics 
new_out <- ode(y = new_init, times = max(times) + times, func = evaluate_base_model,
               parms = list(m = m, P = P, n = n), method = "ode45")

out <- rbind(out, new_out) # combine time-series before and after species removal
colnames(out) <- c("time", 1:(2*n))
tidy_out <- tidy_ode_output(out, n, vacant = FALSE)

pc <- tidy_out %>% ggplot() + # plot combined time-series showing species only
  aes(x = Time, y = Frequency, 
      group = species, color = species) +
  geom_vline(xintercept = 100, size = 5, alpha = 0.2) + # show removal with gray line
  geom_line(size = 1.5) +
  scale_y_continuous(trans = mysqrt) + 
  scale_color_manual(values = dark, guide = FALSE) + 
  theme_classic()

# compose panels and plot
p2 <- ggarrange(pa, pb, pc, ncol = 3, nrow = 1)
show(p2)


##### Figure 3 #####

# 4-panel figure illustrating spectral stability condition for symmetric P

set.seed(14)

n <- 10 # number of species
m <- 0.5 # common local extinction rate
times <- seq(0, 500, by = 1)

## Panel A

# Histogram showing eigenvalues for a matrix with stable coexistence equilibrium

P <- sample_P(n, check.stab = TRUE, symmetric = TRUE, scale.diag = 0.1, min = 0.1)$P # find P matrix meeting stability criterion
                                                                                     # scale diagonal elements to speed up search
eigs <- data.frame(Eigenvalues = eigen(P)$values)
pa <- eigs %>% ggplot() + # histogram of eigenvalues
  aes(x = Eigenvalues) +
  geom_histogram(color = "white", 
                 fill = "#69b3a2", # use a distinct color
                 breaks = seq(floor(min(eigs$Eigenvalues)), ceiling(max(eigs$Eigenvalues)), by = 1)) + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) + # indicate zero with a red line
  ylab("Count") + 
  theme_classic()

## Panel B

# Time-series showing approach to equilibrium from a random point

init <- runif(2*n) # random initial condition
init <- init/sum(init)

# numerically integrate dynamics
out <- ode(init, times, evaluate_base_model, parms = list(m = m, P = P, n = n))

tidy_out <- tidy_ode_output(out, n, vacant = FALSE)

pb <- tidy_out %>% ggplot() + # plot time series showing species only
  aes(x = Time, y = Frequency, 
      group = species, color = species) +
  geom_line(size = 1) +
  scale_y_continuous(trans = mysqrt) + 
  scale_color_brewer(palette = "Paired") + 
  theme_classic() + 
  theme(legend.position = "none")

## Panel C

# Histogram showing eigenvalues for a matrix with unstable coexistence equilibrium

P <- sample_P(n, check.stab = FALSE, symmetric = TRUE, scale.diag = 0.1, min = 0.1)$P # find P matrix not meeting stability criterion
                                                                                      # (occurs with high probablity)
                                                                                      # scale diagonal elements to match distribution in A
eigs <- data.frame(Eigenvalues = eigen(P)$values)
pc <- eigs %>% ggplot() + # histogram of eigenvalues
  aes(x = Eigenvalues) +
  geom_histogram(color = "white", 
                 fill = "#69b3a2", # use a distinct color
                 breaks = seq(floor(min(eigs$Eigenvalues)), ceiling(max(eigs$Eigenvalues)), by = 1)) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) + # indicate zero with a red line
  ylab("Count") + 
  theme_classic()

## Panel D

# Time-series showing some species extinctions

init <- runif(2*n) # random initial condition
init <- init/sum(init)

# numerically integrate dynamics
out <- ode(init, times, evaluate_base_model, parms = list(m = m, P = P, n = n))

tidy_out <- tidy_ode_output(out, n, vacant = FALSE)

pd <- tidy_out %>% ggplot() + # plot time series showing species only
  aes(x = Time, y = Frequency, 
      group = species, color = species) +
  geom_line(size = 1) +
  scale_y_continuous(trans = mysqrt) + 
  scale_color_brewer(palette = "Paired") + 
  theme_classic() + 
  theme(legend.position = "none")

# compose panels and plot
p3 <- ggarrange(pa, pb, pc, pd, ncol = 2, nrow = 2)
show(p3)


##### Figure 4 #####

# 2-panel figure showing transition from stable equilibrium, to limit cycles, to extinctions for nonsymmetric P (cyclic and random)

## Panel A

# Progression for cyclic P

n <- 3 # number of species
P <- matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), # cyclic P matrix
            nrow = n, byrow = T)
ms <- c(0.125, 0.25, 0.50) # three values of local extinction rate, to illustrate qualitatively different dynamics

init <- c(0.14, 0.08, 0.007, 0.14, 0.56, 0.064) # initialize near limit cycle for fast convergence
times <- seq(0, 150, by = 1) 

# integrate dynamics with each value of m and collect output
results <- matrix(c(m = ms[1], 0, init), nrow = 1)
for(m in ms) {
  # numerically integrate dynamics
  out <- ode(init, times, evaluate_base_model, parms = list(m = m, P = P, n = n))
  results <- rbind(results, 
                   cbind(m = rep(m, length(times)), out))
}

tidy_results <- tidy_ode_output(results, n, other.factors = "m")
tidy_results$m <- as.factor(tidy_results$m)
levels(tidy_results$m) <- paste("m =", ms) # change levels to get informative facets

pa <- tidy_results %>% ggplot() + 
  aes(x = Time, y = Frequency, # plot time series for all species and patch types at each level of m
      group = interaction(species, type), color = species, linetype = type) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0,1), trans = mysqrt)+ 
  scale_color_manual(values = dark, guide = FALSE) + 
  scale_linetype_manual(values = c("solid", "21")) +
  theme_classic() + 
  theme(legend.title = element_blank(), # position legend to avoid overplotting on trajectories
        legend.position = c(0.9, 0.25),
        legend.margin = margin(-1, 0, -1, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm")) + 
  facet_wrap(.~m) # facet by m

## Panel B

# Progression for a "random" P (particular seed chosen so that m values illustrate distinct outcomes)

set.seed(2)

n <- 3 # number of species
eps <- 0.01

# Find feasible P stable at low m 
P <- sample_P(n, m = "proportion", eps = eps)$P
max_m <- 1 / sum(rowSums(solve(P)))
ms <- round(max_m * c(0.25, 0.72, 1.01), 2) # three values of local extinction rate, to illustrate qualitatively different dynamics

init <- c(0.01, 0.03, 0.3, 0.06, 0.1, 0.5) # initialize near limit cycle for fast convergence
times <- seq(0, 250, by = 1)

# integrate dynamics with each value of m and collect output
results <- matrix(c(m = ms[1], 0, init), nrow = 1)
for(m in ms) {
  # numerically integrate dynamics
  out <- ode(init, times, evaluate_base_model, parms = list(m = m, P = P, n = n))
  results <- rbind(results, 
                   cbind(m = rep(m, length(times)), out))
}

tidy_results <- tidy_ode_output(results, n, other.factors = "m")
tidy_results$m <- as.factor(tidy_results$m)
levels(tidy_results$m) <- paste("m =", ms) # change levels to get informative facets

pb <- tidy_results %>% ggplot() + # plot time series for all species and patch types at each level of m
  aes(x = Time, y = Frequency, 
      group = interaction(species, type), color = species, linetype = type) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0,1), trans = mysqrt)+ 
  scale_color_manual(values = dark, guide = FALSE) + 
  scale_linetype_manual(values = c("solid", "21"), guide = FALSE) +
  theme_classic() + 
  facet_wrap(.~m) # facet by m

# compose panels and plot
p4 <- ggarrange(pa, pb, nrow = 2, ncol = 1)
show(p4)


##### Figure 5 #####

# 3-panel figure showing parameter space for waning memory model, and illustrating bistability 
# with two simulations from different initial conditions

## Panel A

# Long-term outcomes (assuming stable, symmetric P and P^{-1} 1 > 0 elementwise) across parameter space

m_vals <- c(0.25, 0.5, 1) # values of common local extinction rate
d_vals <- seq(0.01, 1, length.out = 300) # values of common memory decay rate
c_vals <- seq(0.01, 1.5, length.out = 300) # values of common colonization rate for naive patches
q_vals <- c(0.25, 0.5, 1) # values of q = 1^T P^{-1} 1

# pre-allocate matrix for results at all parameter combinations
results <- matrix(nrow = length(m_vals)*length(d_vals)*length(c_vals)*length(q_vals), 
                  ncol = 9)

# loop over all parameter combinations and classify long-term dynamics
i <- 0 # keep a counter to i.d. each combination
for(m in m_vals) {
  for(d in d_vals) {
    for(c in c_vals) {
      for(q in q_vals) {
        
        i <- i + 1
        
        # compute k values in two parts
        a <- d/m + 1/(m * q) - 1 # shared part of roots
        b <- sqrt((1 - d/m - 1/(m*q))^2 - 4 * (d/m) * (1/(c*q) - 1)) # discriminant part
        
        k1 <- (a + b) / 2
        k2 <- (a - b) / 2
        
        feas_thresh <- d / m # threshold for feasiblity (k > d/m)
        stable1 <- if(m / (d + m) < c * q) TRUE else k1 > sqrt((d * (1 - c * q)) / (m * c * q)) # check if equilibrium associated with k1 
                                                                                                # is stable
        stable2 <- if(m / (d + m) < c * q) TRUE else k2 > sqrt((d * (1 - c * q)) / (m * c * q)) # check if equilibrium associated with k2
                                                                                                # is stable (this should never happen)
        results[i, ] <- c(m, d, c, q, k1, k2, feas_thresh, stable1, stable2)
      }
    }
  }
}

results <- results %>% as_tibble()
colnames(results) <- c("m", "d", "c", "q", "k1", "k2", "feas_thresh", "stable1", "stable2")

results$m <- as.factor(results$m)
levels(results$m) <- paste("m =", m_vals) # change levels to get informative facets
results$m <- factor(results$m, levels(results$m)[3:1]) # re-order to get a nicer plot layout

results$q <- as.factor(results$q)
levels(results$q) <- paste("q =", q_vals) # change levels to get informative facets

pa <- results %>% # classify outcome at each parameter combination
  mutate(k1 = ifelse(is.nan(k1), -1, k1), # if k = nan, no coexistence, so set k < 0
         k2 = ifelse(is.nan(k2), -1, k2),
         status = ifelse(k1 > feas_thresh, # k1 feasible
                         ifelse(stable1, # k1 stable
                                ifelse(k2 > feas_thresh, # k2 feasible
                                       ifelse(stable2, 
                                              "two stable", # if k2 also stable (this should never happen)
                                              "bistable"), # if both feasible but only k1 stable
                                       "unique stable"),  # if only k1 feasible (and stable)
                                "unstable"),  # if k1 unstable (k2 must be unstable also)
                         "not feasible")) %>% # if k1 not feasible (k2 must be unfeasible also)
  ggplot() + 
  aes(x = d, y = c) + # plot outcomes across parameter spaces
  geom_tile(aes(fill = status)) + 
  geom_point(data = data.frame(m = "m = 0.5", q = "q = 0.5", d = 0.1, c = 0.25), fill = NA, pch = 5) + # indicate parameter
                                                                                                       # used for panels B and C
  geom_text(data = data.frame(m = "m = 1", q = "q = 0.25", d = 0.25, c = 0.75), aes(label = "bistable")) + 
  geom_text(data = data.frame(m = "m = 1", q = "q = 0.25", d = 0.75, c = 0.25), aes(label = "unfeasible")) + # label colors 
  geom_text(data = data.frame(m = "m = 1", q = "q = 0.25", d = 0.5, c = 1.25), aes(label = "unique stable")) + 
  facet_grid(m ~ q) + 
  theme_classic() + 
  scale_fill_manual(values = c( "#56B4E9", "#D55E00", "#0072B2", "#CC79A7")) + 
  scale_x_continuous(labels = function(x) round(x, 2)) + # round to get cleaner axes
  theme(legend.position = "none")

## Panel B 

# Time series showing coexistence when naive patches not too abundant initially

set.seed(82)

# parameters for B AND C (only initial conditions differ)
n <- 3 # number of species
m <- 0.5 # common local extinction rate
q <- 0.5 # q = 1^T P^{-1} 1 (to re-scale P)
d <- 0.1 # common memory decay rate
c <- 0.25 # common colonization rate of naive patches

P <- matrix(c(2.2, 2.0, 3.2, 2.0, 0.6, 3.6, 3.2, 3.6, 2.6), # match Fig. 1B
            nrow = 3, ncol = 3)
P <- (sum(solve(P)) / q) * P # re-scale P to have desired q

init <- c(runif(2*n), 2) # initialize at a random point with z not too large
init <- init/sum(init)
times <- seq(from = 0, to = 50, by = 0.1) 

# numerically integrate dynamics
out <- ode(y = init, times = times, func = evaluate_waning_model, 
           parms = list(m = m, P = P, n = n, c = c, d = d), method = "ode45")

tidy_out <- tidy_ode_output(out, n)

pb <- tidy_out %>% ggplot() + 
  aes(x = Time, y = Frequency, # plot time series showing all species and patch types
      group = interaction(species, type), color = species, linetype = type) +
  geom_line(lwd = 1.5) + 
  scale_y_continuous(trans = mysqrt)+ 
  scale_color_manual(values = dark) + 
  theme_classic() + 
  scale_linetype_manual(values = c("solid", "21", "11")) + 
  theme(legend.position = "none")

## panel c

init <- c(runif(2*n), 10) # initialize at a random point with z large
init <- init/sum(init)
times <- seq(from = 0, to = 50, by = 0.1) 

# numerically integrate dynamics
out <- ode(y = init, times = times, func = evaluate_waning_model, 
           parms = list(m = m, P = P, n = n, c = c, d = d), method = "ode45")

tidy_out <- tidy_ode_output(out, n)

pc <- tidy_out %>% ggplot() + 
  aes(x = Time, y = Frequency, # plot time series with all species and patch types
      group = interaction(species, type), color = species, linetype = type) +
  geom_line(lwd = 1.5) +
  scale_y_continuous(trans = mysqrt)+ 
  guides(linetype = guide_legend(title = element_blank())) + 
  scale_color_manual(values = dark, guide = "none") + 
  theme_classic() + 
  scale_linetype_manual(values = c("solid", "21", "11")) +
  theme(legend.position = c(0.85, 0.55), # position legend to avoid overplotting on trajectories
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(0, "mm"))

# compose panels and plot
p5 <- ggarrange(pa, 
          ggarrange(pb, pc, nrow = 2), 
          ncol = 2, widths = c(2,1))
show(p5)