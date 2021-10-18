# Functions to simulate, analyze, and plot model dynamics 
# zachmiller@uchicago.edu

library("scales")

THRESH <- 10^-6 # global threshold frequency for extinction

##### Functions for model dynamics #####

evaluate_base_model <- function(t, vars, pars) {
  
  # Calculate derivatives to numerically integrate basic model (Eq. 1 in Main Text) using deSolve::ode()
  #
  # Input
  # t : current time (numeric)
  # vars : numeric vector of length 2n; first n values are interpreted as x_i, 
  #        second n values are interpreted y_i in corresponding order 
  # pars : list of parameters (see below for details)
  #
  # Output
  # list of derivatives
  
  m <- pars$m # local extinction rates (either a numeric vector of length n, or a single value for all species)
  P <- pars$P # colonization rates (numeric n x n matrix)
  n <- pars$n # number of species
  
  x <- vars[1:n]
  y <- vars[(n + 1):(n + n)]
  
  # check if any species have gone extinct (i.e. set species below THRESH to 0)
  if (any(x < THRESH)) {
    ids <- which(x < THRESH)
    x[ids] <- 0
  }
  
  # normalize frequencies to keep dynamics on the simplex
  if (sum(x + y) > 0) {
    x <- x / sum(x + y)
    y <- y / sum(x + y)
  }
  
  # compute derivatives according to Eq 1.
  dxdt <- x * (-m + P %*% y)
  dydt <- x * m - y * t(P) %*% x
  
  return(list(c(dxdt, dydt)))
}

evaluate_waning_model <- function(t, vars, pars) {
  
  # Calculate derivatives to numerically integrate waning memory model (Eq. 58 in SI) using deSolve::ode()
  #
  # Input
  # t : current time (numeric)
  # vars : numeric vector of length 2n + 1; first n values are interpreted as x_i, 
  #        second n values are interpreted y_i in corresponding order, 
  #        final value is interpreted as z, the frequency of naive patches
  # pars : list of parameters (see below for details)
  #
  # Output
  # list of derivatives
  
  m <- pars$m # local extinction rates (either a numeric vector of length n, or a single value for all species)
  P <- pars$P # colonization rates (numeric n x n matrix)
  n <- pars$n # number of species 
  d <- pars$d # memory decay rates (either a numeric vector of length n, or a single value for all species)
  c <- pars$c # colonization rates of naive patches (either a numeric vector of length n, or a single value for all species)

  x <- vars[1:n]
  y <- vars[(n + 1):(n + n)]
  z <- vars[2 * n + 1] 
  
  # check if any species have gone extinct (i.e. set species below THRESH to 0)
  if (any(x < THRESH)) {
    ids <- which(x < THRESH)
    x[ids] <- 0
  }
  
  # normalize frequencies to keep dynamics on the simplex
  if (sum(x + y) + z > 0) {
    x <- x / (sum(x + y) + z)
    y <- y / (sum(x + y) + z)
    z <- z / (sum(x + y) + z)
  }
  
  # compute derivatives according to Eq. SI 58.
  dxdt <- x * (-m + c * z + P %*% y)
  dydt <- x * m - y * t(P) %*% x - d * y
  dzdt <- d * sum(y) - c * z * sum(x) 
  
  return(list(c(dxdt, dydt, dzdt)))
}


##### Functions for plotting #####

tidy_ode_output <- function(out, n, vacant = TRUE, other.factors = NULL) {
  
  # Wrangle time-series data output by deSolve::ode into a convenient form for ggplot
  #
  # Input
  # out : matrix of class deSolve, from integrating evalute_base_model or evaluate_waning_model
  # n : number of species
  # vacant : should vacant patch frequencies be retained? (logical)
  # other.factors : names of any additional columns in out (e.g. classifying multiple simulation runs)
  #
  # Output
  # data in tidy form
  
  out <- out %>% as.data.frame() %>% 
    gather(key = species, value = Frequency, -c(time, other.factors)) %>% # reshape data to long form
    rename(Time = time) %>%
    mutate(species = as.numeric(species), # temporarily convert to numeric for next step
           type = ifelse(species <= n, "x", # classify type (species, vacant patches, naive patches) according to species number
                  ifelse(species == 2 * n + 1, "z", "y")),
           species = ifelse(species <= n, species, # re-number species to x_i matches y_i
                            ifelse(species == 2 * n + 1, "z", species - n)),
           species = as.factor(species), # return species label to factor
           Frequency = ifelse(Frequency < THRESH, 0, Frequency)) # double-check very small frequencies set to 0
  
  if (!vacant) out <- out %>% filter(type == "x") # discard vacant patches if not needed for plotting
  
  return(out)
}

# custom square root transformation for ggplot scales (see discussion: https://github.com/tidyverse/ggplot2/issues/980)
mysqrt <- trans_new("mysqrt", 
                    transform = sqrt,
                    inverse = function(x) ifelse(x < 0, 0, x^2),
                    domain = c(0, Inf))


##### Functions for analysis #####

sample_P <- function(n, check.feas = TRUE, check.stab = TRUE, symmetric = FALSE, 
                     scale.diag = 1, m = 1, eps = 1, distribution = runif, ...) {
  
  # Generate a P matrix at random, and possibly according to some constraints
  # 
  # Input
  # n : number of species
  # check.feas : should the output satisfy P^{-1} m > 0 ? (keep sampling until one is found)
  # check.stab : should the output have a stable coexistence equilibrium? (keep sampling until one is found)
  # symmetric : should the output be symmetric?
  # scale.diag : a numeric value by which diagonal elements are re-scaled
  # m : one of (i) a vector of length n, (ii) a single value to be used for all species, or
  #     (iii) the string "proportional", indicating that a fraction of m_max (given by eps) should 
  #     be used for all species
  # eps : if m = "proportional", the fraction of m_max to be used for m
  # distribution : a function to generate random deviates, from which elements of P are sampled i.i.d.
  #     (subject to symmetry or re-scaling of diagonal)
  # ... : any additional parameters for distribution function 
  #
  # Output
  # an n x n matrix P satisfying input constraints, and logicals feas (does the output have feasible equilibrium)
  #     and stab (does the output have stable equilibrium)
  
  success <- FALSE
  while(!success) { # sample until matrix with desired properties is found
    
    stab <- FALSE
    feas <- FALSE
    
    P <- matrix(distribution(n = n*n, ...), n, n) # draw entries i.i.d. from distribution function
    if (symmetric) P <- P + t(P) # enforce symmetry if required
    diag(P) <- scale.diag * diag(P) # re-scale diagonal (no re-scaling by default)
    
    # if m = "proportion", compute max_m and then m as eps * max_m
    if (m == "proportion") {
      max_m <- 1 / sum(rowSums(solve(P)))
      m <- eps * max_m
    }
    
    # if a single value of m is specified, repeat to obtain a vector of length n
    if (length(m) == 1) m <- rep(m, n)
    
    # check if P and m yield feasible equilibrium
    if (all(solve(P, m) > 0)) feas <- TRUE
    
    # check if P and m yield stable equilibrium
    if (symmetric) { # for symmetric P, just check eigenvalue condition
      eP <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
      if (sum(Re(eP) > 0) == 1) stab <- TRUE # stable iff exactly one positive eigenvalue
    } else { 
      # for nonsymmetric P, build the Jacobian at equilibrium and check eigenvalues
      J <- build_jacobian(P, m)$J
      if (all(Re(eigen(J, only.values = TRUE)$val)[-2*n] < 0)) stab <- TRUE # exclude the (structural) zero eigenvalue
    }
    
    # check that required properties are satisfied; if not, discard and try again
    if ((!check.feas | feas) & (!check.stab | stab)) success <- TRUE
  }
  
  return(list(P = P, feas = feas, stab = stab))
}


build_jacobian <- function(P, m, symmetric = FALSE) {
  
  # Construct the Jacobian matrix evaluated at coexistence equilibrium
  #
  # Input
  # P : n x n numeric matrix of colonization rates
  # m : local extinction rates (vector of length n)
  # symmetric : is the matrix symmetric? (if so, use proportionality to compute xstar)
  #
  # Output
  # Jacobian matrix, J, and equilibrium frequencies xstar and ystar
  
  # get number of species
  n <- nrow(P)
  
  # compute ystar
  ystar <- solve(P, m)
  
  # compute xstar
  if (symmetric) { # if P is symmetric, compute k and use xstar = k * ystar
    k <- 1/sum(ystar) - 1
    xstar <- k * ystar
  } else { # if P is nonsymmetric, compute xstar as an eigenvector
    xstar <- Re(eigen(diag(m) - diag(ystar) %*% t(P))$vectors[,n])
    xstar <- (xstar / sum(xstar)) * (1 - sum(ystar))
  }
  
  # construct J according to Eq. 30 in SI
  J <- matrix(0, 2*n, 2*n)
  J[1:n, (n+1):(2*n)] <- diag(xstar) %*% P
  J[(n+1):(2*n), 1:n] <- diag(m) - diag(ystar) %*% t(P)
  J[(n+1):(2*n), (n+1):(2*n)] <- -diag(as.vector(t(P) %*% xstar))
  
  return(list(J = J, xstar = xstar, ystar = ystar))
} 


simulate_and_classify_dynamics <- function(P, m, start = "nearby", max.time = 5000, max.rounds = 10) {
  
  # Classify the long-term dynamics as (full) equilibrium (based on local stability analysis),  
  # extinction of one or more species, or limit cycles (the latter two distinguished by integrating the dynamics) 
  #
  # Input
  # P : n x n numeric matrix of colonization rates
  # m : local extinction rates (either a numeric vector of length n, or a single value for all species)
  # start : either a vector of initial frequencies (length 2n), or the string "nearby", to initalize
  #         dynamics at a random point near the equilibrium frequencies
  # max.time : length of time to integrate for each round
  # max.rounds : Total number of times the dynamics will be integrated. After each round, check if some species 
  #              have gone extinct. If none have after max.rounds iterations, check if the outcome is a limit 
  #              cycle. If not, the result is indeterminate.
  #
  # Output
  # A flag indicating the type of dynamics: (i) "equilibrium", (ii) "extinction", (iii) "cycle", or (iv) "indeterminate"
  
  # get number of species
  n <- nrow(P)
  
  # build Jacobian and get coexistence equilibrium frequencies 
  J_and_eq <- build_jacobian(P, m)
  J <- J_and_eq$J
  eq <- c(J_and_eq$xstar, J_and_eq$ystar)
  
  # check if equilibrium is locally stable
  eJ <- eigen(J, only.values = TRUE)$val
  if (all(Re(eJ[-2*n]) <= 0)) flag <- "equilibrium" else flag <- "continue"
  
  # prepare to numerically integrate
  pars <- list(n = n, P = P, m = m)
  tt <- seq(0, max.time, by = 1)
  
  # if "nearby" is chosen, initialize frequencies at a random point near equilibrium
  if (start == "nearby") {
    z0 <- eq * runif(2 * n, 0.9, 1.1) # 10% noise around equilibrium
    z0 <- z0 / sum(z0)
  } else { # otherwise, use supplied initial conditions
    z0 <- start
  }
  
  # keep integrating dynamics until an extinction is observed or max.rounds is reached (skip if equilibrium was stable)
  rounds <- 1
  while(flag == "continue") {
    
    # numerically integrate model
    dynamics <- ode(y = z0, times = tt, func = evaluate_base_model, parms = pars, method = "ode45")
    
    # get final frequencies and check for extinctions
    zfinal <- as.numeric(dynamics[nrow(dynamics), -1])
    if (any(zfinal < THRESH)) {
      flag <- "extinction"
    } else { # if max.rounds has been reached, check if the dynamics have approached a limit cycle
      if (rounds >= max.rounds) {
        if (detect_limit_cycle(dynamics[, -1])) {
          flag <- "cycle"
        } else { # if outcome cannot be classified, return "indeterminate" 
          flag <- "indeterminate"
        }
      } else {
        z0 <- zfinal # use final frequencies from last found as initial conditions for next
        rounds <- rounds + 1 
      }
    }
  }
  
  return(flag)
}


detect_limit_cycle <- function(dyn, scale = 10^-2, times = 5) {
  
  
  # Check if dynamics are cycling 
  #
  # Input
  # dyn : ode output with "time" variable discarded (i.e. a matrix where each col is a variable and each row is a time point)
  # scale : a numeric value indicating the fraction of the mean species' standard deviation that should be used as a tolerance
  #         (distance threshold) to call two points identical
  # times : how many times must the last point be repeated in the time series to call the dynamics a limit cycle
  #
  # Output
  # a logical -- are the dynamics cycling stably?
  
  
  tol <- mean(apply(dyn, 2, sd)) * scale # compute a distance threshold for calling two points the same 
                                         # (proportional the mean amplitude of any cycles)
  
  last_val <- dyn[nrow(dyn), ]
  same_pts <- apply(dyn, 1, function(x) sum((x - last_val)^2) < tol) + 0 # classify every point as near the last observed, or not
  same_pts <- diff(same_pts) # take the difference to isolate separate visits to the last point 
  out <- sum(same_pts == 1) > times # is the number of separate visits to the last point greater than times argument?
  
  return(out)
}


build_cyclic_P <- function(n, c) {
  
  # Construct a cyclic permutation matrix of size n, scaled by a value c
  #
  # Input
  # n : number of species
  # c : common colonization rate for all species (non-zero value in P)
  #
  # Output
  # an n x n matrix, P
  
  P <- matrix(0, nrow = n, ncol = n)
  P[1, n] <- 1
  P[2:n, 1:(n-1)] <- diag(n-1)
  P <- c * P
  
  return(P)
}

custom_sqrt <- function(x){
  
  # Evaluate sqrt if input is positive, return 0 otherwise (for computing real parts of eigenvalues)
  
  if (x > 0){
    out <- sqrt(x)
  }else{
    out <- 0
  }
  
  return(out)
}
