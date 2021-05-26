# Run supplemental simulations for Miller and Allesina, "Metapopulations with habitat modification"
# zachmiller@uchicago.edu

source("functions.R")

##### Explore variation in m #####

# Sample many symmetric P matrices meeting (or not) the spectral stability condition. For each, sample
# many random m vectors (with m_i != m_j) and check local stability. Stop and print "bad instance" if
# an example is found where stability type (i.e. stability vs. instability) with varied m does not match 
# stability type with constant m

set.seed(53)

n <- 4 # number of species
num_mats <- 1000 # how many matrices to sample for each stability type
num_reps <- 1000 # how many random m vectors to try for each matrix

# keep count of stable (i) and unstable (j) matrices, and keep sampling until we find num_mats of each type 
i <- 0 
j <- 0
while(i < num_mats | j < num_mats){
  
  print(paste(i, "of", num_mats, "stable matrices and", j, "of", num_mats, "unstable matrices checked."))
  
  success <- FALSE
  while(!success){ # keep the next sampled P if its type is still need
    
    out <- sample_P(n, check.feas = TRUE, check.stab = FALSE, symmetric = TRUE)
    type <- out$stab
    P <- out$P
    
    if(type){ # if stable (assuming constant m)
      if(i < num_mats){
        success <- TRUE
        i <- i + 1
      }
    }else{ # else unstable (assuming constant m)
      if(j < num_mats){
        success <- TRUE
        j <- j + 1
      }
    }
  }
  
  # loop over many random m vectors
  for(nm in 1:num_reps){
    
    # to find m vectors compatible with feasibility, sample random equilibrium (i.e. vector proportional to ystar -- see SI text
    # for details), then compute corresponding m
    
    y <- rexp(n) # sample y (up to rescaling) uniformly from the simplex
    y <- y / sum(y)
    m <- P %*% y # corresponding m
    
    # get m_max and then rescale by a random fraction
    max_scale <- 1 / sum(solve(P) %*% m)
    m <- runif(1) * max_scale * m
    
    # build Jacobian at equilibrium
    J <- build_jacobian(P, as.vector(m), symmetric = TRUE)$J
    
    # check stability
    stable <- all(Re(eigen(J, only.values = TRUE)$val)[-2*n] < 0) # ignore (structural) 0 eigenvalue
    if(type != stable){
      print("bad instance") # if stability type does not match the type without variation in m, break
      
      # to break out of while loop
      i <- Inf
      j <- Inf
      
      break
    }
  }
}


##### Explore variation in c #####

# Sample many symmetric P matrices meeting (or not) the spectral stability condition. For each, sample
# many random c vectors (i.e. species-specific dispersal abilities) and check local stability. Stop and print 
# "bad instance" if an example is found where stability type (i.e. stability vs. instability) with varied c 
# does not match stability type with equal dispersal abilities

set.seed(88)

n <- 4 # number of species
num_mats <- 1000 # how many matrices to sample for each stability type
num_reps <- 1000 # how many random m vectors to try for each matrix

# keep count of stable (i) and unstable (j) matrices, and keep sampling until we find num_mats of each type 
i <- 0
j <- 0
while(i < num_mats | j < num_mats){
  
  print(paste(i, "of", num_mats, "stable matrices and", j, "of", num_mats, "unstable matrices checked."))
  
  success <- FALSE
  while(!success){ # keep the next sampled P if its type is still need
    
    out <- sample_P(n, check.feas = TRUE, check.stab = FALSE, symmetric = TRUE)
    type <- out$stab
    P <- out$P
    
    if(type){ # if stable (assuming constant m)
      if(i < num_mats){
        success <- TRUE
        i <- i + 1
      }
    }else{ # else unstable (assuming constant m)
      if(j < num_mats){
        success <- TRUE
        j <- j + 1
      }
    }
  }
  
  # loop over many random c vectors
  for(nm in 1:num_reps){
    
    # to find c vectors compatible with feasibility, sample random equilibrium (i.e. vector proportional to ystar -- see SI text
    # for details), then compute corresponding c
    
    y <- rexp(n) # sample y (up to rescaling) uniformly from the simplex
    y <- y / sum(y)
    c <- diag(1 / as.vector(P %*% y)) # corresponding c
    
    P <- c %*% P # update P
    m <- runif(1) # random m
    
    # build Jacobian at equilibrium
    J <- build_jacobian(P, rep(m, n), symmetric = FALSE)$J
    
    stable <- all(Re(eigen(J, only.values = TRUE)$val)[-2*n] < 0) # ignore (structural) 0 eigenvalue
    if(type != stable){
      print("bad instance") # if stability type does not match the type without variation in m, break
      
      # to break out of while loop
      i <- Inf
      j <- Inf
      
      break
    }
  }
}


##### Explore global stability #####

# Sample many symmetric P matrices meeting the spectral stability condition, and integrate the model dynamics from
# many random initial conditions. Check if the dynamics always converge to the coexistence equilibrium. Stop and
# print "bad instance" if an example is found where a trajectory does not approach the equilibrium.

set.seed(97)

n <- 4 # number of species
tol <- 10^-4 # distance threshold to determine convergence to equilibrium
max.rounds <- 20 # how many times to keep integrating if equilibrium or extinctions not detected

n_mats <- 1000 # how many matrices to sample
n_ics  <- 100 # how many initial conditions to try for each matrix

for(i in 1:n_mats){
  
  print(paste(i, "of", n_mats, "matrices checked."))
  
  sample_P(n, symmetric = TRUE) # sample symmetric P feasible and stable
  
  # use m = 0.5 * max_m for all species
  max_m <- 1 / sum(rowSums(solve(P)))
  m <- max_m * 0.5
  
  # compute equilibrium frequencies
  ystar <- m * rowSums(solve(P))
  xstar <- ystar * (1 / sum(ystar) - 1)
  zstar <- c(xstar, ystar)
  
  # prepare to integrate
  pars <- list(n = n, P = P, m = m)
  tt <- seq(0, 5000, by = 1)
  
  # loop over random initial conditions
  for(j in 1:n_ics){
    
    # sample initial conditions uniformly on the simplex
    z0 <- rexp(2*n)
    z0 <- z0 / sum(z0)
    
    rounds <- 1
    flag <- "continue"
    while(flag == "continue"){ # keep integrating until dynamics can be classified or max.rounds is reached
      
      # numerically integrate dynamics
      dynamics <- ode(y = z0, times = tt, func = evaluate_base_model, parms = pars, method = "ode45")
      
      # get final frequencies and check for equilibrium and extinctions
      zfinal <- as.numeric(dynamics[nrow(dynamics), -1])
      if(sum((zfinal - zstar)^2) < tol){ # if final frequencies are close enough to equilibrium, call convergence
        flag <- "equilibrium"
      } else {
        if (any(zfinal < THRESH)){
          flag <- "extinction"
        } else { # if max.rounds has been reached, check if the dynamics have approached a limit cycle
          if (rounds >= max.rounds){
            if (detect_limit_cycle(dynamics[, -1])){
              flag <- "cycle"
            } else {
              flag <- "indeterminate"
            }
          } else{
            z0 <- zfinal # use final frequencies from last found as initial conditions for next
            rounds <- rounds + 1 
          }
        }
      }
    }
    # if any outcome except convergence to equilibrium, stop and print
    if(flag != "equilibrium"){
      print("bad instance")
      break
    }
  }
}