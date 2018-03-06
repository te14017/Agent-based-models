
### SIS and SIR model simulation
### Copyright: Te Tan, University of Zurich, March 2018

# model parameters
SIS <- FALSE # whether it runs as a SIS model or SIR
p_beta <- 0.8 # infection rate
p_beta_list <- sort(seq(0.2, 1, by = 0.1), decreasing = FALSE) # try different infection rates
p_gamma <- 1 # recovery rate
n <- 5000 # system size
n0 <- 500 # initial size infected
iter <- 10000 # number of iterations

# initialize the population with size n, and n0 individuals infected
init_population <- function(n, n0) {
  population <- rep(0, n)
  # initial infected individual's indexes
  infect_index <- sample(seq(1, n), n0)
  for(i in seq(1, n0)) {
    population[infect_index[i]] <- 1
  }
  return(population)
}

# compute which case will happen during one iteration: infection or recovery
# return: TRUE if infection will happen, FALSE if recovery will
if_infect_happen <- function(p_beta, p_gamma, papulation) {
  nr_of_suscept <- length(population[population==0])
  n <- length(population)
  prob_infect <- (p_beta*(nr_of_suscept/n))/(p_beta*(nr_of_suscept/n) + p_gamma)
  random <- runif(1)[1]
  if (random <= prob_infect) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# get the index of the randomely infected individual
infect_index <- function(population) {
  index_of_suscept <- which(population == 0)
  if (length(index_of_suscept) == 0) {
    return(0)
  } else {
    index <- sample(index_of_suscept, 1)
    return(index)
  }
}

# get the index of the randomly recovered individual
recover_index <- function(population) {
  index_of_infect <- which(population == 1)
  if (length(index_of_infect) == 0) {
    return(0)
  } else {
    index <- sample(index_of_infect, 1)
    return(index)
  }
}

# run function with the the list of beta values
run_with_beta <- function(SIS) {
  nr_of_infected <- vector() # vector to gather number of infected people of each iteration
  nr_of_recover_stationary <- vector()
  nr_of_infect_stationary <- vector()
  
  population <- init_population(n, n0)
  size <- length(population)
  j <- 0 # number of actual iterations run
    
  for (k in 1:length(p_beta_list)) { # iteration for computing every beta
    for (i in 1:iter) { # iterations of evolution
      if (if_infect_happen(p_beta_list[k], p_gamma, population)) {
        index <- infect_index(population)
        if (index == 0) {
          #print("Program terminate: all people get infected.")
          #break
        }
        population[index] <- 1
      } else {
        index <- recover_index(population)
        if (index == 0) {
          #print("Program terminate: all people get recovered.")
          #break
        }
        if (SIS) {
          population[index] <- 0 # SIS model
        } else {
          population[index] <- 2 # SIR model
        }
      }
      #j <- j + 1
      #nr_of_infected[i] <- sum(population[population==1])/size
    }
    #plot(1:j, nr_of_infected, type = "l",main = "Ratio of infected people over iterations",
      #sub = "with beta equals 0.8", xlab = "Number of iterations", ylab = "Ratio of infected people")
    
    nr_of_infect_stationary[k] <- sum(population[population==1])/size
    if (!SIS) {
      nr_of_recover_stationary[k] <- length(population[population==2])/size
    }
  }
  # plot figures
  if (SIS) {
  plot(p_beta_list, nr_of_infect_stationary, type = "l",main = "Ratio of infected in stationary state over beta",
      xlab = "beta/gamma with gamma as 1", ylab = "Ratio of infected people")
  } else {
    plot(p_beta_list, nr_of_recover_stationary, type = "l",main = "Ratio of recovered in stationary state over beta(SIR)",
         xlab = "beta/gamma with gamma as 1", ylab = "Ratio of recovered people")
  }
}

# running body
run_with_beta(SIS)