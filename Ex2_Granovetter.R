
### Granovetter model simulation for riot problem
### Copyright (C): Te Tan, University of Zurich, March 2018

# model parameters
n <- 5000 # size of the population
n0 <- 100 # initial number of people joining riot
p_mu <- 0.2 # mean value for Gaussian distribution
p_sigma_list <- sort(seq(0.01, 1.0, by = 0.02), decreasing = FALSE) # try different variance

# initialize population
init_population <- function(n, n0) {
  population <- rep(0, n)
  # initial infected individual's indexes
  riot_index <- sample(seq(1, n), n0)
  for(i in seq(1, n0)) {
    population[riot_index[i]] <- 1
  }
  return(population)
}

# initialize threshold vector for individuals
init_threshold <- function(n, p_mu, p_sigma) {
  threshold <- rnorm(n, p_mu, p_sigma)
  return(threshold)
}

# run function
run_granovetter <- function() {
  nr_of_riot_final <- vector() # gather the number of people joining riot finally
  
  for (i in 1:length(p_sigma_list)) {
    population <- init_population(n, n0)
    threshold <- init_threshold(n, p_mu, p_sigma_list[i])
    
    # compute current ratio of people joining riot
    riot_ratio <- n0/n
    previous_ratio <- 0
    while (TRUE) {
      previous_ratio <- riot_ratio
      # get the index of people who will join the riot
      index_join_riot <- which(threshold<=riot_ratio)
      population[index_join_riot] <- 1
      # update ratio of people joining riot
      riot_ratio <- length(population[population==1])/n
      # when no more people join the riot, loop terminates
      if (previous_ratio == riot_ratio) {
        break
      }
    }
    nr_of_riot_final[i] <- sum(population)/n
  }
  plot(p_sigma_list, nr_of_riot_final, type = "l", main = "Ratio of people joint riot finally",
       xlab = "sigma of threshold Gaussian distribution", ylab = "Ratio of people joint riot")
}


# run body
run_granovetter()