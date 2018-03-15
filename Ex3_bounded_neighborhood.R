
### Bounded-neighbourhood model simulation (RED and BLUE population)
### Copyright (C): Te Tan, University of Zurich, March 2018

# model parameters
n <- 10000 # population size
m <- 1000 # neighborhood size
r <- 0.5 # initial rate of RED people, represented as 1. BLUE as 0
Rr <- 2 # tolerance for RED people
Rb <- sort(seq(0, 4.0, by = 0.1), decreasing = FALSE) # tolerance rate for BLUE people

# initialize population
init_population <- function(n, r) {
  # initialize all people as BLUE(0)
  population <- rep(0, n)
  # initial infected individual's indexes
  red_index <- sample(seq(1, n), r*n)
  for(i in seq(1, r*n)) {
    population[red_index[i]] <- 1
  }
  return(population)
}

# initialize an empty neighborhood vector: 1 is in neighborhood, 0 is not
init_neighborhood <- function(n) {
  return(rep(0, n))
}

# whether neighborhood is full
is_full <- function(neigh, m) {
  return(sum(neigh) >= m)
}

# return current ratio of RED in neighborhood
red_ratio_in_neigh <- function(popu, neigh, m) {
  if (length(which(neigh==1)) == 0) return(0)
  else {
    in_neigh_index <- which(neigh==1)
    return(sum(popu[in_neigh_index])/m)
  }
}

# return current ratio of RED in neighborhood
blue_ratio_in_neigh <- function(popu, neigh, m) {
  if (length(which(neigh==1)) == 0) return(0)
  else {
    in_neigh_index <- which(neigh==1)
    return((length(in_neigh_index)-sum(popu[in_neigh_index]))/m)
  }
}

# function returns the tolerance
tolerance <- function(R, in_ratio) {
  # slove the line equation to get coefficient
  a <- solve(100, -R)
  return((a*in_ratio+R)*in_ratio)
}

# run function
run_bounded_neigh <- function() {
  red_ratio_equil <- vector()
  blue_ratio_equil <- vector()
  
  for (i in 1:length(Rb)) {
    popu <- init_population(n, r)
    neigh <- init_neighborhood(n)
    
    red_in <- red_ratio_in_neigh(popu, neigh, m)
    blue_in <- blue_ratio_in_neigh(popu, neigh, m)
    previous_red_in <- 0
    previous_blue_in <- 0
    
    while (TRUE) {
      previous_red_in <- red_in
      previous_blue_in <- blue_in
      
      for (j in 1:n) {
        # BLUE in or out according to current ratio of BLUE
        if (popu[j]==0 & neigh[j]==0 & red_in <= tolerance(Rb[i], blue_in) & !is_full(neigh, m)) {
          neigh[j] <- 1
          next
        }
        if (popu[j]==0 & neigh[j]==1 & red_in > tolerance(Rb[i], blue_in)) {
          neigh[j] <- 0
          next
        }
        # RED in or out according to current ratio of BLUE
        if (popu[j]==1 & neigh[j]==0 & blue_in <= tolerance(Rr, red_in) & !is_full(neigh, m)) {
          neigh[j] <- 1
          next
        }
        if (popu[j]==1 & neigh[j]==1 & blue_in > tolerance(Rr, red_in)) {
          neigh[j] <- 0
          next
        }
      }
      red_in <- red_ratio_in_neigh(popu, neigh, m)
      blue_in <- blue_ratio_in_neigh(popu, neigh, m)
      if (red_in == previous_red_in & blue_in == previous_blue_in) break
    }
    
    red_ratio_equil[i] <- red_in
    blue_ratio_equil[i] <- blue_in
  }
  # plot figure
  plot(Rb, red_ratio_equil*m, type = "l", col = "red", main = "RED and BLUE in neighborhood as a function of Rb",
       xlab = "Tolerance of BLUE", ylab = "Number of people in neighborhood")
  lines(Rb, blue_ratio_equil*m, type = "l", col = "blue")
}

# run body
run_bounded_neigh()