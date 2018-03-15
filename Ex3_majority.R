
### Majority model simulation for consensus problem (Plus or Minus)
### Copyright (C): Te Tan, University of Zurich, March 2018

# model parameters
n <- 5000 # system size
n_plus <- 0.5 # ratio of popultion with + opinion
n_plus_list <- sort(seq(0.1, 1.0, by = 0.1), decreasing = FALSE)
a2 <- 0.25 # probability of groug being size of 2
a3 <- 0.5
a4 <- 0.25
runs <- 100 # number of runs for computing fraction
stop_threshold <- 10000

# initialize the system
init_population <- function(n, n_plus) {
  population <- rep(0, n)
  
  plus_index <- sample(seq(1, n), n_plus*n)
  for(i in seq(1, n_plus*n)) {
    population[plus_index[i]] <- 1
  }
  return(population)
}

# assign groups
assign_groups <- function(n, a2, a3, a4) {
  groups_vect <- rep(0, n)
  
  counter <- 1
  group_index <- 0
  while (counter <= n) {
    random <- runif(1)[1]
    group_index <- group_index + 1
    
    if (random < a2) {
      groups_vect[seq(counter, min(counter+1, n))] <- group_index
      counter <- counter + 2
    } else if (random >= a2 & random < a2 + a3) {
      groups_vect[seq(counter, min(counter+2, n))] <- group_index
      counter <- counter + 3
    } else {
      groups_vect[seq(counter, min(counter+3, n))] <- group_index
      counter <- counter + 4
    }
  }
  
  return(sample(groups_vect))
}

# run plotting time evolution
run_majority_evolution <- function(plot = TRUE, n_plus = n_plus) {
  population <- init_population(n, n_plus)
  ratio_records <- vector() # collecting ratio evolution
  
  counter <- 1
  previous_ratio <- 0
  ratio_track <- n_plus
  evol_counts <- 0
  
  while (TRUE) {
    previous_ratio <- ratio_track
    groups_vect <- assign_groups(n, a2, a3 ,a4)
    
    for (i in 1:length(unique(groups_vect))) {
     member_index <- which(groups_vect == i)
     if (sum(population[member_index])/length(member_index) >= 0.5) {
       population[member_index] <- 1
     } else {
       population[member_index] <- 0
     }
    }
    
    ratio_track <- sum(population)/n
    ratio_records[counter] <- ratio_track
    counter <- counter + 1
    evol_counts <- evol_counts + 1
    if (ratio_track == previous_ratio | evol_counts > stop_threshold) break
  }
  
  if (plot == TRUE) {
  plot(1:length(ratio_records), ratio_records, type = "l", main = "Ratio of Plus evolution",
       xlab = "Number of iteration", ylab = "Ratio of Plus")
  }
  
  return(sum(population)/n)
}

# compute the fraction of Plus over many runs
run_compute_fraction <- function(n_plus = n_plus) {
  num_of_plus <- 0
  for (i in 1:runs) {
    ratio_stable <- run_majority_evolution(plot = FALSE, n_plus)
    if (ratio_stable == 1.0) num_of_plus <- num_of_plus + 1
  }
  
  return(num_of_plus/runs)
}

# run body
#run_majority_evolution()
ratio_plus_stable <- vector()
for (i in 1:length(n_plus_list)) {
  ratio_plus <- run_compute_fraction(n_plus = n_plus_list[i])
  ratio_plus_stable[i] <- ratio_plus
}
plot(n_plus_list, ratio_plus_stable, type = "l", main = "Stable ratio of Plus over initial ratio",
     xlab = "Initial ratio of Plus", ylab = "Stable ratio of Plus")

