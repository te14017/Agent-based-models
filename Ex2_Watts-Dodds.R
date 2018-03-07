
### Watts-Dodds model simulation for riot problem
### Copyright (C): Te Tan, University of Zurich, March 2018

library(igraph)

# model parameters
p_list <- sort(seq(0.004, 0.2, by = 0.01), decreasing = FALSE) # probabilities of degree in random graph
n <- 5000 # size of population
n0 <- 200 # initial number of people joining riot
tau <- 0.12 # the homogeneous threshold of all individuals

# initialize Erdos-Renyi network
init_net <- function(n, n0, p) {
  net <- erdos.renyi.game(n, p, type = "gnp")
  V(net)$value <- 0
  # initially joint people's indexes
  riot_index <- sample(seq(1, n), n0)
  V(net)[riot_index]$value <- 1
  
  return(net)
}

# run function trying different probabilities
run_watts <- function() {
  nr_of_riot_final <- vector()
  
  for (i in 1:length(p_list)) {
    net <- init_net(n, n0, p_list[i])
    
    riot_ratio <- n0/n
    previous_ratio <- 0
    while (TRUE) {
      previous_ratio <- riot_ratio
      if (length(V(net)[V(net)$value==0]) == 0) {
        break
      }
      # iterate over neighbors of each vertice who do not yet join riot
      for (j in 1:n) {
        if (V(net)[j]$value != 0) next
        Neighbors <- neighbors(net, V(net)[j])
        if (length(Neighbors) == 0) next
        
        # if the ratio of neighbors who joined riots are bigger than threshold, then he joins
        if (sum(Neighbors$value)/length(Neighbors) >= tau) {
          V(net)[j]$value <- 1
        }
      }
      # if previous ratio equals current one, evolution terminates
      riot_ratio <- sum(V(net)$value)/n
      if (previous_ratio == riot_ratio) {
        break
      }
    }
    nr_of_riot_final[i] <- sum(V(net)$value)/n
  }
  plot(p_list, nr_of_riot_final, type = "l", main = "Ratio of people joint riot over P value of network",
       xlab = "P value of Erdos-Renyi network", ylab = "Ratio of people joint riot")
}

# run body
run_watts()