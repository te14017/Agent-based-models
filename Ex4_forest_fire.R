
### Forest fire model simulation
### Copyright (C): Te Tan, University of Zurich, March 2018

# model parameters
n <- 200 # size of the forest
p <- sort(seq(0.1, 1.0, by = 0.1), decreasing = FALSE) # density of trees
burn_times <- 10 # how many burns we consider for each p

# initialize the forest
init_forest <- function(n, p) {
  forest <- matrix(0, nrow = n, ncol = n)
  
  # place trees randomly inside the forest
  tree_index <- sample(seq(1, n*n), p*n*n)
  for (i in 1:length(tree_index)) {
    forest[tree_index[i]] <- 1
  }
  return(forest)
}

# if index is legal and is a tree
is_tree <- function(i, forest) {
  if (i < 1 | i > n*n) { return(FALSE) }
  return(forest[i] == 1)
}

# get indexes of neighbor trees
get_neighbors_index <- function(i, forest) {
  indexes <- vector()
  j <- 0
  up_left <- i-n-1
  if (is_tree(up_left, forest)) {
    j <- j + 1
    indexes[j] <- up_left
  }
  up <- i-n
  if (is_tree(up, forest)) {
    j <- j + 1
    indexes[j] <- up
  }
  up_right <- i-n+1
  if (is_tree(up_right, forest)) {
    j <- j + 1
    indexes[j] <- up_right
  }
  left <- i-1
  if (is_tree(left, forest)) {
    j <- j + 1
    indexes[j] <- left
  }
  right <- i+1
  if (is_tree(right, forest)) {
    j <- j + 1
    indexes[j] <- right
  }
  down_left <- i+n-1
  if (is_tree(down_left, forest)) {
    j <- j + 1
    indexes[j] <- down_left
  }
  down <- i+n
  if (is_tree(down, forest)) {
    j <- j + 1
    indexes[j] <- down
  }
  down_right <- i+n+1
  if (is_tree(down_right, forest)) {
    j <- j + 1
    indexes[j] <- down_right
  }
  return(indexes)
}

# get burn target tree index
get_target_index <- function(forest) {
  tree_index <- which(forest == 1)
  if (length(tree_index) == 0) { return(0) }
  random_pick <- sample(seq(1, length(tree_index)), 1)
  return(tree_index[random_pick])
}

# light one tree to burn its area, return the forest after burnt
one_burn <- function(forest, i) {
  forest[i] <- 0
  neighbors <- get_neighbors_index(i, forest)
  if (length(neighbors) == 0) { return(forest) }
  
  # burn its neighbors
  forest[neighbors] <- 0
  
  previous_sum <- 0
  sum_forest <- sum(forest)
  while (previous_sum != sum_forest) {
    previous_sum <- sum_forest
    if (length(neighbors) == 0) { break }
    new_expand_neighbors <- get_neighbors_index(neighbors[1], forest)
    if (length(new_expand_neighbors) > 0) { forest[new_expand_neighbors] <- 0 }
    
    if (length(neighbors) > 1) {
      for (j in 2:length(neighbors)) {
        expand_neighbors <- get_neighbors_index(neighbors[j], forest)
        if (length(expand_neighbors) > 0) {
          forest[expand_neighbors] <- 0
          new_expand_neighbors[(length(new_expand_neighbors)+1):(length(new_expand_neighbors)+length(expand_neighbors))] <- expand_neighbors
        }
      }
    }
    
    neighbors <- new_expand_neighbors
    sum_forest <- sum(forest)
  }
  
  return(forest)
}

# run function
run_forest_fire <- function() {
  average_num_burnt <- vector()
  
  for(i in 1:length(p)) {
    forest <- init_forest(n, p[i])
    # record total burnt number of multiple burns to compute the average
    burnt_num_total <- 0
    for(j in 1:burn_times) {
      forest_rep <- rep(forest)
      light_target <- get_target_index(forest_rep)
      forest_after_burn <- one_burn(forest_rep, light_target)
      burnt_num_total <- burnt_num_total + (sum(forest) - sum(forest_after_burn))
    }
    average_num_burnt[i] <- burnt_num_total/burn_times
  }
  
  plot(p, average_num_burnt, type = "l", main = "average number of trees burnt over P",
       xlab = "The density of trees: P", ylab = "average number of burnt trees")
}

# run body
run_forest_fire()

