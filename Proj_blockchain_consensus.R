
### The simulation of consensus process of Blockchain
### Copyright (C): Te Tan, University of Zurich, March 2018

library(igraph)

# System parameters
n <- 200 # network size
net_type <- 0 # network type: 0 is Erdos-Renyi, 1 is Scale-free
k <- c(2, 4, 8, 16, 32) # average degree of Erdos-Renyi network
pa <- c(1, 2, 5, 8, 10) # power of preferential attachment for Scale-free network

iter <- 50 # number of iterations, as well as number of blocks
t_b <- 1000 # block generation time units
lamda_t_power <- sort(seq(-1, 5, by = 0.2), decreasing = FALSE)
power_10 <- function(a) { return(10^a) }
lamda_t <- lapply(lamda_t_power, power_10) # Block transmission time parameter
lamda_c <- sort(seq(0.5, 3, by = 0.5), decreasing = FALSE) # power distribution parameter


# initialize the network environment
init_env <- function(n, net_type, para, lamda_c, lamda_t) {
  if (net_type == 0) {
    net <- erdos.renyi.game(n, (n*para)/2, type = "gnm")
  } else {
    net <- sample_pa(n, power = para)
  }
  # initialize nodes attributes
  c_power <- rexp(n, rate = lamda_c) # sample from Exp to get computing power
  c_accum <- cumsum(c_power) # compute accumulation to support selecting which node generate block
  V(net)$c_power <- c_power
  V(net)$c_accum <- c_accum
  net$max_accum <- c_accum[n]
  V(net)$remain_t <- 0 # remaining time for each node to propagate current chain
  V(net)$height <- 0 # the height of blocktree
  
  E(net)$time <- rexp(length(E(net)), rate = lamda_t^-1)
  
  net$chain <- matrix(0, n, iter)
  net$branches <- 1
  net$off_block <- 0
  
  return(net)
}

# select one node to min a new block accordign to its computing power
generate_new_block <- function(net) {
  rand <- runif(1, min = 0, max = net$max_accum) # generate a random seed
  # select the node according to computing power accumulation, add the new block to its chain
  index <- V(net)[V(net)$c_accum > rand][1]
  
  return(index)
}

# propagate the blocktree from one node to its neighbor
propagate <- function(net, from, to) {
  branch <- FALSE
  if (length(net$chain[to,][net$chain[to,]!=0]) == 0) return(list(branch=branch,off_block=0))
  fork_position <- which(net$chain[to,]!=0)[1]
  
  for (i in 1:min(length(net$chain[from,][net$chain[from,]!=0]), length(net$chain[to,][net$chain[to,]!=0]))) {
    if (net$chain[from, i] == net$chain[to, i]) next
    branch <- TRUE
    fork_position <- i-1
    break
  }
  off_block <- length(net$chain[to,][net$chain[to,]!=0]) - fork_position
  
  return(list(branch=branch,off_block=off_block))
}

# run blockchain on one network
run_on_net <- function(net) {
  # every iteration will generate a new block in the network
  for (i in 1:iter) {
    # generate a new block with id of i
    index <- generate_new_block(net)
    new_block_position <- which(net$chain[index,]==0)[1]
    net$chain[index, new_block_position] <- i
    V(net)[index]$height <- new_block_position
    V(net)$remain_t <- V(net)$remain_t + t_b # add new remaining transmission time for each node
    
    while (TRUE) {
      # iterate over each node to do the block propagation process
      propagate_happened <- FALSE
      for (j in 1:n) {
        height <- V(net)[j]$height
        if (height == 0) next
        Neighbors <- neighbors(net, V(net)[j])
        if (length(Neighbors) == 0) next
        Edges <- E(net)[from(V(net)[j])]
        
        max_trans_time <- 0
        for (k in 1:length(Neighbors)) {
          nei_index <- Neighbors[k]
          # if local blocktree is higher and remaining time is enough, propagate it to neighbor
          if (height > V(net)[nei_index]$height & V(net)[j]$remain_t >= E(net)[V(net)[j] %--% V(net)[nei_index]]$time) {
            propagate_happened <- TRUE
            
            result <- propagate(net, V(net)[j], V(net)[nei_index]) # propagate the block
            if (result$branch) {
              net$branches <- net$branches + 1
              net$off_block <- net$off_block + result$off_block
            }
            net$chain[nei_index,] <- net$chain[j,] # update the neighbor's local blockchain
            V(net)[nei_index]$remain_t <- V(net)[j]$remain_t - E(net)[V(net)[j] %--% V(net)[nei_index]]$time
            V(net)[nei_index]$height <- height
            
            if (E(net)[V(net)[j] %--% V(net)[nei_index]]$time > max_trans_time) {
              max_trans_time <- E(net)[V(net)[j] %--% V(net)[nei_index]]$time
            }
          }
        }
        if (propagate_happened == TRUE) {
          V(net)[j]$remain_t <- V(net)[j]$remain_t - max_trans_time
        }
      }
      # if no more propagation happened, go to next iteration
      if (propagate_happened == FALSE) break
    }
  }
  
  return(net)
}

# run experiemtns with parameters
run_experiments <- function() {
  net <- init_env(n, 0, 8, 1, 10^2)
  net <- run_on_net(net)
  print(paste("max blockchain height: ", max(V(net)$height)))
  print(paste("branches during consensus: ", net$branches))
  print(paste("blocks discarded: ", net$off_block))
}

# start
run_experiments()

x <- c(10^2, 10^2.2, 10^2.4, 10^2.5, 10^2.7, 10^2.8, 10^3, 10^3.1, 10^3.2, 10^3.3, 10^3.4, 10^3.5, 10^3.6, 10^3.7, 10^3.8)
y1 <- matrix(c(1,2,38,209,536,536,468,440,418,293,250,179,164,116,120,1,1,3,6,230,281,261,312,482,544,744,637,576,372,332,1,1,1,1,2,3,57,103,318,403,722,1140,715,664,443),15,3)
matplot(x, y1, main = expression(paste('Nr of branches as a function of ', lambda[t])),type = "b", pch = 15:17, col = 1:3, log = "x", ylim = c(1,1200), xlab = expression(lambda[t]), ylab = 'Nr of branches')
legend("topleft", legend = c('k = 2', 'k = 4', 'k = 8'),inset = 0.05, pch = 15:17, col = 1:3, bg = ("white"), horiz = F)

y2 <- matrix(c(0,1,37,210,711,560,488,452,479,378,311,242,236,142,160,0,2,5,229,280,268,322,497,590,1045,832,785,699,405,0,0,0,1,2,56,112,342,436,832,1198,894,859,460),15,3)
matplot(x, y2, main = expression(paste('Nr of discarded blocks as a function of ', lambda[t])),type = "b", pch = 12:14, col = 1:3, log = "x", ylim = c(1,1250), xlab = expression(lambda[t]), ylab = 'Nr of discarded blocks')
legend("topleft", legend = c('k = 2', 'k = 4', 'k = 8'),inset = 0.05, pch = 12:14, col = 1:3, bg = ("white"), horiz = F)

y3 <- matrix(c(0.84,0.8,0.64,0.56,0.48,0.42,0.34,0.32,0.3,0.24,0.2,0.18,0.16,0.14,0.1,0.98,0.98,0.96,0.9,0.72,0.7,0.56,0.46,0.44,0.4,0.32,0.28,0.26,0.24,0.2,1.0,1.0,1.0,1.0,0.96,0.9,0.78,0.74,0.7,0.62,0.56,0.48,0.38,0.34,0.3),15,3)
matplot(x, y3, main = expression(paste('Finalized block ratio as a function of ', lambda[t])),type = "b", pch = 20:22, col = 1:3, log = "x", ylim = c(0,1.2), xlab = expression(lambda[t]), ylab = 'Finalized block ratio')
legend("topright", legend = c('k = 2', 'k = 4', 'k = 8'),inset = 0.05, pch = 20:22, col = 1:3, bg = ("white"), horiz = F)

