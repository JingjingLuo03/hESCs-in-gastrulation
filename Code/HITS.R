data<-linklist
nodes <- unique(c(data$regulatoryGene, data$targetGene))
num_nodes <- length(nodes)
adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes)
for (i in 1:nrow(data)) {
  node_index <- match(data[i, "regulatoryGene"], nodes)
  target_index <- match(data[i, "targetGene"], nodes)
  adj_matrix[node_index, target_index] <- 1
}
authority <- rep(1, num_nodes)
hub <- rep(1, num_nodes)
max_iter <- 20
tolerance <- 1e-2
for (iter in 1:max_iter) {
  prev_authority <- authority
  prev_hub <- hub
  for (i in 1:num_nodes) {
    authority[i] <- sum(hub[adj_matrix[, i] == 1])
  }
  
  for (i in 1:num_nodes) {
    hub[i] <- sum(authority[adj_matrix[i, ] == 1])
  }
  
  authority <- authority / sqrt(sum(authority^2))
  hub <- hub / sqrt(sum(hub^2))
  
  if (max(abs(authority - prev_authority)) < tolerance && max(abs(hub - prev_hub)) < tolerance) {
    break
  }
}
node_data <- data.frame(Node = nodes, Authority = authority, Hub = hub)






