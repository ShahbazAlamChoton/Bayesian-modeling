#'[Red]
#'@Purple    
#'*Green*

###### data setup
get_row_wise_adjacency_location <- function(row) {
  if (any(row == 1)) {
    return(which(row == 1)) 
  }
}

unitization <- function(x) {
  
  min.x = min(x, na.rm = T)
  max.x = max(x, na.rm = T)
  
  return((x - min.x)/(max.x - min.x))
}

geweke_stat <- function(x){
  geweke.z = unname(coda::geweke.diag(x)$z)
  p.val = 2*(1-pnorm(abs(geweke.z)))
  return(c(geweke.z = geweke.z, p.val = p.val))
}

prob_matrix <- function(grid_length, m) {
  
  prob_mat <- matrix(0, nrow = grid_length, ncol = grid_length)
  
  for(i in 1 : m) {
    prob_mat[i, setdiff(1:(m+i), i)] = round((m + (i-1))^-1, 3)
  }
  
  for(i in (m+1) : (grid_length - m)) {
    prob_mat[i, setdiff((i-m):(i+m), i)] = round((2 * m)^-1, 3)
  }
  
  
  for( i in (grid_length - m + 1) : grid_length) {
    prob_mat[i, setdiff((i-m):grid_length, i)] = round(length(setdiff((i-m):grid_length, i))^-1, 3)
  }
  
  return(prob_mat)
}
