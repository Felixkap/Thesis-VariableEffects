how_many <- function(scenarios){
  
  n_vars <- scenarios[['n_vars']]
  th <- scenarios[['thrs']]
  dat <- rmvnorm(n, mean = rep(0, n_vars), sigma = CSgenerate(n_vars, 0))
  logical_list <- vector(mode = 'list', length = n_vars)
  for (j in 1:n_vars) {
    logical_vec <- (-th < dat[,j]) & (dat[,j] < th)
    logical_list[[j]] <- logical_vec
  }
  
  that_many <- c(that_many = sum(Reduce("&", logical_list)))
  
  return(that_many)
}

reps <- 1e3
n_vars <- seq(2, 12, 4)
thrs <- c(0.4, 1)
scenarios <- data.frame(expand.grid(n_vars, thrs))
colnames(scenarios) = c("n_vars", "thrs")
scenarios <- split(scenarios, seq(nrow(scenarios)))
result <- lapply(X = scenarios, FUN = how_many)


Map(c, scenarios, result)