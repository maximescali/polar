# Load and setup rstan
library(rstan)
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# Paths to the stan files
eigenmodel_file <- '/Users/maximescali/Desktop/Research/MScthesis/sampling/polar/network_eigenmodel.stan'

# Expose the user defined functions so that we have access to lower_tri in R.
expose_stan_functions(eigenmodel_file, show_compiler_warnings = FALSE)

# Load the network adjacency matrix   
Y <- as.matrix(read.table('/Users/maximescali/Desktop/Research/MScthesis/sampling/polar/hoff.dat', sep=' '))

# Create the variables that will be passed into Stan.
p <- 270 # Number of nodes in the network
y <- lower_tri(Y) # The lower triangular entries of the adjacency matrix Y
k <- 3 # Latent dimension
sd_lam <- sqrt(p) # Prior standard deviation for lam parameter
sd_c <- 10 # Prior standard deviation for c parameter

# The following list stores all the data we need to pass to the Stan sampler. 
dat <- list('p'=p, 'y'=y, 'k'=k, 'sd_lam'=sd_lam, 'sd_c'=sd_c)

chains <- 1 # Number of Markov chains run
warmup <- 2500 # Number of warmup iterations per chain 
num_iter <- 10000 # Number of post warmup iterations per chain 

# Do posterior sampling for the standard network eigenmodel
eigenmodel_fit <- stan(file = eigenmodel_file, data = dat, 
                      iter = num_iter + warmup, chains = chains, refresh=1, 
                      thin=1, warmup=warmup)

eigenmodel_draws <- extract(eigenmodel_fit)



















# # OLD CODE: 
# 
# 
# # Produce posterior samples and then posterior means of the connection
# # probabilities for the held out dyads. 
# 
# ind_exclude <- setdiff(1:length(y), ind_include) 
# 
# eigenmodel_prob_ex_draws <- 
#   matrix(data=rep(0, length(y)*num_iter), nrow=num_iter)
# 
# sparse_eigenmodel_prob_ex_draws <- 
#   matrix(data=rep(0, length(y)*num_iter), nrow=num_iter)
# 
# # Degree distribution of posterior predictive samples
# eigenmodel_pp_degree <- matrix(data=rep(0, 80*num_iter), nrow=num_iter)
# sparse_eigenmodel_pp_degree <- matrix(data=rep(0, 80*num_iter), nrow=num_iter)
# 
# 
# for(i in 1:num_iter){
#   eigenmodel_prob_ex_draws[i,] <- 
#     lower_tri(pnorm(eigenmodel_draws$c[i] + eigenmodel_draws$Q[i,,] %*% diag(eigenmodel_draws$lam[i,]) %*% t(eigenmodel_draws$Q[i,,])))
#   
#   Ytmp <- diag(p)
#   Ytmp[lower.tri(Ytmp, diag=FALSE)] <- rbinom(n=p*(p-1)/2, p=eigenmodel_prob_ex_draws[i,], size=1)
#   Ytmp <- Ytmp + t(Ytmp) - diag(p)
#   ind_nonzero <- as.integer(names(table(rowSums(Ytmp) - 1))) + 1
#   eigenmodel_pp_degree[i, ind_nonzero] <- table(rowSums(Ytmp) - 1)
#   
#   sparse_eigenmodel_prob_ex_draws[i,] <- 
#     lower_tri(pnorm(sparse_eigenmodel_draws$c[i] + sparse_eigenmodel_draws$Q[i,,] %*% diag(sparse_eigenmodel_draws$lam[i,]) %*% t(sparse_eigenmodel_draws$Q[i,,])))
#   
#   Ytmp <- diag(p)
#   Ytmp[lower.tri(Ytmp, diag=FALSE)] <- rbinom(n=p*(p-1)/2, p=sparse_eigenmodel_prob_ex_draws[i,], size=1)
#   Ytmp <- Ytmp + t(Ytmp) - diag(p)
#   ind_nonzero <- as.integer(names(table(rowSums(Ytmp) - 1))) + 1
#   sparse_eigenmodel_pp_degree[i, ind_nonzero] <- table(rowSums(Ytmp) - 1)
#   
#   
#   print(i)
# }
# 
# eigenmodel_prob_ex_estimates <- 
#   apply(X=eigenmodel_prob_ex_draws, MARGIN=2, FUN=mean)
# 
# sparse_eigenmodel_prob_ex_estimates <- 
#   apply(X=sparse_eigenmodel_prob_ex_draws, MARGIN=2, FUN=mean)
# 
# # plot ideas: 
# # - degree distribution 
# # - histogram comparison of entries of Q 
# # - latent space plots
# 
# quants <- c(.025, .5, .975)
# eigenmodel_pp_degree_quantiles <- apply(X=eigenmodel_pp_degree, MARGIN=2, quantile , probs = quants, na.rm = TRUE)
# sparse_eigenmodel_pp_degree_quantiles <- apply(X=sparse_eigenmodel_pp_degree, MARGIN=2, quantile , probs = quants, na.rm = TRUE)
# 
# 
# # Can plot the observed degree distribution as follows: 
# plot(table(rowSums(Y) - 1), lwd=5)
# 
# points(x=0:79, y=eigenmodel_pp_degree_quantiles[1,], col="red")
# points(x=0:79, y=eigenmodel_pp_degree_quantiles[2,], col="red")
# points(x=0:79, y=eigenmodel_pp_degree_quantiles[3,], col="red")
# 
# 
# points(x=0:79, y=sparse_eigenmodel_pp_degree_quantiles[1,], col="blue")
# points(x=0:79, y=sparse_eigenmodel_pp_degree_quantiles[2,], col="blue")
# points(x=0:79, y=sparse_eigenmodel_pp_degree_quantiles[3,], col="blue")
# 
# 
