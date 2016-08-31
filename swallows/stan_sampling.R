# 1. Turn processed data into a form compatible with Stan
# 2. Sample from Stan





#Load processed data and relevant libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("processed_data/pre_swallow.RData")





#Reformat data for Stan
#######################
start = strptime("2014-07-19", format = "%Y-%m-%d")

N = sapply(times,length)
ints_mat = matrix(data = 0, nrow = choose(length(tag_ids),2), ncol = max(N))
dur_mat = matrix(data = 0, nrow = choose(length(tag_ids),2), ncol = max(N))
active_mat = matrix(data = 0, nrow = choose(length(tag_ids),2), ncol = 2)
for(i in 1:choose(length(tag_ids),2)) {
   active_mat[i,] = as.numeric(difftime(active_dyads[[i]],start,unit="hours"))
   if(N[i] != 0) {
      ints_mat[i,1:N[i]] = as.numeric(difftime(times[[i]],start,unit="hours"))
      dur_mat[i,1:N[i]] = duration[[i]]
   }
}

#Gender pairings: FF = 0, MF = 1, MM = 2
gender_stan = (covariates[match(indices[,1],covariates$Tag),5] == "M") + 
   (covariates[match(indices[,2],covariates$Tag),5] == "M")










#Stan Sampler
#############

pre_load = stanc(file = "swallow.stan")
barn_model = stan_model(stanc_ret = pre_load)
expose_stan_functions(barn_model)



#Setup for Stan
data_list <- list(D = choose(length(tag_ids),2), N = N, durations = dur_mat, 
                  interactions = ints_mat, active_times = active_mat, 
                  gender = gender_stan)
param_init <- list(c = rep(0.5,3), k = rep(1,8), sparsity = 0.1, t_prob = 0.001)
stan_init <- list(param_init)



#MCMC Sampling
start_t <- Sys.time()
fit <- sampling(barn_model, data = data_list, init = stan_init, iter = 10000, 
                  warmup = 1000, chains = 1)
stop_t <- Sys.time()
dur_t <- stop_t - start_t



#Save results
#############
#print(fit, pars = c("c","k","sparsity","t_prob"))

save(fit, data_list, param_init, stan_init,
     file = "processed_data/stan_swallow_full.RData")

temp <- extract(fit)
est_params = list(k = temp$k, c = temp$c, sparsity = temp$sparsity, 
                  t_prob = temp$t_prob)
sample_ships = t(sapply(1:136, function(i) sapply(1:45, function(j) 
   quantile(temp$ships[,i,j],0.5))))
q05_ships = t(sapply(1:136, function(i) sapply(1:45, function(j) 
   quantile(temp$ships[,i,j],0.05))))
q95_ships = t(sapply(1:136, function(i) sapply(1:45, function(j) 
   quantile(temp$ships[,i,j],0.95))))

int_fracs <- matrix(data = 0, nrow = choose(length(tag_ids),2), ncol = max(N))
for(d in 1:length(N)) {
   if(N[d] != 0) {
      if(N[d] != 1) {
         int_fracs[d,1:N[d]] = diff(calc_breaks(N[d], ints_mat[d,1:(N[d])], 
                                       dur_mat[d,1:N[d]], active_mat[d,]))
      } else {
         int_fracs[d,1] = diff(calc_breaks(N[d], ints_mat[d,1:(N[d])], 
                                           dur_mat[d,1:N[d]], active_mat[d,]))
      }
   }
}
normalized <- function(f_probs, weights) {
   if(sum(weights) != 0) {
      sum(f_probs*weights)/sum(weights)
   } else {
      0
   }
}
#Probability of friendship integrated over entire time interval
dyad_friends <- sapply(1:136, function(i) 
   normalized(sample_ships[i,1:N[i]],int_fracs[i,1:N[i]]))

#Midpoints of intervals
ship_mids <- list()
for(d in 1:length(N)) {
   if(N[d] != 0) {
      ship_mids[[d]] = ship_times(N[d], ints_mat[d,1:N[d]], 
                                  dur_mat[d,1:N[d]], active_mat[d,])
   } else {
      ship_mids[[d]] = 0
   }
}

save(data_list, sample_ships, q05_ships, q95_ships, est_params, ship_mids,
     dyad_friends, file = "processed_data/stan_swallow_simp.RData")
