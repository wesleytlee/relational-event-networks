# 1. Sample from Stan





#Load processed data and relevant libraries
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("processed_data/stan_pre.RData")



#Sampling
#########
pre_load = stanc(file = "students.stan")
mit_model = stan_model(stanc_ret = pre_load)
start.t = Sys.time()
fit <- sampling(mit_model, data = data_list, init = list(param_init), 
                iter = 10000, warmup = 1000, chains = 1, verbose = TRUE)
stop.t <- Sys.time()
stop.t - start.t
output <- extract(fit)
expose_stan_functions(mit_model)



#Estimate network using mean of estimated parameter values
est_ships <- list()
for(d in 1:length(data_list$N)) {
   if(data_list$N[d] > 1) {
      est_ships[[d]] = infer_ships(data_list$N[d], 
         data_list$interactions[d,1:data_list$N[d]], 
         data_list$durations[d,1:data_list$N[d]], 
         data_list$active_times[d,], data_list$same_floor[d], 
         data_list$same_year[d], data_list$sun_adjust, data_list$term_dates, 
         data_list$k_2, colMeans(output$c_0), mean(output$c_1), 
         colMeans(output$c_2), mean(output$c_3), mean(output$k_0), 
         colMeans(output$k_1), colMeans(output$k_3), mean(output$t_prob), 
         mean(output$s_0), colMeans(output$s))
   } else {
      est_ships[[d]] = 0
   }
}


#Midpoint of intervals 
ship_ints <- list()
for(d in 1:(data_list$D)) {
   if(data_list$N[d] != 0) {
      ship_ints[[d]] = ship_times(data_list$N[d], 
                        data_list$interactions[d,1:(data_list$N[d])], 
                        data_list$durations[d,1:(data_list$N[d])], 
                        data_list$active_times[d,])
   }else {
      ship_ints[[d]] = 0
   }
}

#Length of each interval
int_lengths <- matrix(data = 0, nrow = data_list$D, ncol = max(data_list$N))
for(d in 1:(data_list$D)) {
   if(data_list$N[d] != 0) {
      if(data_list$N[d] != 1) {
         int_lengths[d,1:(data_list$N[d])] = 
            diff(calc_breaks(data_list$N[d], 
                             data_list$interactions[d,1:(data_list$N[d])], 
                             data_list$durations[d,1:(data_list$N[d])], 
                             data_list$active_times[d,]))
      } else {
         int_lengths[d,1] = 
            diff(calc_breaks(data_list$N[d], 
                             data_list$interactions[d,1:(data_list$N[d])], 
                             data_list$durations[d,1:(data_list$N[d])], 
                             data_list$active_times[d,]))
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

#Probability of a connection integrated over active time interval
dyad_friends <- sapply(1:(data_list$D), function(i) 
   normalized(est_ships[[i]],int_lengths[i,1:(data_list$N[i])]))

save(est_ships,fit,dyad_friends,ship_ints, 
     file="processed_data/stan_output.RData")
