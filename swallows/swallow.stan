functions {

   //Returns a mod b
   real mod(real a, int b) {
      return a - floor(a/b)*b;
   }

   //Returns min(a,b)
   real minR(real a, real b) {
      return (b<a)*b + (b>=a)*a;
   }

   //Returns max(a,b)
   real maxR(real a, real b) {
      return (b<a)*a + (b>=a)*b;
   }

   //Return period/session associated with a given time
   int period(real t) {
      if(t < 12) {
         //Garbage case
         return 0;
      } else if(t < 24) {
         return 1;
      } else if (t < 36) {
         return 2;
      } else if (t < 48) {
         return 3;
      } else if (t < 60) {
         return 4;
      } else if (t < 72) {
         return 5;
      } else if (t < 84) {
         return 6;
      } else if (t < 96) {
         return 7;
      } else {
         return 8;
      }
   }

   //Amount of time in each period elapsed between t1 and t2
   row_vector p_time(real t1, real t2) {
      row_vector[8] times;
      //17-20 is Day 1 PM
      times[1] = maxR(minR(t2,20) - maxR(t1,17),0);
      //30-33 is Day 2 AM
      times[2] = maxR(minR(t2,33) - maxR(t1,30),0);
      //41-44 is Day 2 PM
      times[3] = maxR(minR(t2,44) - maxR(t1,41),0);
      //54-57 is Day 3 AM
      times[4] = maxR(minR(t2,57) - maxR(t1,54),0);
      //65-68,89-92 is Day 3 PM
      times[5] = maxR(minR(t2,68) - maxR(t1,65),0);
      //78-81 is Day 4 AM
      times[6] = maxR(minR(t2,81) - maxR(t1,78),0);
      //89-92 is Day 4 PM
      times[7] = maxR(minR(t2,92) - maxR(t1,89),0);
      //102-105 is Day 5 AM
      times[8] = maxR(minR(t2,105) - maxR(t1,102),0);
      return times;
   }

   //Calculate breaks that partitions full time interval into smaller intervals
   real[] calc_breaks(int N, real[] interactions, real[] durations,
                        real[] active_times) {
      real brk_points[max(N+1,2)];
      brk_points[1] = active_times[1]-2.0/6; //Buffer
      if(N > 1) {
         for(n in 2:N) {
            brk_points[n] = (interactions[n-1] + durations[n-1]
                              + interactions[n])/2;
         }
      }
      brk_points[N+1] = active_times[2]+-2.0/6;
      return brk_points;
   }

  //Calculate midpoints of partitioned intervals
   real[] ship_times(int N, real[] interactions, real[] durations,
                     real[] active_times) {
      real brk_points[N+1];
      real times[N];
      if(N != 0) {
         brk_points = calc_breaks(N, interactions, durations, active_times);
         for(n in 1:N) {
         times[n] = (brk_points[n] + brk_points[n+1])/2;
         }
      } 
      return times;
   }
}

data {
   int D; //Number of dyads
   int<lower=0> N[D]; //Number of interactions for each dyad
   real<lower=0> interactions[D,max(N)]; //Interaction times (HOURS)
   real<lower=0> durations[D,max(N)]; //Duration of interactions (HOURS)
   real<lower=0> active_times[D,2]; //Period for which dyad is active
   int<lower=0,upper=2> gender[D]; //Dyadic gender pairing (0: FF, 1: MF, 2: MM)
}

transformed data {
   int<lower=0,upper=8> periods[D,max(N)]; //Period of each interaction
   real<lower=-0.01> brk_points[D,max(N)+1]; //Break times
   real<lower=0> diffHours[D,max(N)]; //Difference between interaction times
   //Difference between break times in each period
   row_vector<lower=0>[8] period_time[D,max(N)]; 

   //Initialize arrays
   periods = rep_array(0,D,max(N));
   brk_points = rep_array(0.0,D,max(N)+1);
   diffHours = rep_array(0.0,D,max(N));
   period_time = rep_array(rep_row_vector(0.0,8),D,max(N));

   //For each dyad
   for (d in 1:D) {
      brk_points[d,1] = active_times[d,1]-2.0/6; //Buffer
    
      if(N[d] != 0) {       
         //Calculate breaks
         for(n in 2:(N[d]+1)) {
            if (n <= N[d]) {
               brk_points[d,n] = (interactions[d,n-1] + durations[d,n-1]
                                 + interactions[d,n])/2;
            } else {
               brk_points[d,n] = active_times[d,2]+2.0/6;
            }
         }
      
         //Calculate characteristics about interactions and intervals
         for(n in 1:N[d]) {
            periods[d,n] = period(interactions[d,n]);
            if(n == 1) {
               diffHours[d,n] = (brk_points[d,2] - brk_points[d,1])/2;
            } else {
               diffHours[d,n] = (brk_points[d,n+1] - brk_points[d,n-1])/2;
            }
            //Ignore duration while interaction occurs
            period_time[d,n] = p_time(brk_points[d,n],interactions[d,n]) 
               + p_time(interactions[d,n]+durations[d,n],brk_points[d,n+1]);
         }
       
      //If no interactinos
      } else {
         brk_points[d,2] = active_times[d,2]+2.0/6;
         period_time[d,1] = p_time(brk_points[d,1],brk_points[d,2]);
      }
    
   }
}

parameters {
   vector<lower=0>[8] k; //Baseline interaction rates
   real<lower=0> c[3]; //Relationship multiplier
   real<lower=0,upper=1> sparsity; //CTMC sparsity
   real<lower=0> t_prob; //CTMC transition probability ("q")
}

model {
   real log_obs[D]; //Vector of log-likelihoods
   row_vector[2] probs_0; //Probabilities to transition to 0
   row_vector[2] probs_1; //Probabilities to transition to 1

   real integral; //Placeholder variable for calculating various integrals

   log_obs = rep_array(0.0,D);

   //Priors
   c[1] ~ exponential(2);
   c[2] ~ exponential(2);
   c[3] ~ exponential(2);
   for(i in 1:8) {
      k[i] ~ exponential(1);
   }
   sparsity ~ beta(1,9);
   t_prob ~ exponential(1000);
  
   for(d in 1:D) {
      if(N[d] != 0) {
         //Forward variables from forward-backward algorithm
         row_vector[2] forward[N[d]];
        
         //Transition probabilities between intervals
         probs_0[1] = log_sum_exp(log1m(sparsity), 
                        log(sparsity)-t_prob*diffHours[d,1]);
         probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[d,1]);
         probs_1[1] = log1m_exp(probs_0[1]);      
         probs_1[2] = log1m_exp(probs_0[2]);

         //Integral of rate over period
         integral = period_time[d,1]*k;

         //Calculate forward variables
         forward[1][1] = log_sum_exp(log1m(sparsity)+probs_0[1], 
                                       log(sparsity)+probs_0[2])
                           + log(k[periods[d,1]]) - integral;
         forward[1][2] = log_sum_exp(log1m(sparsity) + probs_1[1],
                                       log(sparsity) + probs_1[2])
                           + log((1+c[gender[d]+1])*k[periods[d,1]])
                           - (1+c[gender[d]+1])*integral;
         
         if(N[d] > 1) {
            for (n in 2:N[d]) {
               probs_0[1] = log_sum_exp(log1m(sparsity), 
                                          log(sparsity)-t_prob*diffHours[d,n]);
               probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[d,n]);
               probs_1[1] = log1m_exp(probs_0[1]);      
               probs_1[2] = log1m_exp(probs_0[2]);

               integral = period_time[d,n]*k;

               forward[n][1] = log_sum_exp(forward[n-1] + probs_0) 
                                 + log(k[periods[d,n]]) - integral;
               forward[n][2] = log_sum_exp(forward[n-1] + probs_1)
                                 + log((1+c[gender[d]+1])*k[periods[d,n]])
                                 - (1+c[gender[d]+1])*integral;
            }
         }
         
         //Probability of observed emissions (interactions) sequence
         log_obs[d] = log_sum_exp(forward[N[d]]);
         
      } else {
         //Assumes if no observations then they are not friends
         log_obs[d] = -period_time[d,1]*k;
      }
   }
   
   //Posterior
   target += sum(log_obs);
}

generated quantities {
   //Probability of friendship at selected times (midpoint of intervals)
   real<lower=-1,upper=1> ships[D,max(N)];
  
   ships = rep_array(-1.0,D,max(N));
  
   for(d in 1:D) {
      if(N[d] != 0) {
         //Forward and backward variables
         row_vector[2] forward[N[d]];
         row_vector[2] backward[N[d]];
         real integral;
         
         row_vector[2] probs_0;
         row_vector[2] probs_1;
         row_vector[2] intensities;
         row_vector[2] integrals;
         int m;

         probs_0[1] = log_sum_exp(log1m(sparsity), 
                                    log(sparsity) - t_prob*diffHours[d,1]);
         probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[d,1]);
         probs_1[1] = log1m_exp(probs_0[1]);      
         probs_1[2] = log1m_exp(probs_0[2]);

         integral = period_time[d,1]*k;

         //Calculate forward variables
         forward[1][1] = log_sum_exp(log1m(sparsity) + probs_0[1],
                                       log(sparsity) + probs_0[2])
                           + log(k[periods[d,1]]) - integral;
         forward[1][2] = log_sum_exp(log1m(sparsity) + probs_1[1],
                                       log(sparsity)+probs_1[2]) 
                           + log((1+c[gender[d]+1])*k[periods[d,1]]) 
                           - (1+c[gender[d]+1])*integral;
       
         if(N[d] > 1) {
           for (n in 2:N[d]) {
             probs_0[1] = log_sum_exp(log1m(sparsity), 
                                       log(sparsity)-t_prob*diffHours[d,n]);
             probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[d,n]);
             probs_1[1] = log1m_exp(probs_0[1]);      
             probs_1[2] = log1m_exp(probs_0[2]);


             integral = period_time[d,n]*k;

             forward[n][1] = log_sum_exp(forward[n-1] + probs_0) 
                              + log(k[periods[d,n]]) - integral;
             forward[n][2] = log_sum_exp(forward[n-1] + probs_1) 
                              + log((1+c[gender[d]+1])*k[periods[d,n]]) 
                              - (1+c[gender[d]+1])*integral;
           }
         }
         
         //Calculate backward variables
         backward[N[d]][1] = 0;
         backward[N[d]][2] = 0;
         
         if(N[d] > 1) {
            for (n in 2:N[d]) {
               m = N[d] - n + 1;

               //Probabilities to transition from 0 and 1
               probs_0[1] = log_sum_exp(log1m(sparsity), 
                                          log(sparsity) - t_prob*diffHours[d,m]);
               probs_0[2] = log1m_exp(probs_0[1]);
               probs_1[1] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[d,m]);    
               probs_1[2] = log1m_exp(probs_1[1]);   

               intensities[1] = k[periods[d,m+1]];
               intensities[2] = (1+c[gender[d]+1])*k[periods[d,m+1]];
               integrals[1] = period_time[d,m+1]*k;
               integrals[2] = (1+c[gender[d]+1])*period_time[d,m+1]*k;

               backward[m][1] = log_sum_exp(probs_0 + backward[m+1] 
                                             + log(intensities) - integrals);
               backward[m][2] = log_sum_exp(probs_1 + backward[m+1] 
                                             + log(intensities) - integrals);
            }
         }
         
         for(n in 1:N[d]) {
           ships[d,n] = inv_logit(forward[n][2] + backward[n][2] 
                                    - forward[n][1] - backward[n][1]);
         }

      //If no interactions, probability of a relationship is very low (~ zero) 
      } else {
         ships[d,1] = 0;
      }
   }
}