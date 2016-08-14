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
   
   //Indicator function for whether or not an interaction occurs during the day
   //(defined to be 8 AM - 5 PM)
   int isDay(real t) {
      return mod(t,24) >= 8 && mod(t,24) <= 17;
   }
   
   //Which term did an interaction occur in (1-3)
   int period(real t,row_vector term_dates) {
      if(t <= term_dates[2]) {
         return 1;
      } else if (t >= term_dates[5]) {
         return 3;
      } else {
         return 2;
      }
   }
  
   //Amount of daytime elapsed between times t1 and t2
   real daytime(real t1, real t2) {
      real time;
      real t3;
      
      time = floor((t2-t1)/24)*9; //Daytime from full days
      
      //Now calculate daytime between t1 and t3, which are within 24 hours
      t3 = t2 - time*24/9;
      
      //If on the same day
      if (floor((t1-8)/24) == floor((t3-8)/24)) {
         if (isDay(t3)) {
            time = time + t3 - t1;
         } else if (isDay(t1)) {
            time = time + 9 - mod(t1-8,24);
         }
      
      //If on seperate days
      } else {
         if (isDay(t3)) {
            time = time + mod(t3-8,24);
         } else {
            time = time + 9;
         }
         if (isDay(t1)) {
            time = time + 9 - mod(t1-8,24);
         }
      }
      return time;
   }
   
   //Calculate time within each term between t1 and t2
   row_vector p_times(real t1, real t2, row_vector term_dates) {
      row_vector[3] times;
      times[1] = maxR(minR(t2,term_dates[2]) - maxR(t1,term_dates[1]),0);
      times[2] = maxR(minR(t2,term_dates[4]) - maxR(t1,term_dates[3]),0);
      times[3] = maxR(minR(t2,term_dates[6]) - maxR(t1,term_dates[5]),0);
      return times;
   }
  
   //Daytime elapsed within each term between t1 and t2
   row_vector day_terms(real t1, real t2, row_vector term_dates) {
      row_vector[3] ptimes;
      row_vector[3] times;
      ptimes = p_times(t1,t2,term_dates);
      times = rep_row_vector(0.0,3);
      if(ptimes[1] > 0)
         times[1] = daytime(maxR(t1,term_dates[1]),minR(t2,term_dates[2]));
      if(ptimes[2] > 0)
         times[2] = daytime(maxR(t1,term_dates[3]),minR(t2,term_dates[4]));
      if(ptimes[3] > 0)
         times[3] = daytime(maxR(t1,term_dates[5]),minR(t2,term_dates[6]));
      return times;
   }
  
   //Split period between t1 and t2 into subperiods corresponding to each term
   //with weekly adjustment for use with sine terms
   row_vector rebreak(real t1, real t2, int sun_adjust, row_vector term_dates) {
      row_vector[3] ptimes;
      row_vector[6] times;
      ptimes = p_times(t1,t2,term_dates);
      times = rep_row_vector(0.0,6);
      if(ptimes[1] > 0) {
         times[1] = maxR(t1,term_dates[1]) + sun_adjust;
         times[2] = minR(t2,term_dates[2]) + sun_adjust;
      }
      if(ptimes[2] > 0) {
         times[3] = maxR(t1,term_dates[3]) + sun_adjust;
         times[4] = minR(t2,term_dates[4]) + sun_adjust;
      }
      if(ptimes[3] > 0) {
         times[5] = maxR(t1,term_dates[5]) + sun_adjust;
         times[6] = minR(t2,term_dates[6]) + sun_adjust;    
      }
      return times;
   }
  
   //Calculate breaks that partitions full time interval into smaller intervals
   real[] calc_breaks(int N, real[] interactions, real[] durations, 
                        real[] active_times) {
      real brk_points[max(N+1,2)];
      brk_points[1] = active_times[1]-24; //Buffer
      if(N > 1) {
         for(n in 2:N) {
            brk_points[n] = (interactions[n-1] + durations[n-1]
                              + interactions[n])/2;
         }
      }
      brk_points[N+1] = active_times[2]+24;
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
  
   //Given parameter values, use forward-backward algorithm to estimate dyadic
   //edge at desired times. To be used in R after calling expose_stan_fucntions
   real[] infer_ships(int N, real[] interactions, real[] durations, 
                        real[] active_times, int same_floor, int same_year, 
                        int sun_adjust, row_vector term_dates, vector k_2, 
                        vector c_0, real c_1, vector c_2, real c_3, real k_0, 
                        vector k_1, vector k_3, real t_prob, real s_0, 
                        vector s) {
                           
      real ships[N]; //Vector of estimated network connection at various times
      
      real mult;
      real sparsity;
      
      real intensity[4];
      real integral[10]; //1+3*3
      real intensity_sum;
      real integral_sum;
      row_vector[2] probs_0;
      row_vector[2] probs_1;
      row_vector[2] intensities;
      row_vector[2] integrals;
      vector[3] baseline;
      int m;

      baseline[1] = 1;
      baseline[2] = c_0[1];
      baseline[3] = c_0[2];
      sparsity = (1+s[1]*same_floor)*(1+s[2]*same_year)*s_0;
    
      ships = rep_array(-1.0,max(N,1));
    
      if(N != 0) {
         real scale3h;
         real scale12h;
         real scalecos;
         real brk_points[N+1];
         real diffHours[N];
         row_vector[3] periodTimes[N];
         row_vector[3] dayTimes[N];
         real areDays[N];
         int periods[N];
         row_vector[6] scaledBreaks1[N];
         row_vector[6] scaledBreaks2[N];
         real scaledInts[N];
         
         row_vector[2] forward[N];
         row_vector[2] backward[N];
      
         scale3h = 1.0/3;
         scale12h = 1.0/12;
         scalecos = pi()/(7*12);  
         
         brk_points = rep_array(0.0,N+1);
         diffHours = rep_array(0.0,N);
         periodTimes = rep_array(rep_row_vector(0.0,3),N);
         dayTimes = rep_array(rep_row_vector(0.0,3),N);
         areDays = rep_array(0.0,N);
         periods = rep_array(0,N);
         scaledBreaks1 = rep_array(rep_row_vector(0.0,6),N);
         scaledBreaks2 = rep_array(rep_row_vector(0.0,6),N);
         scaledInts = rep_array(0.0,N);

         brk_points[1] = active_times[1]-24; //Buffer
      
         scaledInts[1] = scalecos*(interactions[1]+sun_adjust);
         for(n in 2:N) {
            brk_points[n] = (interactions[n-1] + interactions[n])/2;
            scaledInts[n] = scalecos*(interactions[n]+sun_adjust);
         }
         brk_points[N+1] = active_times[2]+24;
         diffHours[1] = (brk_points[2] - brk_points[1])/2;
         areDays[1] = scale12h*isDay(interactions[1]);
         periods[1] = period(interactions[1],term_dates);

         scaledBreaks1[1] = scalecos*rebreak(brk_points[1],
                              interactions[2],sun_adjust,term_dates);
         scaledBreaks2[1] = scalecos*rebreak(interactions[1]+durations[1],
                              brk_points[2],sun_adjust,term_dates);
         periodTimes[1] = p_times(brk_points[1],interactions[1],term_dates) 
            + p_times(interactions[1]+durations[1],brk_points[2],term_dates);
         dayTimes[1] = scale12h*
            (day_terms(brk_points[1],interactions[1],term_dates)
            + day_terms(interactions[1]+durations[1],brk_points[2],term_dates));
            
         if(N > 1) {
            for(n in 2:N) {
               diffHours[n] = (brk_points[n+1] - brk_points[n-1])/2;
               scaledBreaks1[n] = scalecos*rebreak(brk_points[n],
                                    interactions[n],sun_adjust,term_dates);
               scaledBreaks2[n] = scalecos*rebreak(interactions[n]+durations[n],
                                    brk_points[n+1],sun_adjust,term_dates);
               areDays[n] = scale12h*isDay(interactions[n]);
               periods[n] = period(interactions[n],term_dates);
               periodTimes[n] = 
                  p_times(brk_points[n],interactions[n],term_dates) 
                  + p_times(interactions[n]+durations[n],
                              brk_points[n+1],term_dates);
               dayTimes[n] = scale12h*
                  (day_terms(brk_points[n],interactions[n],term_dates)
                  + day_terms(interactions[n]+durations[n],brk_points[n+1],
                     term_dates));
            }
         }

         //Probability from 0 and from 1
         probs_0[1] = log_sum_exp(log1m(sparsity), 
                                    log(sparsity)-t_prob*diffHours[1]);
         probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[1]);
         probs_1[1] = log1m_exp(probs_0[1]);      
         probs_1[2] = log1m_exp(probs_0[2]);

         mult = 1+same_floor*c_1;

         intensity = rep_array(0.0,4);
         integral = rep_array(0.0,10);
         intensity[1] = scale3h*k_0;
         integral[1] = scale3h*k_0*(periodTimes[1]*baseline);
         for (i in 1:3) {
            intensity[i+1] = scale12h*k_1[i]*cos(k_2[i]*scaledInts[1]+k_3[i]);
            for(j in 1:3) {
               if(periodTimes[1,j] > 0) {
                  integral[i*3+j-2] = 7/pi()*baseline[j]*k_1[i]/k_2[i]*(
                  sin(k_2[i]*scaledBreaks1[1][2*j]+k_3[i]) - 
                  sin(k_2[i]*scaledBreaks1[1][2*j-1]+k_3[i]) + 
                  sin(k_2[i]*scaledBreaks2[1][2*j]+k_3[i]) - 
                  sin(k_2[i]*scaledBreaks2[1][2*j-1]+k_3[i]));
               }
            }
         }
         intensity_sum = sum(intensity)*baseline[periods[1]];
         integral_sum = sum(integral);

         forward[1][1] = log_sum_exp(log1m(sparsity)+probs_0[1],
                           log(sparsity)+probs_0[2]) + log(mult*intensity_sum) -
                           mult*integral_sum;
         forward[1][2] = log_sum_exp(log1m(sparsity)+probs_1[1],
                           log(sparsity)+probs_1[2]) + 
                           log(mult*(1+c_2[same_floor+1])*(intensity_sum +
                              c_3*areDays[1]*baseline[periods[1]])) -
                           mult*(1+c_2[same_floor+1])*
                              (integral_sum + c_3*dayTimes[1]*baseline);
        
         if(N > 1) {
            for (n in 2:N) {
               probs_0[1] = log_sum_exp(log1m(sparsity), 
                                          log(sparsity)-t_prob*diffHours[n]);
               probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[n]);
               probs_1[1] = log1m_exp(probs_0[1]);      
               probs_1[2] = log1m_exp(probs_0[2]);

               intensity = rep_array(0.0,4);
               integral = rep_array(0.0,10);
               intensity[1] = scale3h*k_0;
               integral[1] = scale3h*k_0*(periodTimes[n]*baseline);
               for (i in 1:3) {
                  intensity[i+1] = scale12h*k_1[i]*
                                    cos(k_2[i]*scaledInts[n]+k_3[i]);
                  for(j in 1:3) {
                     if(periodTimes[n,j] > 0) {
                        integral[i*3+j-2] = 7/pi()*baseline[j]*k_1[i]/k_2[i]*(
                        sin(k_2[i]*scaledBreaks1[n][2*j]+k_3[i]) - 
                        sin(k_2[i]*scaledBreaks1[n][2*j-1]+k_3[i]) +
                        sin(k_2[i]*scaledBreaks2[n][2*j]+k_3[i]) - 
                        sin(k_2[i]*scaledBreaks2[n][2*j-1]+k_3[i]));
                     }
                  }
               }
               intensity_sum = sum(intensity)*baseline[periods[n]];
               integral_sum = sum(integral);

               forward[n][1] = log_sum_exp(forward[n-1] + probs_0) + 
                                 log(mult*intensity_sum) - mult*integral_sum;
               forward[n][2] = log_sum_exp(forward[n-1] + probs_1) + 
                                 log(mult*(1+c_2[same_floor+1])*(intensity_sum +
                                    c_3*areDays[n]*baseline[periods[n]])) -
                                 mult*(1+c_2[same_floor+1])*
                                 (integral_sum + c_3*dayTimes[n]*baseline);
            }
         }

         backward[N][1] = 0;
         backward[N][2] = 0;

         if(N > 1) {
            for (n in 2:N) {
               m = N - n + 1;

               //Probabilities to transition from 0 and 1
               probs_0[1] = log_sum_exp(log1m(sparsity), 
                                          log(sparsity)-t_prob*diffHours[m]);
               probs_0[2] = log1m_exp(probs_0[1]);
               probs_1[1] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[m]);    
               probs_1[2] = log1m_exp(probs_1[1]);   
             
               intensity = rep_array(0.0,4);
               integral = rep_array(0.0,10);
               intensity[1] = scale3h*k_0;
               integral[1] = scale3h*k_0*periodTimes[m+1]*baseline;
               for (i in 1:3) {
                  intensity[i+1] = scale12h*k_1[i]*
                                    cos(k_2[i]*scaledInts[m+1]+k_3[i]);
                  for(j in 1:3) {
                     if(periodTimes[m+1,j] > 0) {
                        integral[i*3+j-2] = 7/pi()*baseline[j]*k_1[i]/k_2[i]*(
                        sin(k_2[i]*scaledBreaks1[m+1][2*j]+k_3[i]) - 
                        sin(k_2[i]*scaledBreaks1[m+1][2*j-1]+k_3[i]) +
                        sin(k_2[i]*scaledBreaks2[m+1][2*j]+k_3[i]) - 
                        sin(k_2[i]*scaledBreaks2[m+1][2*j-1]+k_3[i]));
                     }
                  }
               }
               intensity_sum = sum(intensity)*baseline[periods[m+1]];
               integral_sum = sum(integral);
                    
               intensities[1] = mult*intensity_sum;
               intensities[2] = (1+c_2[same_floor+1])*
                  (intensities[1]+mult*c_3*areDays[m+1]*baseline[periods[m+1]]);

               integrals[1] = mult*integral_sum;
               integrals[2] = (1+c_2[same_floor+1])*
                  (integrals[1] + mult*c_3*dayTimes[m+1]*baseline);

               backward[m][1] = log_sum_exp(probs_0 + backward[m+1] + 
                                             log(intensities) - integrals);
               backward[m][2] = log_sum_exp(probs_1 + backward[m+1] +
                                             log(intensities) - integrals);
            }
         }      
         for(n in 1:N) {
            ships[n] = inv_logit(forward[n][2] + backward[n][2] 
                                    - forward[n][1] - backward[n][1]);
         }

      //If no interactions, probability of a relationship is very low (~ zero) 
      } else {
         ships[1] = 0;
      }
      return ships;
   }
}

data {
   int D; //Number of dyads
   int<lower=0> N[D]; //Number of interactions for each dyad
   real<lower=0> interactions[D,max(N)]; //Interaction times (in hours)
   real<lower=0> durations[D,max(N)]; //Duration of interactions (in hours)
   real<lower=0> active_times[D,2]; //Active tiems for each dyad (in hours)
   int<lower=0,upper=1> same_floor[D]; //Indicator for same floor
   int<lower=0,upper=1> same_year[D]; //Indicator for same year
   //Adjustment to Sunday (beginning of week) - used for weekly periodicity
   int sun_adjust;
   row_vector[6] term_dates; //Important academic calendar times (in hours)
   positive_ordered[3] k_2;
}

transformed data {
   //Often used constants
   real<lower=0> scale3h;
   real<lower=0> scale12h;
   real<lower=0> scalecos;
   
   real<lower=-24> brk_points[D,max(N)+1]; //Break times
   real<lower=0> diffHours[D,max(N)]; //Difference between event bins in hours
   row_vector<lower=0>[3] periodTimes[D,max(N)]; //Term time between break times
   row_vector<lower=0>[3] dayTimes[D,max(N)]; //Day time between break times
   real<lower=0> areDays[D,max(N)]; //Is this event time at day?
   int<lower=0,upper=3> periods[D,max(N)]; //Which term did this event occur in?
   row_vector<lower=0>[6] scaledBreaks[D,max(N),2]; //Breaks scaled by scalecos
   real<lower=0> scaledInts[D,max(N)]; //Interaction times scaled by scalecos

   scale3h = 1.0/3;
   scale12h = 1.0/12;
   scalecos = pi()/(7*12);  

   //Initialize arrays
   brk_points = rep_array(0.0,D,max(N)+1);
   diffHours = rep_array(0.0,D,max(N));
   periodTimes = rep_array(rep_row_vector(0.0,3),D,max(N));
   dayTimes = rep_array(rep_row_vector(0.0,3),D,max(N));
   areDays = rep_array(0.0,D,max(N));
   periods = rep_array(0,D,max(N));
   scaledBreaks = rep_array(rep_row_vector(0.0,6),D,max(N),2);
   scaledInts = rep_array(0.0,D,max(N));

   //For each dyad
   for (d in 1:D) {
      brk_points[d,1] = active_times[d,1]-24; //Buffer
    
      if(N[d] != 0) {
         
         //Calculate breaks and scaled interactions
         scaledInts[d,1] = scalecos*(interactions[d,1]+sun_adjust);
         for(n in 2:(N[d]+1)) {
            if (n <= N[d]) {
               brk_points[d,n] = (interactions[d,n-1] + durations[d,n-1] 
                                 + interactions[d,n])/2;
               scaledInts[d,n] = scalecos*(interactions[d,n]+sun_adjust);
            } else {
               brk_points[d,n] = active_times[d,2]+24;
            }
         }

         //Calculate characteristics about interactions and intervals
         for(n in 1:N[d]) {
            if(n == 1) {
               diffHours[d,1] = (brk_points[d,2] - brk_points[d,1])/2;               
            } else {
               diffHours[d,n] = (brk_points[d,n+1] - brk_points[d,n-1])/2;               
            }
            //Ignore duration while interaction occurs
            scaledBreaks[d,n,1] = scalecos*
               rebreak(brk_points[d,n],interactions[d,n],sun_adjust,term_dates);
            scaledBreaks[d,n,2] = scalecos*
               rebreak(interactions[d,n] + durations[d,n],brk_points[d,n+1],
                        sun_adjust,term_dates);
            areDays[d,n] = scale12h*isDay(interactions[d,n]);
            periods[d,n] = period(interactions[d,n],term_dates);
            periodTimes[d,n] = p_times(brk_points[d,n],interactions[d,n],
                                       term_dates)
                                 + p_times(interactions[d,n] + durations[d,n],
                                             brk_points[d,n+1],term_dates);
            dayTimes[d,n] = scale12h*
               (day_terms(brk_points[d,n],interactions[d,n],term_dates)
                  + day_terms(interactions[d,n] + durations[d,n],
                              brk_points[d,n+1],term_dates));
         }
         
      //If no interactions
      } else {
         brk_points[d,2] = active_times[d,2]+24;
         scaledBreaks[d,1,1] = scalecos*rebreak(brk_points[d,1],brk_points[d,2],
                                                sun_adjust,term_dates);
         periodTimes[d,1] = p_times(brk_points[d,1],brk_points[d,2],term_dates);
      }
   }
}

parameters {
   vector<lower=0>[2] c_0; //Multipliers for term (first term implicitly = 1)
   real<lower=0> c_1; //Floor multiplier
   vector<lower=0>[2] c_2; //Friendship interactivity multiplier
   real<lower=0> c_3; //Daytime adjustment to baseline for friends
   real<lower=0> k_0; //Weekly baseline
   vector<lower=0>[3] k_1; //Sine weights
   vector<lower=0,upper=(2*pi())>[3] k_3; //Sine offsets
   real<lower=0> t_prob; //CTMC transition probability "q"
   real<lower=0,upper=1> s_0; //Base level of sparsity; CTMC initial probability
   vector<lower=0>[2] s; //Sparsity multipliers for floor/year
}

model {
   real log_obs[D]; //Vector of log-likelihoods

   row_vector[2] probs_0; //Log probabilities to transition to 0
   row_vector[2] probs_1; //Log probabilities to transition to 1
   
   vector[3] baseline; //Baseline intensity based on terms
      
   //Assorted placeholer variables
   real mult;
   real sparsity;
   real intensity[4];
   real integral[10];
   real intensity_sum;
   real integral_sum;

   //Priors
   c_0[1] ~ exponential(1);
   c_0[2] ~ exponential(1);
   c_1 ~ exponential(1);
   c_2[1] ~ exponential(100);
   c_2[2] ~ exponential(100);
   c_3 ~ exponential(1);
   k_0 ~ exponential(1);
   k_1[1] ~ exponential(1);
   k_1[2] ~ exponential(1);
   k_1[3] ~ exponential(1);
   t_prob ~ exponential(1e10);
   s_0 ~ beta(1,9);
   s[1] ~ exponential(1);
   s[2] ~ exponential(1);

   baseline[1] = 1;
   baseline[2] = c_0[1];
   baseline[3] = c_0[2];
  
   for(d in 1:D) {
      if(N[d] != 0) {
         //Forward variables from forward-backward algorithm         
         row_vector[2] forward[N[d]];
      
         //CTMC sparsity term
         sparsity = (1+s[1]*same_floor[d])*(1+s[2]*same_year[d])*s_0;
         //Multiplier due to floor
         mult = 1+same_floor[d]*c_1;
         
         //Transition probabilities between intervals
         probs_0[1] = log_sum_exp(log1m(sparsity), 
                                    log(sparsity)-t_prob*diffHours[d,1]);
         probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[d,1]);
         probs_1[1] = log1m_exp(probs_0[1]);      
         probs_1[2] = log1m_exp(probs_0[2]);


         //Calculate intensity and integral of intensity over period
         intensity = rep_array(0.0,4);
         integral = rep_array(0.0,10);
         intensity[1] = scale3h*k_0;
         integral[1] = scale3h*k_0*(periodTimes[d,1]*baseline);
         
         //Sine terms
         for (i in 1:3) {
            intensity[i+1] = scale12h*k_1[i]*cos(k_2[i]*scaledInts[d,1]+k_3[i]);
            for(j in 1:3) {
               if(periodTimes[d,1,j] > 0) {
                  integral[i*3+j-2] = 7/pi()*baseline[j]*k_1[i]/k_2[i]*(
                     sin(k_2[i]*scaledBreaks[d,1,1,2*j]+k_3[i]) -
                     sin(k_2[i]*scaledBreaks[d,1,1,2*j-1]+k_3[i]) +
                     sin(k_2[i]*scaledBreaks[d,1,2,2*j]+k_3[i]) -
                     sin(k_2[i]*scaledBreaks[d,1,2,2*j-1]+k_3[i]));
               }
            }
         }
         intensity_sum = sum(intensity)*baseline[periods[d,1]];
         integral_sum = sum(integral);

         //Calculate forward variables
         forward[1][1] = log_sum_exp(log1m(sparsity)+probs_0[1],
                           log(sparsity)+probs_0[2])
            + log(mult*intensity_sum) - mult*integral_sum;
         forward[1][2] = log_sum_exp(log1m(sparsity)+probs_1[1],
                           log(sparsity)+probs_1[2])
            + log(mult*(1+c_2[same_floor[d]+1])*
               (intensity_sum + c_3*areDays[d,1]*baseline[periods[d,1]]))
            - mult*(1+c_2[same_floor[d]+1])*
               (integral_sum + c_3*dayTimes[d,1]*baseline);
	  
         if(N[d] > 1) {
            for (n in 2:N[d]) {
               probs_0[1] = log_sum_exp(log1m(sparsity), 
                              log(sparsity)-t_prob*diffHours[d,n]);
               probs_0[2] = log1m(sparsity) + log1m_exp(-t_prob*diffHours[d,n]);
               probs_1[1] = log1m_exp(probs_0[1]);      
               probs_1[2] = log1m_exp(probs_0[2]);

               intensity = rep_array(0.0,4);
               integral = rep_array(0.0,10);
               intensity[1] = scale3h*k_0;
               integral[1] = scale3h*k_0*(periodTimes[d,n]*baseline);
               for (i in 1:3) {
                  intensity[i+1] = scale12h*k_1[i]
                                    *cos(k_2[i]*scaledInts[d,n]+k_3[i]);
                  for(j in 1:3) {
                     if(periodTimes[d,n,j] > 0) {
                        integral[i*3+j-2] = 7/pi()*baseline[j]*k_1[i]/k_2[i]*(
                        sin(k_2[i]*scaledBreaks[d,n,1,2*j]+k_3[i]) -
                        sin(k_2[i]*scaledBreaks[d,n,1,2*j-1]+k_3[i]) +
                        sin(k_2[i]*scaledBreaks[d,n,2,2*j]+k_3[i]) -
                        sin(k_2[i]*scaledBreaks[d,n,2,2*j-1]+k_3[i]));
                     }
                  }
               }
               intensity_sum = sum(intensity)*baseline[periods[d,n]];
               integral_sum = sum(integral);

               forward[n][1] = log_sum_exp(forward[n-1] + probs_0) + 
                                 log(mult*intensity_sum) - mult*integral_sum;
               forward[n][2] = log_sum_exp(forward[n-1] + probs_1) + 
                  log(mult*(1+c_2[same_floor[d]+1])
                     *(intensity_sum + c_3*areDays[d,n]*baseline[periods[d,n]]))
                  - mult*(1+c_2[same_floor[d]+1])
                     *(integral_sum + c_3*dayTimes[d,n]*baseline);
            }
         }
         log_obs[d] = log_sum_exp(forward[N[d]]);
         
      } else {
         //Assumes if no observations then they are not friends         
         integral = rep_array(0.0,10);
         integral[1] = scale3h*k_0*(periodTimes[d,1]*baseline);
         for (i in 1:3) {
            for(j in 1:3) {
               if(periodTimes[d,1,j] > 0) {
                  integral[i*3+j-2] = 7/pi()*baseline[j]*k_1[i]/k_2[i]*
                     (sin(k_2[i]*scaledBreaks[d,1,1,2*j]+k_3[i]) +
                     sin(k_2[i]*scaledBreaks[d,1,1,2*j-1]+k_3[i]));
               }
            }
         }
         mult = 1+same_floor[d]*c_1;
         log_obs[d] = -mult*sum(integral);
      }
   }
 
   //Posterior
   target += sum(log_obs);
}
