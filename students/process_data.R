# 1. Process/clean raw data
# 2. Convert to form usable by Stan




#Read-in raw data
#################

#Restrict users: Must be freshman to senior and report floor of residence
subjects <- read.csv("raw_data/Subjects.csv")
subjects[subjects == ""] = NA

#Reported relastionship information
zz = bzfile("raw_data/RelationshipsFromSurveys.csv.bz2")
reported_ships <- read.csv(zz)
survey_dates = unique(reported_ships$survey.date)
r_ship_types = unique(reported_ships$relationship)

#Restricted users
r_users = subjects$user_id
restrict = 
   (subjects$year_school %in% c("Freshman","Sophomore","Junior","Senior")) & 
   (subjects$floor %in% levels(subjects$floor)[2:9])
r_users = r_users[r_users %in% which(restrict)]





as.time = function(t) as.POSIXct(t, origin = "1970-01-01")

#Load in proximity data
zz = bzfile("raw_data/Proximity.csv.bz2")
proximity = read.csv(zz)

year = subjects$year_school
floor = subjects$floor

term_dates = as.POSIXct(c("2008-09-03","2008-12-13","2009-01-05",
                          "2009-03-20","2009-03-30","2009-05-15"))


#Restrict to proximity observations with prob. phones on same floor > 0.2
proximity = proximity[is.na(proximity$prob2) | proximity$prob2 > 0.2,]
proximity$time = as.time(proximity$time)
#Restrict to observations during the active school year
proximity = proximity[
   (proximity$time > term_dates[1]) & (proximity$time < term_dates[2]) |
   (proximity$time > term_dates[3]) & (proximity$time < term_dates[4]) |
   (proximity$time > term_dates[5]) & (proximity$time < term_dates[6]),]



#Restrict to users with one or more proximity entries
int_counts <- sapply(r_users, function(r) 
   sum(proximity$user.id == r | proximity$remote.user.id.if.known == r))
r_users = r_users[int_counts != 0]
indices = t(combn(r_users,2))



#Calculate time intervals over which each user is active
time_intervals = matrix(data = 0, nrow = length(r_users), ncol = 2)
time_intervals[,1] = sapply(1:length(r_users),function(u) 
   min(proximity[proximity$user.id == r_users[u] | 
         proximity$remote.user.id.if.known == r_users[u],]$time))
time_intervals[,2] = sapply(1:length(r_users),function(u) 
   max(proximity[proximity$user.id == r_users[u] | 
         proximity$remote.user.id.if.known == r_users[u],]$time))



#Reduce data to relevant users
proximity = proximity[proximity$user.id %in% r_users,]
proximity = proximity[proximity$remote.user.id.if.known %in% r_users,]
proximity = proximity[proximity$user.id != proximity$remote.user.id.if.known,]

#Series of function that remove duplicate interactions by merging interactions
# within 30 min; also deals with overlapping and unreciprocated interactions
expand_secs <- function(i, j) {
   if(sum(proximity$user.id == j & proximity$remote.user.id.if.known == i)>0) {
      #Pings every 6 minutes
      t1 = proximity$time[proximity$user.id == j & 
                          proximity$remote.user.id.if.known == i]-30*60/2
      t2 = proximity$time[proximity$user.id == j & 
                          proximity$remote.user.id.if.known == i]+30*60/2
      as.vector(sapply(1:length(t1), function(x) seq(t1[x],t2[x],by=1)))
   }
}

thin <- function(time_seq) {
   if(sum(diff(time_seq)>1) > 0) {
      end_points <- (diff(time_seq) > 1)*(1:(length(time_seq)-1))
      end_points <- c(end_points[end_points != 0],length(time_seq))
      start_points <- c(1,end_points[1:(length(end_points)-1)]+1)
      cbind(time_seq[start_points],time_seq[end_points])
   } else {
      cbind(time_seq[1],time_seq[length(time_seq)])
   }
}

merge <- function(i,j) {
   times1 <- expand_secs(i,j)
   times2 <- expand_secs(j,i)
   if(is.null(times1) && !is.null(times2)) {
      thin(times2)
   } else if(!is.null(times1) && is.null(times2)) {
      thin(times1)
   } else if(!is.null(times1) && !is.null(times2)) {
      thin(sort(c(times1,times2)))
   }
}

#New dataset after cleaning
proximity2 <- apply(indices,1,function(x) merge(x[1],x[2]))
durations <- sapply(proximity2, function(x) (x[,2] - x[,1])/(60*60))
proximity2 <- sapply(proximity2, function(x) x[,1])
N = sapply(proximity2,length)

#Calculate interval over which each dyad is active
calc_dyad_active <- function(index) {
   interval = c(max(time_intervals[r_users == index[1],1],
                    time_intervals[r_users == index[2],1]),
                min(time_intervals[r_users == index[1],2],
                    time_intervals[r_users == index[2],2]))
   interval[2] = max(interval) #In case of non-overlapping time_intervals
   interval
}
active_dyads <- t(apply(indices,1,calc_dyad_active))





#Calculate some baseline parameters values

#Aggregated interactions
ptimes = as.POSIXlt(unlist(proximity2), origin = "1970-01-01")
weekly = ptimes$wday*48+ptimes$hour*2+(ptimes$min>=30)

#Interaction patterns in a week (by half-hour)
w_counts = table(weekly)





#Process data for use by Stan
#############################

start = term_dates[1]
D = length(N)
ints_mat = matrix(data = 0, nrow = D, ncol = max(N))
dur_mat = matrix(data = 0, nrow = D, ncol = max(N))
active_mat = matrix(data = 0, nrow = D, ncol = 2)
for(i in 1:D) {
   active_mat[i,] = 
      as.numeric(difftime(as.time(active_dyads[i,]),start,unit="hours"))
   if(N[i] != 0) {
      ints_mat[i,1:N[i]] = 
         as.numeric(difftime(as.time(proximity2[[i]]),start,unit="hours"))
      dur_mat[i,1:N[i]] = durations[[i]]
   }
}

same_floor <- apply(indices,1,function(x) floor[x[1]] == floor[x[2]])
same_year <- apply(indices,1,function(x) year[x[1]] == year[x[2]])

#Estimate patterns of weekly periodicity
weekly_vol = as.numeric(table(weekly))
weekly_vol = weekly_vol/mean(weekly_vol)
k_2 = which(abs(fft(weekly_vol))[1:168] %in% tail(sort(abs(fft(weekly_vol))[2:168]),n = 3)) - 1
cs = -Re(sqrt(as.complex(-1))*log(fft(weekly_vol)[k_2+1]/abs(fft(weekly_vol)[k_2+1])))


#Lists used for Stan
data_list <- list(D = D, N = N, interactions = ints_mat, durations = dur_mat, 
                  active_times = active_mat, same_floor = same_floor, 
                  same_year = same_year, 
                  sun_adjust = as.POSIXlt(term_dates[1])$wday*24, 
                  term_dates = as.numeric(difftime(term_dates,start,
                                 unit = "hours") + c(0,24,-24,24,-24,24)),  
                  k_2 = k_2)
param_init <- list(c_0 = c(1.5,2), c_1 = 1, c_2 = c(1,1), c_3 = 1/10, k_0 = 0.1, 
                   k_1 = 2*abs(fft(weekly_vol)[k_2+1])/336*12/30, 
                   k_3 = cs %% (2*pi), t_prob = 1e-10, s_0 = 0.1, s = c(1,1))





#Construct reported friendship using desired measure
calc_reported_ships <- function(types, reciprocal) {
   reported_net <- lapply(1:length(data_list$N), function(i) rep(0,5))
   for(t in 1:5) {
      edges = reported_ships[reported_ships$survey.date == survey_dates[t] & 
               (reported_ships$relationship %in% r_ship_types[types]),1:2]
      
      #Non-reciprocal
      if(reciprocal) {
         e_indices = sapply(1:dim(indices)[1], function(i) 
            (sum(edges[,1] == indices[i,1] & edges[,2] == indices[i,2]) > 0) & 
            (sum(edges[,2] == indices[i,1] & edges[,1] == indices[i,2]) > 0)) 
      } else {
         e_indices = sapply(1:dim(indices)[1], function(i) 
            sum(edges[,1] == indices[i,1] & edges[,2] == indices[i,2]) + 
            sum(edges[,2] == indices[i,1] & edges[,1] == indices[i,2]) > 0)      
      }
      
      edges = indices[e_indices,]
      
      for(e in 1:dim(edges)[1]) {
         index = which(indices[,1] == edges[e,1] & indices[,2] == edges[e,2])
         if(length(index)) {
            reported_net[[index]][t] = 1
         }
      }
   }
   reported_net  
}

#We will use (possibly unreciprocated) "close friends" as our measure of
#reported friendship
reported_net <- calc_reported_ships(1,FALSE)




#Save processing
save(data_list,param_init,proximity2,weekly_vol,indices,r_users,time_intervals,
     reported_net,weekly,floor,start,survey_dates, 
     file="processed_data/stan_pre.RData")
