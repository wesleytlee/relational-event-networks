# 1. Process raw data





#Read-in raw data
#################

tag_ids <- c(2,5,6,8,9,10,11,12,15,17,18,19,23,51,52,53,55)
N = length(tag_ids)

raw_data <- list()
for(i in 1:N) {
  raw_data[[i]] <- read.table(
    paste("raw_data/140726_115712 encountered logs for ID ",
          tag_ids[i],".txt",sep=""),
    sep = "\t", header = TRUE)  
}

covariates <- read.csv("raw_data/Levin attribute file.csv", 
                       stringsAsFactors = FALSE)
covariates$Tag <- as.numeric(substring(covariates$Tag,4))

data <- rbind(raw_data[[1]],raw_data[[2]],raw_data[[3]],raw_data[[4]],
              raw_data[[5]],raw_data[[6]],raw_data[[7]],raw_data[[8]],
              raw_data[[9]],raw_data[[10]],raw_data[[11]],raw_data[[12]],
              raw_data[[13]],raw_data[[14]],raw_data[[15]],raw_data[[16]])
data$first.time <- strptime(data$first.time, format = "%Y/%m/%d %H:%M:%S")
data$last.time <- strptime(data$last.time, format = "%Y/%m/%d %H:%M:%S")










#Data cleaning
##############

#Restrict to interactions between swallows
data = data[data$this.ID %in% tag_ids,]
data = data[data$enc.ID %in% tag_ids,]

#Calculate periods of activity
indices = t(combn(tag_ids,2))
active_periods <- sapply(tag_ids, function(i) 
  c(min(data$first.time[data$this.ID == i | data$enc.ID == i]),
    max(data$last.time[data$this.ID == i | data$enc.ID == i])))
dyad_active <- function(index) {
  c(max(active_periods[tag_ids == index[1]][[1]][1],
        active_periods[tag_ids == index[2]][[1]][1]),
    min(active_periods[tag_ids == index[1]][[1]][2],
        active_periods[tag_ids == index[2]][[1]][2]))
}
active_dyads <- apply(indices, 1, dyad_active)



#Merge overlapping interactions (within 30 seconds)
expand_secs <- function(i, j) {
  if(sum(data$this.ID == j & data$enc.ID == i) > 0) {
    t1 = data[data$this.ID == j & data$enc.ID == i,]$first.time-15
    t2 = data[data$this.ID == j & data$enc.ID == i,]$last.time+15
    unlist(sapply(1:length(t1), function(x) seq(t1[x],t2[x],by=1)))
  }
}
merged <- function(i,j) {
  times1 <- expand_secs(i,j)
  times2 <- expand_secs(j,i)
  if(is.null(times1) && !is.null(times2)) {
    times2
  } else if(!is.null(times1) && is.null(times2)) {
    times1
  } else if(!is.null(times1) && !is.null(times2)) {
    sort(c(times1,times2))
  }
}
times <- apply(indices,1,function(x) merged(x[1],x[2]))



#Thin to intervals, then check RSSI
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
times <- sapply(times, thin)
duration <- sapply(times, function(x) (x[,2] - x[,1])/(60*60))
rssi_max <- sapply(sapply(times,length)/2, function(x) rep(0,x))
data = data[data$RSSI.max > 0,]
for(i in 1:dim(data)[1]) {
  index = which(rowSums(t(apply(indices,1,function(x) x == range(data[i,1:2]))))
                == 2)
  index2 = which((times[[index]][,1] <= data[i,5]) &
                   (times[[index]][,2] >= data[i,5]))
  rssi_max[[index]][index2] = 1
}
duration <- sapply(1:136, function(i) 
  if(length(times[[i]])) {duration[[i]][rssi_max[[i]] == 1]})
times <- sapply(1:136, function(i) 
  if(length(times[[i]])) {as.POSIXct(times[[i]][rssi_max[[i]] == 1,1],
                                     origin = "1970-01-01")})



#Restrict to observations within designated periods
restrict <- function(time_seq) {
  if(length(time_seq) > 0) {
    (as.POSIXlt(time_seq)$hour >= 6 & as.POSIXlt(time_seq)$hour <= 8) | 
      (as.POSIXlt(time_seq)$hour >= 17 & as.POSIXlt(time_seq)$hour <= 19)
  }
}
subset <- sapply(times, restrict)
times <- sapply(1:136, function(i) times[[i]][subset[[i]]])
duration <- sapply(1:136, function(i) duration[[i]][subset[[i]]])

save(times,duration,tag_ids,active_periods,active_dyads,covariates,indices, 
     file = "processed_data/pre_swallow.RData")
