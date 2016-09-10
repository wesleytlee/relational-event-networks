# 1. Create assorted figures and animations for swallow data





#Save figures as .eps (TRUE) or .pdf (FALSE)
#Issue with grid_arrange_shared_legend() and .pdf format
printEPS = TRUE

#Load in requisite libraries and processed daa
##############################################

if(printEPS) {
   #Requires package "Cairo"
   setHook(packageEvent("grDevices", "onLoad"),
           function(...) grDevices::X11.options(type='cairo'))
   options(device='x11')
}

library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)
library(igraph)
library(rstan)
library(Matrix)

load("processed_data/pre_swallow.RData")
load("processed_data/stan_swallow_simp.RData")

pre_load = stanc(file = "swallow.stan")
barn_model = stan_model(stanc_ret = pre_load)
expose_stan_functions(barn_model)

start = strptime("2014-07-19", format = "%Y-%m-%d")





#Descriptive plots
##################

#Multi-plot function from https://stackoverflow.com/questions/30611474/
grid_arrange_shared_legend <- function(...) {
   plots <- list(...)
   g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
   lheight <- sum(legend$height)
   p <- grid.arrange(
      #do.call(arrangeGrob, lapply(plots, function(x)
      #  x + theme(legend.position="none"))),
      arrangeGrob(plots[[1]] + theme(legend.position="none"),plots[[2]] 
                  + theme(legend.position="none"),nrow = 1),
      legend, nrow = 2, heights = unit.c(unit(1, "npc") - lheight, lheight))
   return(p)
}

#Histogram of Interaction Counts per Dyad
gender <- as.factor(sapply(data_list$gender, function(i) 
   if(i == 0) {"FF"} else if(i == 1) {"MF"} else {"MM"}))
gender.cols <- brewer.pal(3, "Set1")[c(1,3,2)]

p1 <- ggplot(NULL, aes(data_list$N, fill = gender)) +
   geom_histogram(binwidth = 2) +
   scale_fill_manual(values = gender.cols, 
                     guide = guide_legend(title = "Sex Pairing")) +
   ggtitle("Interaction Counts Per Dyad") +
   ylab("Number of Dyads") +
   xlab("Interaction Counts") +
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12))

#Histogram of Interactions over Time
#Desire slightly non-linear scale for better use of space
map_times <- function(t, skip) {
  if(t >= 102) {
    (3+skip)*7 + 1 + t - 102
  } else if(t >= 89) {
    (3+skip)*6 + 1 + t - 89
  } else if(t >= 78) {
    (3+skip)*5 + 1 + t - 78
  } else if(t >= 65) {
    (3+skip)*4 + 1 + t - 65
  } else if(t >= 54) {
    (3+skip)*3 + 1 + t - 54
  } else if(t >= 41) {
    (3+skip)*2 + 1 + t - 41
  } else if(t >= 30) {
    (3+skip)*1 + 1 + t - 30
  } else {
    1 + t - 17
  }
}

all_times <- difftime(as.POSIXlt(unlist(times), origin = "1970-01-01"),
                      start, unit = "hours")
all_gender <- as.factor(unlist(sapply(1:data_list$D, function(i) 
                                 rep(as.character(gender[i]),data_list$N[i]))))
all_dyads <- as.factor(unlist(sapply(1:data_list$D, function(i) 
                                 rep(i,data_list$N[i]))))
n_dyads_active <- rowSums(sapply(1:data_list$D, function(i) 
   p_time(data_list$active_times[all_dyads[i],1],
          data_list$active_times[all_dyads[i],2])))/3
all_weights <- sapply(1:length(all_times), function(i) 
                        n_dyads_active[period(all_times[i])])
all_times <- sapply(floor(all_times), function(i) map_times(i,2))

p2 <- ggplot(NULL, aes(all_times, fill = all_gender)) +
         geom_histogram(aes(weight = 1/all_weights), binwidth = 1) +
         scale_fill_manual(values = gender.cols, 
            guide = guide_legend(title = "Sex Pairing")) +
         ggtitle("Interactions Over Time") +
         ylab("Avg. Interactions Per Dyad") +
         scale_x_continuous(breaks = seq(5,35,10), 
            labels = c("July 20","July 21","July 22","July 23"), name = "Time") +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

if(printEPS) {
   cairo_ps("figures/fig8_descriptive_swallow.eps", 
            width=9, height=3, fallback_resolution = 800)
} else {
   pdf("figures/fig8_descriptive_swallow.pdf", width=9, height=3)
}
p <- grid_arrange_shared_legend(p1, p2)
dev.off()





#Trace plots
############

p1 <- ggplot(NULL, aes(x=1:9000, y=est_params$c[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(c[FF])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p2 <- ggplot(NULL, aes(x=1:9000, y=est_params$c[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(c[MF])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p3 <- ggplot(NULL, aes(x=1:9000, y=est_params$c[,3])) + geom_line() +
         xlab("Iterations") + ylab(expression(c[MM])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p4 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[1])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p5 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,2])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[2])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p6 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,3])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[3])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p7 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,4])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[4])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p8 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,5])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[5])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p9 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,6])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[6])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p10 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,7])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[7])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p11 <- ggplot(NULL, aes(x=1:9000, y=est_params$k[,8])) + geom_line() +
         xlab("Iterations") + ylab(expression(k[8])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p12 <- ggplot(NULL, aes(x=1:9000, y=est_params$sparsity)) + geom_line() +
         xlab("Iterations") + ylab(expression(s)) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p13 <- ggplot(NULL, aes(x=1:9000, y=est_params$t_prob)) + geom_line() +
         xlab("Iterations") + ylab(expression(q)) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

if(printEPS) {
   cairo_ps("figures/fig14_swallow_traces.eps", 
            width=9, height=10, fallback_resolution = 800)
} else {
   pdf("figures/fig14_swallow_traces.pdf", width=9, height=10)
}
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,ncol = 3)
dev.off()





#Network layouts
################

#Circular network layout
#Seperates swallows by gender and order males by tail streamer length
circ_order <- invPerm(c(which(covariates$Sex == "M")
                      [order(covariates$MeanTS[covariates$Sex == "M"])], 
                which(covariates$Sex == "F")))
m <- matrix(data = 0, nrow = length(tag_ids), ncol = length(tag_ids))
colnames(m) <- tag_ids
rownames(m) <- tag_ids
net <- graph.adjacency(m,mode = "undirected", weighted = TRUE, diag = FALSE)
circ_layout = layout.circle(net)[circ_order,]

#FR layout
m <- matrix(data = 0, nrow = length(tag_ids), ncol = length(tag_ids))
colnames(m) <- tag_ids
rownames(m) <- tag_ids
m[lower.tri(m, diag = FALSE)] <- dyad_friends
m <- m + t(m)
m[m < 0.1] <- 0
net <- graph.adjacency(m,mode = "undirected", weighted = TRUE, diag = FALSE)

set.seed(3)
fr_layout = layout.fruchterman.reingold(net)
fr_layout[17,] = fr_layout[17,] + c(0,0.5)






#Daily snapshots
################

#Daily snapshots
days <- function(t) {
   1*(t < start + 24*60*60) + 2*(t > start + 24*60*60 & t < start + 48*60*60) +
      3*(t > start + 48*60*60 & t < start + 72*60*60) +
      4*(t > start + 72*60*60 & t < start + 96*60*60) +
      5*(t > start + 96*60*60 & t < start + 108*60*60)
}

int_counts_by_day <- t(sapply(times, function(t) 
   table(factor(as.character(days(t)), levels = 1:5))))
day_names <- c("7/19", "7/20", "7/21", "7/22", "7/23")

nodes.cols = covariates$Sex
nodes.cols = gsub("F",rgb(1,0,0,0.2),nodes.cols)
nodes.cols = gsub("M",rgb(0,0,1,0.4),nodes.cols)

if(printEPS) {
   cairo_ps("figures/fig1_daily_snapshots.eps", 
            width=9, height=2, fallback_resolution = 800)
} else {
   pdf("figures/fig1_daily_snapshots2.pdf", width=9, height=2)
}
par(mar = c(0,0.5,1,0.5), mfrow = c(1,5))
for(day in 1:5) {
    m <- matrix(data = 0, nrow = length(tag_ids), ncol = length(tag_ids))
    m[lower.tri(m, diag = FALSE)] <- int_counts_by_day[,day]
    m <- m + t(m)
    net <- graph.adjacency(m,mode = "undirected", weighted = TRUE, diag = FALSE)
    V(net)$color = nodes.cols
    V(net)$color[sapply(active_periods, function(x) days(x[1])) > day |
    sapply(active_periods, function(x) days(x[2])) < day] = rgb(1,1,1,0.4)
    plot.igraph(net,vertex.label="",
                layout=circ_layout,edge.width=E(net)$weight/3, 
                main = day_names[day], vertex.size = 20)
}
dev.off()






#Function to networks at various times
######################################
est_snapshots <- function(time_hrs) {
   estd_ships <- t(sapply(1:data_list$D, function(i) rep(NA,length(time_hrs))))
   for(d in 1:data_list$D) {
      if(data_list$N[d] > 1) {
         for(t in 1:length(time_hrs)) {
            #If outside active interval, NA
            if(time_hrs[t] < data_list$active_times[d,1] | 
               time_hrs[t] > data_list$active_times[d,2]) {
               estd_ships[d,t] = NA
            
            #If time is before/after first/last midpoint, take midpoint value
            }else if(time_hrs[t] <= ship_mids[[d]][1]) {
               estd_ships[d,t] = sample_ships[d,1]
               
            }else if(time_hrs[t] >= ship_mids[[d]][data_list$N[d]]) {
               estd_ships[d,t] = sample_ships[d,data_list$N[d]]
            
            #Otherwise interpolate
            } else {
               ob = which(time_hrs[t] <= ship_mids[[d]][2:data_list$N[d]] & 
                           time_hrs[t] > ship_mids[[d]][2:data_list$N[d]-1])
               weight = (ship_mids[[d]][ob+1] - time_hrs[t])/
                           (ship_mids[[d]][ob+1] - ship_mids[[d]][ob])
               estd_ships[d,t] = weight*sample_ships[d,ob] + 
                                    (1-weight)*sample_ships[d,ob+1]
            }
         } 
      } else if (data_list$N[d] == 1) {
         for(t in 1:length(time_hrs)) {
            if(time_hrs[t] >= data_list$active_times[d,1] &
               time_hrs[t] <= data_list$active_times[d,2]) {
                  estd_ships[d,t] = sample_ships[d,1]
            }
         }
      } else {
         for(t in 1:length(time_hrs)) {
            if(time_hrs[t] >= data_list$active_times[d,1] & 
               time_hrs[t] <= data_list$active_times[d,2]) {
               estd_ships[d,t] = 0
            }
         }
      }
   }
   estd_ships
}





#Cluster paths by gender
########################
time_hrs = seq(17,105,0.1)
estd_ships <- est_snapshots(time_hrs)
network_paths <- data.frame(time <- rep(time_hrs*60*60+start,136), 
                            probs <- as.numeric(t(estd_ships)), 
                            dyad <- as.numeric(sapply(1:(data_list$D), 
                                       function(i) rep(i,length(time_hrs)))))
colnames(network_paths) <- c("time","probs","dyads")
network_paths$dyads <- as.factor(network_paths$dyads)

p1 <- ggplot(network_paths[network_paths$dyads %in% which(gender == "FF"),], 
             aes(x=time,y=probs,group=dyads)) +
   geom_line(colour = gender.cols[1],size = 1, alpha = 0.25) +
   ggtitle("Female/Female Pairs") +
   ylab("Probability of a Social Relation") +
   xlab("Time") +
   xlim(17,105) +
   ylim(0,1) +
   scale_x_datetime() +
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12),
         legend.position="none")

p2 <- ggplot(network_paths[network_paths$dyads %in% which(gender == "MF"),], 
             aes(x=time,y=probs,group=dyads)) +
   geom_line(colour = gender.cols[2],size = 1, alpha = 0.25) +
   ggtitle("Male/Female Pairs") +
   ylab("Probability of a Social Relation") +
   xlab("Time") +
   xlim(17,105) +
   ylim(0,1) +
   scale_x_datetime() +
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12),
         legend.position="none")

p3 <- ggplot(network_paths[network_paths$dyads %in% which(gender == "MM"),], 
             aes(x=time,y=probs,group=dyads)) +
   geom_line(colour = gender.cols[3],size = 1, alpha = 0.25) +
   ggtitle("Male/Male Pairs") +
   ylab("Probability of a Social Relation") +
   xlab("Time") +
   xlim(17,105) +
   ylim(0,1) +
   scale_x_datetime() + 
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12),
         legend.position="none")

if(printEPS) {
   cairo_ps("figures/fig9_swallow_paths.eps", 
            width=9, height=3, fallback_resolution = 800)
} else {
   pdf("figures/fig9_swallow_paths.pdf", width=9, height=3)
}
grid.arrange(p1,p2,p3,nrow = 1)
dev.off()





#Sexual Selection - Tail Streamer Length in Males
#################################################
if(printEPS) {
   cairo_ps("figures/fig11_swallow_selection.eps", 
            width=9, height=3, fallback_resolution = 800)
} else {
   pdf("figures/fig11_swallow_selection.pdf", width=9, height=3)
}

m <- matrix(c(1,2,3),nrow = 1,ncol = 3)
layout(mat = m,widths = c(0.45,0.45,0.1))

m <- matrix(data = 0, nrow = length(tag_ids), ncol = length(tag_ids))
colnames(m) <- tag_ids
rownames(m) <- tag_ids
m[lower.tri(m, diag = FALSE)] <- ifelse(gender == "MF",dyad_friends,0)
m <- m + t(m)
m[m < 0.1] <- 0
net <- graph.adjacency(m,mode = "undirected", weighted = TRUE, diag = FALSE)
V(net)$sex = covariates$Sex
V(net)$color = V(net)$sex
V(net)$color = gsub("F",rgb(1,0,0,0.2),V(net)$color)

#Tail streamer length
for(s in which(covariates$Sex == "M")) {
   if(is.na(covariates$MeanTS[s])) {
      V(net)$color[s] = rgb(0,0,0,1)
   } else {
      V(net)$color[s] = rgb(0,0,1,(covariates$MeanTS[s]-80)/20)   
   }
}
V(net)$label.color = "black"
V(net)$label.color[is.na(covariates$MeanTS)] = rgb(1,1,1)
par(mar = c(0,0,1,0))
plot.igraph(net,vertex.label=V(net)$name,layout=circ_layout, 
            edge.color = "grey", edge.width=E(net)$weight*5, 
            main = "Tail Streamer Length", vertex.size = 30)

#Ventral plumage coloration
for(s in which(covariates$Sex == "M")) {
  V(net)$color[s] = rgb(0,0,1,(-covariates$Color[s]+2)/7)   
}
V(net)$label.color = "black"
par(mar = c(0,0,1,0))
plot.igraph(net,vertex.label=V(net)$name,layout=circ_layout, 
            edge.color = "grey", edge.width=E(net)$weight*5, 
            main = "Ventral Plumage Coloration", vertex.size = 30)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="center",inset = 0, title = "Estimated\nProbability\nof a Relation",
       legend = c("0.25","0.5","0.75","1.0"), 
       col="grey", lwd=5*c(0.25,0.5,0.75,1), 
       cex=1, horiz = FALSE, bty = "n")
dev.off()





#Sub-network of male swallows
#############################
if(printEPS) {
   cairo_ps("figures/fig10_swallow_males.eps", 
            width=5, height=3, fallback_resolution = 800)
} else {
   pdf("figures/fig10_swallow_males.pdf", width=5, height=3)
}

par(mar = c(0,0,1,0), mfrow = c(1,1))
m <- matrix(data = 0, nrow = length(tag_ids), ncol = length(tag_ids))
m[lower.tri(m, diag = FALSE)] <- dyad_friends
colnames(m) <- tag_ids
rownames(m) <- tag_ids
m <- m + t(m)
m[m < 0.8] <- 0
net <- graph.adjacency(m[covariates$Sex == "M",covariates$Sex == "M"],
                        mode = "undirected", weighted = TRUE, diag = FALSE)
V(net)$color = rgb(0,0,1,0.4)
V(net)$label.color = "black"
set.seed(2)
layout <- layout.fruchterman.reingold(net)
layout[4,] = layout[4,] - c(4,4)
plot.igraph(net,vertex.label=V(net)$name,layout=layout,
            edge.width=E(net)$weight*5, main = "", vertex.size = 20)
dev.off()





#Network movie
##############
nodes.cols.bold = covariates$Sex
nodes.cols.bold = gsub("F",rgb(1,0,0),nodes.cols.bold)
nodes.cols.bold = gsub("M",rgb(0,0,1),nodes.cols.bold)

break_times = strptime(c("2014-07-19 17:00:00","2014-07-19 20:00:00",
                         "2014-07-20 6:00:00","2014-07-20 9:00:00",
                         "2014-07-20 17:00:00","2014-07-20 20:00:00",
                         "2014-07-21 6:00:00","2014-07-21 9:00:00",
                         "2014-07-21 17:00:00","2014-07-21 20:00:00",
                         "2014-07-22 6:00:00","2014-07-22 9:00:00",
                         "2014-07-22 17:00:00","2014-07-22 20:00:00",
                         "2014-07-23 6:00:00","2014-07-23 9:00:00",
                         "2014-07-24 23:00:00"),format="%Y-%m-%d %H:%M:%S")
breaks_in_hours <- as.numeric(difftime(break_times,start,unit="hours"))

#Remove hours outside of active intervals
merge_hours <- function(t) {
   rowSums(sapply(seq(1,15,2), function(i) 
      (t>=breaks_in_hours[i] & t<breaks_in_hours[i+2])*
      ((t-breaks_in_hours[i]>=3)*3+(t-breaks_in_hours[i]<3)*
          (t-breaks_in_hours[i])+3/2*i-3/2)))
}

active_pds_hrs = matrix(data = 0, nrow = length(tag_ids), ncol = 2)
for(i in 1:length(tag_ids)) {
   active_pds_hrs[i,] = 
      as.numeric(difftime(active_periods[[i]],start,unit="hours"))
}

active_consec <- cbind(merge_hours(active_pds_hrs[,1]),
                       merge_hours(active_pds_hrs[,2]))
ints_consec <- t(sapply(1:(data_list$D), function(i) 
                  merge_hours(data_list$interactions[i,])))

#Rounded to 5 minute timestamps
active_m5 <- active_consec
active_m5[,1] <- floor(active_consec[,1]*12)
active_m5[,2] <- ceiling(active_consec[,2]*12)
ints_m5 <- floor(ints_consec*12)

possible.times <- c(seq(break_times[1],break_times[2]-1,60*5),
                    seq(break_times[3],break_times[4]-1,60*5),
                    seq(break_times[5],break_times[6]-1,60*5),
                    seq(break_times[7],break_times[8]-1,60*5),
                    seq(break_times[9],break_times[10]-1,60*5),
                    seq(break_times[11],break_times[12]-1,60*5),
                    seq(break_times[13],break_times[14]-1,60*5),
                    seq(break_times[15],break_times[16],60*5))
num.time.steps <- length(possible.times)
possible_m5 <- 0:(num.time.steps-1)

look.back.default <- 7
alpha.vec.default <- seq(0, 1, length=look.back.default + 1)
alpha.vec.length <- length(alpha.vec.default)

#Modifed function to estimate network at times (taking into account rounding)
ship_snap <- function(i, t) {
   time_snap = possible.times[t]
   if(length(times[[i]]) == 0) {
      0
   } else if (max(active_m5[tag_ids %in% indices[i,],1]) > possible_m5[t] | 
              min(active_m5[tag_ids %in% indices[i,],2]) < possible_m5[t]) {
      0 #NA
   } else if (max(times[[i]]) <= time_snap) {
      sample_ships[i,data_list$N[i]]
   } else if (time_snap <= min(times[[i]])) {
      sample_ships[i,1]
   } else {
      index <- sum(difftime(times[[i]],time_snap) < 0)
      w1 <- as.numeric(difftime(time_snap,times[[i]][index], unit = "hours"))
      w2 <- as.numeric(difftime(times[[i]][index+1],time_snap, unit = "hours"))
      w2/(w1+w2)*sample_ships[i,index] + w1/(w1+w2)*sample_ships[i,index+1]
   }
}
snapshots <- sapply(1:num.time.steps, function(i) 
                    sapply(1:(data_list$D), ship_snap, t=i))

#Fade connections as nodes de-activate
for(i in 1:(data_list$D)) {
   last_index = min(active_m5[tag_ids %in% indices[i,],2])
   if(last_index <= max(possible_m5)) {
      fade_len = 1:(min(last_index+look.back.default,
                        max(possible_m5))-last_index)
      snapshots[i,fade_len+last_index] = (1-fade_len/look.back.default)*
         snapshots[i,last_index]
   }
}  

m <- matrix(data = 0, nrow = length(tag_ids), ncol = length(tag_ids))
net_inferred <- graph.adjacency(m,mode = "undirected", weighted = TRUE)
V(net_inferred)$color <- nodes.cols
V(net_inferred)$label.color <- "black"
V(net_inferred)$red = gsub("M",0,gsub("F",1,covariates$Sex))
V(net_inferred)$blue = gsub("M",1,gsub("F",0,covariates$Sex))
V(net_inferred)$name = tag_ids
net_raw <- net_inferred

net_inferred <- add.edges(net_inferred, 
                 as.character(sapply(rep(1:(data_list$D),num.time.steps), 
                                function(i) indices[i,1:2])))
E(net_inferred)$m5 = sort(rep(possible_m5+1,data_list$D))
E(net_inferred)$alpha = as.numeric(snapshots)
E(net_inferred)$color = "black"

net_raw <- add.edges(net_raw,
                     as.character(sapply(rep(1:(data_list$D),data_list$N), 
                                         function(i) indices[i,1:2])))
E(net_raw)$m5 = unlist(sapply(1:(data_list$D), function(i) 
   if(data_list$N[i]) ints_m5[i,1:data_list$N[i]]))
E(net_raw)$color = "black"
E(net_raw)$green = unlist(sapply(1:(data_list$D), function(i) if(data_list$N[i]) 
   snapshots[i,ints_m5[i,1:(data_list$N[i])]+1] > 0.5))*1
E(net_raw)$alpha = 0


for(i in 1:num.time.steps) {
   look.back <- i - look.back.default
   if (look.back < 0) {
      look.back <- 1
   }
   date.fade.index <- (look.back-1):(i-1)
   date.fade.index.length <- length(date.fade.index)
   alpha.vec <- alpha.vec.default
   if ((alpha.vec.length - date.fade.index.length) > 0) {
      alpha.vec <- alpha.vec[-(1:(alpha.vec.length - date.fade.index.length))]
   }
   for (j in 1:date.fade.index.length) {
      active.edges <- which(E(net_raw)$m5 %in% date.fade.index[j])
      if (length(active.edges) > 0) {
         E(net_raw)[active.edges]$alpha <- alpha.vec[j]
      }
   }  
   
   V(net_raw)$color <- rgb(0,0,0,0.4)
   V(net_raw)$label.color <- rgb(0,0,0)
   isActive = ((i-1)>=active_m5[,1] & (i-1)<=active_m5[,2])
   V(net_raw)$color[isActive] <- nodes.cols[isActive]
   V(net_raw)$label.color[isActive] <- nodes.cols.bold[isActive]
   nowInactive = (i-1)-active_m5[,2] <= 
      look.back.default & (i-1)-active_m5[,2] > 0
   smooth = (1-(i-1-active_m5[nowInactive,2])/look.back.default)
   if(length(smooth)) {
      V(net_raw)$color[nowInactive] <- 
         rgb(smooth*as.numeric(V(net_raw)$red[nowInactive]), 0, 
             smooth*as.numeric(V(net_raw)$blue[nowInactive]), 
             0.4 - 0.2*as.numeric(V(net_raw)$red[nowInactive]))
      V(net_raw)$label.color[nowInactive] <- 
         rgb(smooth*as.numeric(V(net_raw)$red[nowInactive]), 0, 
             smooth*as.numeric(V(net_raw)$blue[nowInactive]))
   }
   V(net_inferred)$color = V(net_raw)$color
   V(net_inferred)$label.color = V(net_raw)$label.color
   E(net_raw)$color = rgb(0,E(net_raw)$green,0,E(net_raw)$alpha)
   E(net_inferred)[E(net_inferred)$m5 != i]$color = rgb(0,0,0,0)
   E(net_inferred)[E(net_inferred)$m5 == i]$color = 
      rgb(0,0,0,E(net_inferred)[E(net_inferred)$m5 == i]$alpha)
   
   png(paste("figures/animated_net/NetAnimation",
             sprintf("%03d", i),".png",sep=""),
       7.5,3,units="in",pointsize=6,res = 800)
   par(mfrow = c(1,2), mar = c(1,3,2,3), oma = c(1,0,0,0))
   plot.igraph(net_raw,layout=fr_layout, edge.curved = 0, edge.width = 2, 
               main = "Encounters")
   mtext(strftime(possible.times[i]), side = 1, adj = 1.4, cex = 1.2)
   plot.igraph(net_inferred,layout=fr_layout, edge.curved = 0, edge.width = 2, 
               main = "Inferred Network")
   dev.off()
}

#Transform sequence of .png into .avi
#Requires ffmpeg, which can be found here: https://ffmpeg.org/
#ffmpeg -framerate 10 -i NetAnimation%03d.png -c:v libx264 -r 60 -pix_fmt yuv420p barn_swallow_network.mp4
