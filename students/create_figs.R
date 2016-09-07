# 1. Create assorted figures and animations for MIT Social Evolution data




#Load in requisite libraries and processed daa
##############################################

#Requires package "Cairo"
setHook(packageEvent("grDevices", "onLoad"),
        function(...) grDevices::X11.options(type='cairo'))
options(device='x11')

library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(rstan)
library(igraph)
library(Matrix)

load("processed_data/stan_pre.RData")
load("processed_data/stan_output.RData")

pre_load = stanc(file = "students.stan")
barn_model = stan_model(stanc_ret = pre_load)
expose_stan_functions(barn_model)

as.time = function(t) as.POSIXct(t, origin = "1970-01-01")





#Descriptive plots
##################

#Multi-plot function adapted from https://stackoverflow.com/questions/30611474/
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  p <- grid.arrange(
    do.call(arrangeGrob, c(lapply(2:length(plots), function(x)
      plots[[x]] + theme(legend.position="none")),list(nrow = 1))),
    legend,
    nrow = 2,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
  return(p)
}

#Basic descriptives of proximity interaction data
during_night <- !(weekly %% 48 >= 8*2 & weekly %% 48 < 17*2)
months <- as.POSIXlt(unlist(proximity2), origin = "1970-01-01")$mon

p1 <- ggplot(NULL, aes(as.time(unlist(proximity2)), fill = during_night)) +
   geom_histogram(bins = 60) +
   scale_fill_discrete(labels = c(TRUE,FALSE),
                       guide = guide_legend(title = "During Day?")) +
   ggtitle("Interactions Over Time") +
   ylab("Number of Interactions") +
   xlab("Date") +
   scale_x_datetime(date_minor_breaks = "1 month", 
      breaks = as.POSIXct(c("2008-09-01","2008-10-01","2008-11-01","2008-12-01",
         "2009-01-01","2009-02-01", "2009-03-01","2009-04-01","2009-05-01",
         "2009-06-01"), tz = "GMT"), 
      labels = c("","Oct","","","Jan","","","Apr","","")) +
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12))
p2 <- ggplot(NULL, aes(weekly[months >= 9], fill = during_night[months >= 9])) +
   geom_histogram(bins = 336/2) +
   ggtitle("Weekly Periodicity (Oct to Dec)") +
   ylab("Number of Interactions") +
   scale_x_continuous(breaks = seq(0,336,24*2), name = "Time",
      labels = c("Su","M","T","W","Th","F","Sa","Su")) +
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12))
p3 <- ggplot(NULL, aes(weekly[months %in% c(3,4)], 
                       fill = during_night[months %in% c(3,4)])) +
   geom_histogram(bins = 336/2) +
   ggtitle("Weekly Periodicity (Apr to May)") +
   ylab("Number of Interactions") +
   scale_x_continuous(breaks = seq(0,336,24*2), name = "Time",
                      labels = c("Su","M","T","W","Th","F","Sa","Su")) +
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12))

cairo_ps("figures/fig3_descriptive_mit.eps", 
         width=9, height=3, fallback_resolution = 800)
p <- grid_arrange_shared_legend(p1, p1, p2, p3)
dev.off()






#Compare weekly periodicity patterns across friend and floor combinations
non_friends <- sapply(reported_net,sum) <= 1
friends <- sapply(reported_net,sum) >= 4
active_len <- data_list$active_times[,2] - data_list$active_times[,1]
prox_to_weekly <- function(subset,normalization) {
   ptimes = as.POSIXlt(unlist(proximity2[subset]), origin = "1970-01-01")
   weekly = ptimes$wday*48+ptimes$hour*2+(ptimes$min>=30)
   as.numeric(table(weekly))/sum(active_len[subset])/normalization
}

wv <- predict(smooth.spline(
               prox_to_weekly(non_friends & !data_list$same_floor,1e-5), 
               spar = 0.2))$y
wv2 <- predict(smooth.spline(
               prox_to_weekly(friends & !data_list$same_floor,1e-5), 
               spar = 0.2))$y
wv3 <- predict(smooth.spline(
               prox_to_weekly(non_friends & data_list$same_floor,1e-5), 
               spar = 0.2))$y
wv4 <- predict(smooth.spline(
               prox_to_weekly(friends & data_list$same_floor,1e-5), 
               spar = 0.2))$y

labels = c("Non-friends / Different Floors","Friends / Different Floors",
           "Non-friends / Same Floor","Friends / Same Floor")
weekly_patterns <- data.frame(x = rep(1:336,4), ints = c(wv,wv2,wv3,wv4), 
   factor = factor(sapply(labels, function(i) rep(i,336)), levels = labels))



factor.cols <- brewer.pal(4, "Paired")
#Rectangles to highlight days
rect <- data.frame(xmin=8*2+seq(0,335,24*2), xmax=17*2+seq(0,335,24*2),
                   ymin=rep(-Inf,7), ymax=rep(Inf,7))

p0 <- ggplot(weekly_patterns, aes(x,ints, colour = factor)) +
   geom_line(aes(size = factor)) +
   scale_size_manual(values = c(1, 1, 1.5, 1.5)) +
   scale_colour_manual(values = factor.cols[c(2,1,4,3)]) +
   guides(size=guide_legend(title = "Legend",keywidth = 2, keyheight = 1),
          colour=guide_legend(title = "Legend")) +
   theme_bw() +
   theme(plot.title = element_text(size=12), text = element_text(size=12))
p1 <- ggplot(weekly_patterns[c(1:336,336+1:336),], 
             aes(x,ints,group = factor,colour = factor)) +
   geom_line(aes(size = factor), alpha=0.8) +
   ylab("Normalized Intensity") +
   scale_x_continuous(breaks = seq(0,336,24*2), 
      labels = c("Su","M","T","W","Th","F","Sa","Su"), name = "") +
   scale_size_manual(values = c(1, 1)) +
   scale_colour_manual(values = factor.cols[c(2,1)]) +
   ylim(0,17.5) +
   theme_bw() +
   geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
               color="grey", linetype = 0, alpha=0.2, inherit.aes = FALSE) +
   theme(plot.title = element_text(size=12), text = element_text(size=12))
p2 <- ggplot(weekly_patterns[c(336*2+1:336,336*3+1:336),], 
             aes(x,ints,group = factor, colour = factor)) +
   geom_line(aes(size = factor), alpha=0.8) +
   ylab("Normalized Intensity") +
   scale_size_manual(values = c(1.5, 1.5)) +
   scale_x_continuous(breaks = seq(0,336,24*2), 
      labels = c("Su","M","T","W","Th","F","Sa","Su"), name = "") +
   scale_colour_manual(values = factor.cols[c(4,3)]) +
   ylim(0,17.5) +
   theme_bw() +
   geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
             color="grey", linetype = 0, alpha=0.2, inherit.aes = FALSE) +
   theme(plot.title = element_text(size=12), text = element_text(size=12))
p3 <- ggplot(weekly_patterns[c(1:336,336*2+1:336),], 
             aes(x,ints,group = factor, colour = factor)) +
   geom_line(aes(size = factor), alpha=0.8) +
   ylab("Normalized Intensity") +
   scale_x_continuous(breaks = seq(0,336,24*2), 
      labels = c("Su","M","T","W","Th","F","Sa","Su"), name = "") +
   scale_size_manual(values = c(1, 1.5))   +
   scale_colour_manual(values = factor.cols[c(2,4)]) +
   ylim(0,17.5) +
   theme_bw() +
   geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
             color="grey", linetype = 0, alpha=0.2, inherit.aes = FALSE) +
   theme(plot.title = element_text(size=12), text = element_text(size=12))

cairo_ps("figures/fig4_mit_periodicity.eps", 
         width=9, height=3, fallback_resolution = 800)
p <- grid_arrange_shared_legend(p0, p1, p2, p3)
dev.off()






#Estimating CTMC parameters from survey data
#Estimate proportion of connected dyads
estd_sparsity <- rbind(rowMeans(sapply(reported_net, identity)
                                [,data_list$same_floor & data_list$same_year]),
                  rowMeans(sapply(reported_net, identity)
                           [,data_list$same_floor & !data_list$same_year]),
                  rowMeans(sapply(reported_net, identity)
                           [,!data_list$same_floor & data_list$same_year]),
                  rowMeans(sapply(reported_net, identity)
                           [,!data_list$same_floor & !data_list$same_year]))
s_hat <- rowMeans(estd_sparsity)

prop_changed <- function(x) {
  x_len = length(x)
  abs(x[2:x_len] - x[2:x_len-1])
}

#Estimate probability of changes
prop_changes <- rbind(rowMeans(sapply(reported_net, prop_changed)
                               [,data_list$same_floor & data_list$same_year]),
                      rowMeans(sapply(reported_net, prop_changed)
                               [,data_list$same_floor & !data_list$same_year]),
                      rowMeans(sapply(reported_net, prop_changed)
                               [,!data_list$same_floor & data_list$same_year]),
                      rowMeans(sapply(reported_net, prop_changed)
                               [,!data_list$same_floor & !data_list$same_year]))

survey_difft <- diff(as.numeric(difftime(survey_dates[1:5],start,unit="hours")))

est_q <- function(subset,period) {
  temp = (s_hat[subset]*(1-estd_sparsity[subset,period]) 
          + (1-s_hat[subset])*estd_sparsity[subset,period])
  -log(1 - prop_changes[subset,period]/temp)/survey_difft[period]
}

estd_q <- sapply(1:4, function(p) sapply(1:4, function(s) est_q(s,p)))

#Convert to data.frame for ggplot
labels = c("Same Floor / Same Year", "Same Floor / Different Years",
           "Different Floors / Same Year", "Different Floors / Different Years")
sparsity_df <- data.frame(x = sort(rep(1:5,4)), level = rep(labels,5),
                          sparsity = as.numeric(estd_sparsity))
q_df <- data.frame(x = sort(rep(1:4,4)), level = rep(labels,4),
                   q = as.numeric(estd_q))

p1 <- ggplot(sparsity_df, aes(x=x, y=sparsity, group = level, color = level)) +
         geom_line() +
         geom_point() +
         ylim(0,1) +
         ylab(paste("Empirical",expression(s))) +
         scale_x_continuous(breaks = seq(1:5), 
         labels = c("9/9","10/19","12/13","3/5","4/17"), name = "Date") +
         guides(colour=guide_legend(title = "")) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

p2 <- ggplot(q_df, aes(x=x, y=q, group = level, color = level)) +
         geom_line() +
         geom_point() +
         ylim(0,1e-3) +
         ylab(paste("Empirical",expression(q))) +
         scale_x_continuous(breaks = seq(1:4), name = "Period of Transition",
         limits = c(1,4.25),
         labels = c("9/9 to 10/19","10/19 to 12/13",
                    "12/13 to 3/5","3/5 to 4/17")) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

cairo_ps("figures/fig5_survey_ctmc.eps", 
         width=9, height=3, fallback_resolution = 800)
p <- grid_arrange_shared_legend(p1,p1,p2)
dev.off()





#Comparing survey and interaction data
breaks <- c(0,10,50,100,200,400,600)
survey_means <- sapply(reported_net, mean)

p1 <- ggplot(NULL, aes(x = sqrt(data_list$N[survey_means == 0]))) + 
        geom_histogram(bins = 30, aes(y=..density..)) +
         scale_x_continuous(breaks = sqrt(breaks), labels = breaks, 
                            limits = c(NA,26)) +
         ylim(0,0.6) +
         ylab("Density") +
         xlab("Number of Interactions") +
         ggtitle("Reported Friendship 0/5") +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

p2 <- ggplot(NULL, aes(x = sqrt(data_list$N[survey_means == 0.2]))) + 
        geom_histogram(bins = 30, aes(y=..density..)) +
         scale_x_continuous(breaks = sqrt(breaks), labels = breaks, 
                            limits = c(NA,26)) +
         ylim(0,0.6) +
         ylab("Density") +
         xlab("Number of Interactions") +
         ggtitle("Reported Friendship 1/5") +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

p3 <- ggplot(NULL, aes(x = sqrt(data_list$N[survey_means == 0.4]))) + 
         geom_histogram(bins = 30, aes(y=..density..)) +
         scale_x_continuous(breaks = sqrt(breaks), labels = breaks, 
                            limits = c(NA,26)) +
         ylim(0,0.6) +
         ylab("Density") +
         xlab("Number of Interactions") +
         ggtitle("Reported Friendship 2/5") +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

p4 <- ggplot(NULL, aes(x = sqrt(data_list$N[survey_means == 0.6]))) + 
         geom_histogram(bins = 30, aes(y=..density..)) +
         scale_x_continuous(breaks = sqrt(breaks), labels = breaks, 
                            limits = c(NA,26)) +
         ylim(0,0.6) +
         ylab("Density") +
         xlab("Number of Interactions") +
         ggtitle("Reported Friendship 3/5") +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

p5 <- ggplot(NULL, aes(x = sqrt(data_list$N[survey_means == 0.8]))) + 
         geom_histogram(bins = 30, aes(y=..density..)) +
         scale_x_continuous(breaks = sqrt(breaks), labels = breaks, 
                            limits = c(NA,26)) +
         ylim(0,0.6) +
         ylab("Density") +
         xlab("Number of Interactions") +
         ggtitle("Reported Friendship 4/5") +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

p6 <- ggplot(NULL, aes(x = sqrt(data_list$N[sapply(reported_net, mean) == 1]))) + 
         geom_histogram(bins = 30, aes(y=..density..)) +
         scale_x_continuous(breaks = sqrt(breaks), labels = breaks, 
                            limits = c(NA,26)) +
         ylim(0,0.6) +
         ylab("Density") +
         xlab("Number of Interactions") +
         ggtitle("Reported Friendship 5/5") +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

cairo_ps("figures/fig12_ints_vs_surveys.eps", 
         width=9, height=4, fallback_resolution = 800)
grid.arrange(p1,p2,p3,p4,p5,p6,nrow = 2)
dev.off()






#Trace plots
############
output <- extract(fit)

p1 <- ggplot(NULL, aes(x=1:9000, y=output$c_0[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(c["0,t2"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p2 <- ggplot(NULL, aes(x=1:9000, y=output$c_0[,2])) + geom_line() +
         xlab("Iterations") + ylab(expression(c["0,t3"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p3 <- ggplot(NULL, aes(x=1:9000, y=output$c_1)) + geom_line() +
         xlab("Iterations") + ylab(expression(c[1])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p4 <- ggplot(NULL, aes(x=1:9000, y=output$c_2[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(c["2,1"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p5 <- ggplot(NULL, aes(x=1:9000, y=output$c_2[,2])) + geom_line() +
         xlab("Iterations") + ylab(expression(c["2,2"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p6 <- ggplot(NULL, aes(x=1:9000, y=output$c_3)) + geom_line() +
         xlab("Iterations") + ylab(expression(c[3])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p7 <- ggplot(NULL, aes(x=1:9000, y=output$k_0)) + geom_line() +
         xlab("Iterations") + ylab(expression(k[0])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p8 <- ggplot(NULL, aes(x=1:9000, y=output$k_1[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(k["1,1"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p9 <- ggplot(NULL, aes(x=1:9000, y=output$k_1[,2])) + geom_line() +
         xlab("Iterations") + ylab(expression(k["1,2"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p10 <- ggplot(NULL, aes(x=1:9000, y=output$k_1[,3])) + geom_line() +
         xlab("Iterations") + ylab(expression(k["1,3"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p11 <- ggplot(NULL, aes(x=1:9000, y=output$k_3[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(k["3,1"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p12 <- ggplot(NULL, aes(x=1:9000, y=output$k_3[,2])) + geom_line() +
         xlab("Iterations") + ylab(expression(k["3,2"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p13 <- ggplot(NULL, aes(x=1:9000, y=output$k_3[,3])) + geom_line() +
         xlab("Iterations") + ylab(expression(k["3,3"])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p14 <- ggplot(NULL, aes(x=1:9000, y=output$s_0)) + geom_line() +
         xlab("Iterations") + ylab(expression(s[0])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p15 <- ggplot(NULL, aes(x=1:9000, y=output$s[,1])) + geom_line() +
         xlab("Iterations") + ylab(expression(s[1])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p16 <- ggplot(NULL, aes(x=1:9000, y=output$s[,2])) + geom_line() +
         xlab("Iterations") + ylab(expression(s[2])) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))
p17 <- ggplot(NULL, aes(x=1:9000, y=output$t_prob)) + geom_line() +
         xlab("Iterations") + ylab(expression(q)) +
         theme_bw() +
         theme(plot.title = element_text(size=12), text = element_text(size=12))

cairo_ps("figures/fig13_mit_traces.eps", 
         width=9, height=12, fallback_resolution = 800)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,
             ncol = 3)
dev.off()






#Function to estimate networks at various times
###############################################
est_snapshots <- function(time_hrs) {
   estd_ships <- lapply(1:(data_list$D), function(i) rep(0,length(time_hrs)))
   for(d in 1:data_list$D) {
      if(data_list$N[d] > 1) {
         for(t in 1:length(time_hrs)) {
            #If outside active interval, NA
            if(time_hrs[t] < data_list$active_times[d,1] | 
               time_hrs[t] > data_list$active_times[d,2]) {
               estd_ships[[d]][t] = NA
               
            #If time is before/after first/last midpoint, take midpoint value
            }else if(time_hrs[t] <= ship_ints[[d]][1]) {
               estd_ships[[d]][t] = est_ships[[d]][1]
               
            }else if(time_hrs[t] >= ship_ints[[d]][data_list$N[d]]) {
               estd_ships[[d]][t] = est_ships[[d]][data_list$N[d]]
               
            #Otherwise interpolate
            } else {
               ob = which(time_hrs[t] <= ship_ints[[d]][2:data_list$N[d]] & 
                             time_hrs[t] > ship_ints[[d]][2:data_list$N[d]-1])
               weight = (ship_ints[[d]][ob+1] - time_hrs[t])/
                  (ship_ints[[d]][ob+1] - ship_ints[[d]][ob])
               estd_ships[[d]][t] = weight*est_ships[[d]][ob] + 
                  (1-weight)*est_ships[[d]][ob+1]
            }
         } 
      } else if (data_list$N[d] == 1) {
         for(t in 1:length(time_hrs)) {
            if(time_hrs[t] >= data_list$active_times[d,1] &
               time_hrs[t] <= data_list$active_times[d,2]) {
               estd_ships[[d]][t] = est_ships[[d]][1]
            }
         }
      } else {
         for(t in 1:length(time_hrs)) {
            if(time_hrs[t] >= data_list$active_times[d,1] & 
               time_hrs[t] <= data_list$active_times[d,2]) {
               estd_ships[[d]][t] = 0
            }
         }
      }
   }
   estd_ships
}





#Networks highlighting individual #2
####################################

#Circle layout for all snapshots
circ_order <- invPerm(order(floor[r_users], decreasing = TRUE))
m <- matrix(data = 0, nrow = length(r_users), ncol = length(r_users))
colnames(m) <- r_users
rownames(m) <- r_users
net <- graph.adjacency(m,mode = "undirected", weighted = TRUE, diag = FALSE)
circ_layout = layout.circle(net)[circ_order,]
circ_layout[1,] = circ_layout[1,] + 0.2*c(-1,1)

#Estimate network at survey dates
time_hrs = as.numeric(difftime(survey_dates[1:5],start, unit = "hours"))
estd_ships <- est_snapshots(time_hrs)

#Calculate periods when indivduals are active
active_hours <- difftime(as.POSIXct(time_intervals, origin = "1970-01-01"),
                         start, unit = "hours")
active_users <- sapply(time_hrs, function(x) 
   active_hours[,1] < x & active_hours[,2] > x)

floor_colors = brewer.pal(8, name = "Set1")
floors = levels(floor[r_users])



#Convert list of relations to an adjacency matrix
ships_to_adj <- function(ships, snap) {
   m <- matrix(data = 0, nrow = length(r_users), ncol = length(r_users))
   m[lower.tri(m)] = sapply(ships, function(x) x[snap])
   m = m + t(m)
   m
}



#Function to plot snapshots from adjacency matrices
plot_snapshot <- function(adj_mat, layout, active, title = "") {
   net = graph.adjacency(adj_mat, mode = "undirected", 
                         weighted = TRUE, diag = FALSE)
   V(net)$name = c("2",rep("",length(r_users)-1))
   V(net)$color = as.character(floor[r_users])
   for(i in 1:8) {
      V(net)$color = gsub(floors[i+1],floor_colors[i], V(net)$color)
   }
   V(net)$color[!active] = rgb(1,1,1)
   E(net)$weight[is.na(E(net)$weight)] = 0
   E(net)$color = rgb(0.7,0.7,0.7,E(net)$weight/2)
   highlighted <- attr(E(net),"vnames") == "2|"
   E(net)$color[highlighted] = rgb(0,0,0,E(net)$weight[highlighted])
   E(net)$lwd = 1
   E(net)$lwd[highlighted] = 2
   plot.igraph(net, layout = circ_layout, edge.width = E(net)$lwd,
               vertex.label="", vertex.size = c(15,rep(10,length(r_users)-1)),
               main = title)
}


cairo_ps("figures/fig6_mit_snapshots.eps", 
         width=9, height=2, fallback_resolution = 800)
par(mfrow = c(1,4), mar = c(0,0,1,0))
plot_snapshot(ships_to_adj(reported_net,3), 
              layout, rep(TRUE,57), "12/13/08 - Reported Network")
plot_snapshot(ships_to_adj(reported_net,4), 
              layout, rep(TRUE,57), "03/05/09 - Reported Network")
plot_snapshot(ships_to_adj(estd_ships,3), 
              layout, active_users[,3], "12/13/08 - Inferred Network")
plot_snapshot(ships_to_adj(estd_ships,4), 
              layout, active_users[,4], "03/05/09 - Inferred Network")
dev.off()





#Examine egocentric network around user #2
##########################################
#Restrict to users with relations to individual #2 (and individual 2)
users2 <- unique(c(r_users[!!ships_to_adj(estd_ships,3)[r_users == 2,]>0.5],
                   r_users[!!ships_to_adj(estd_ships,4)[r_users == 2,]>0.5]))
users2[is.na(users2)] = 2
users2 <- sort(users2)


days_of_interest = c("2008-12-13", "2008-12-24","2009-01-05","2009-01-17",
                     "2009-01-28","2009-02-03","2009-02-09","2009-03-05")
time_hrs = as.numeric(difftime(as.time(days_of_interest),start, unit = "hours"))
estd_ships <- est_snapshots(time_hrs)

active_users <- sapply(time_hrs, function(x) 
   active_hours[,1] < x & active_hours[,2] > x)

ego_order <- invPerm(order(floor[users2], decreasing = TRUE))
m <- matrix(data = 0, nrow = length(users2), ncol = length(users2))
net <- graph.adjacency(m,mode = "undirected", weighted = TRUE, diag = FALSE)
ego_layout = layout.circle(net)[ego_order,]
ego_layout[1,] = c(0,0)




plot_snapshot <- function(adj_mat, subset, layout, active, title = "") {
   net = graph_from_adjacency_matrix(adj_mat[subset,subset], 
                                     mode = "undirected", weighted = TRUE)
   V(net)$name = users2
   V(net)$label.cex = 1.25
   V(net)$label.color = "white"
   V(net)$color = as.character(floor[users2])
   for(i in 1:8) {
      V(net)$color = gsub(floors[i+1],floor_colors[i], V(net)$color)
   }
   V(net)$color[!active[subset]] = "white"
   V(net)$label.color[!active[subset]] = "black"
   E(net)$weight[is.na(E(net)$weight)] = 0
   E(net)$color = rgb(0.7,0.7,0.7,(E(net)$weight))
   plot.igraph(net, vertex.size = 30, layout = layout,
               edge.width=E(net)$weight*5, main = title)
}



cairo_ps("figures/fig7_mit_ego_snapshots.eps", 
         width=9, height=4, fallback_resolution = 800)

m <- matrix(c(1,2,3,4,5,6,7,8,9,9,9,9),nrow = 3,ncol = 4,byrow = TRUE)
layout(mat = m,heights = c(0.425,0.425,0.15))

for(i in 1:8) {
   par(mar = c(0,0,1,0))
   plot_snapshot(ships_to_adj(estd_ships,i),r_users %in% users2, ego_layout, 
              active_users[,i], days_of_interest[i])
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0, title = "Estimated Probability of a Relation",
       legend = c("0.25","0.5","0.75","1.0"), 
       col=rgb(0.7,0.7,0.7,(c(0.25,0.5,0.75,1))), lwd=5*c(0.25,0.5,0.75,1), 
       cex=1, horiz = TRUE, bty = "n")

dev.off()





#Movie of egocentric network
look.back.default <- 3
alpha.vec.default <- seq(0, 1, length=look.back.default + 1)
alpha.vec.length <- length(alpha.vec.default)

#Calculate inferred network at desired times
possible.times <- seq(as.POSIXct("2008-12-13"),as.POSIXct("2009-03-05"), 
                      by = 24*60*60)
num.time.steps <- length(possible.times)
indices2 <- which(indices[,1] %in% users2 & indices[,2] %in% users2)

time_hrs = difftime(possible.times,start, unit = "hours")
estd_ships <- est_snapshots(time_hrs)
estd_ships <- lapply(indices2, function(x) estd_ships[[x]])

time_to_stamp <- function(t) {
   floor((t - as.numeric(time_hrs[1]))/24)
}

active_stamps <- time_to_stamp(active_hours[r_users %in% users2,])



#Basic structure of network of raw data
net_d <- graph.adjacency(matrix(data=0,length(users2),length(users2)),
                         mode = "undirected", weighted = TRUE, diag = FALSE)
V(net_d)$name = users2
V(net_d)$label.cex = 1.25
V(net_d)$label.color = "white"
nodes.cols = as.character(floor[users2])
for(i in 1:8) {
   nodes.cols = gsub(floors[i+1],floor_colors[i], nodes.cols)
}
nodes.rgb <- t(col2rgb(nodes.cols))/255

#Basic structure and edges of inferred network
net_i <- net_d

net_i <- add.edges(net_i, as.character(
   sapply(rep(indices2,num.time.steps), function(i) indices[i,1:2])))
E(net_i)$timestamp = sort(rep(1:num.time.steps-1,choose(length(users2),2)))                                            
E(net_i)$alpha = as.numeric(t(sapply(estd_ships, function(x) x)))
E(net_i)$alpha[is.na(E(net_i)$alpha)] = 0
E(net_i)$color = grey(0.7)


#Add edges to data network (w/ color highlighting corresponding to the
#inferred network)
net_d <- add.edges(net_d, as.character(
   sapply(rep(indices2,data_list$N[indices2]), function(i) indices[i,1:2])))
E(net_d)$alpha = 0
E(net_d)$color = grey(0.7)
E(net_d)$timestamp = unlist(sapply(indices2, function(i) if(data_list$N[i]) 
   time_to_stamp(data_list$interactions[i,1:data_list$N[i]])))
E(net_d)$dyad = unlist(sapply(indices2, function(i) rep(i,data_list$N[i])))
E(net_d)$green = rep(0,length(E(net_d)))
for(i in 1:length(E(net_d))) {
   if(E(net_d)$timestamp[i] >= 0 & E(net_d)$timestamp[i] <= 
         as.numeric(time_to_stamp(time_hrs[length(time_hrs)]))) {
      E(net_d)$green[i] = (estd_ships[[which(indices2 %in% E(net_d)$dyad[i])]]
                              [E(net_d)$timestamp[i]] > 0.8)*1
   }
}


#Create a series of snapshots
for(i in 1:num.time.steps) {
   #Create a fade-out for data network
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
      active.edges <- which(E(net_d)$timestamp %in% date.fade.index[j])
      if (length(active.edges) > 0) {
         E(net_d)[active.edges]$alpha <- alpha.vec[j]
      }
   }  
   
   #Color indices as active/inactive
   V(net_i)$color <- "white"
   V(net_i)$label.color <- "white"
   isActive = ((i-1)>=active_stamps[,1] & (i-1)<=active_stamps[,2])
   V(net_i)$color[isActive] <- nodes.cols[isActive]
   nowInactive = (i-1)-active_stamps[,2] <= look.back.default & 
      (i-1)-active_stamps[,2] > 0
   smooth = (1-(i-1-active_stamps[nowInactive,2])/look.back.default)
   if(length(smooth)) {
      V(net_i)$color[nowInactive] <- rgb(smooth*nodes.rgb[nowInactive,1], 
                                         smooth*nodes.rgb[nowInactive,2], 
                                         smooth*nodes.rgb[nowInactive,3], 0.4)
   }
   V(net_d)$color = V(net_i)$color
   
   
   #Color edges (green in interaction network if inferred > 0.8)
   E(net_d)$color = rgb(0,E(net_d)$green,0,E(net_d)$alpha)
   E(net_i)[E(net_i)$timestamp != i-1]$color = rgb(0,0,0,0)
   E(net_i)[E(net_i)$timestamp == i-1]$color = 
      rgb(0,0,0,E(net_i)[E(net_i)$timestamp == i-1]$alpha)
   
   png(paste("figures/animated_net/NetAnimation",
             sprintf("%03d", i),".png",sep=""),
       7.5,3,units="in",pointsize=6,res = 800)
   par(mfrow = c(1,2), mar = c(1,3,2,3), oma = c(1,0,0,0))
   plot.igraph(net_d,layout=ego_layout, edge.curved = 0, edge.width = 5, 
               main = "Proximity Pings")
   mtext(strftime(possible.times[i]), side = 1, adj = 1.4, cex = 1.2)
   plot.igraph(net_i,layout=ego_layout, edge.curved = 0, edge.width = 5, 
               main = "Inferred Network")
   dev.off()
}

#Transform sequence of .png into .mp4
#Requires ffmpeg, which can be found here: https://ffmpeg.org/
#ffmpeg -framerate 5 -i NetAnimation%03d.png -c:v libx264 -r 60 -pix_fmt yuv420p egocentric_mit_network.mp4
