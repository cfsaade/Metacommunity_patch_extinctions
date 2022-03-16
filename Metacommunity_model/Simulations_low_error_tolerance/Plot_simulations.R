rm(list = ls())

## Libraries and themes ########################################################
library(ggplot2)
library(ggpubr)
library(data.table)
theme_set(theme_classic(base_size = 14))

## Reading data ################################################################
community_list = c("fit", "trade_off", "no_interaction", "randomized")

full_data = fread("./out/simulations/full_data.csv", header = T)
full_data$correlation_label[full_data$n_extinctions == 0] = "No extinction"

landscape_table = read.csv("./misc/SimulationsLandscapeTable.csv")

#changing labels of scenario
full_data$community[full_data$community == "fit"] = "emp."
full_data$community[full_data$community == "randomized"] = "rand."
full_data$community[full_data$community == "no_interaction"] = "no_int."
full_data$community[full_data$community == "trade_off"] = "comp.-col."

community_list = c("emp.", "comp.-col.", "no_int.", "rand.")

# changing label of correlation
full_data$correlation_label[full_data$correlation_label == "Clumped extinctions"] = "Clustered extinctions"

full_data$correlation_label = as.factor(full_data$correlation_label)

full_tube_table = full_data[full_data$time == 0,]
full_tube_table$recovery_time = 0

#computing return times

steady_state_bioarea = full_data$bioarea[full_data$time < 300]

## setting the threshold for a return time as the 2.5% percentile
bioarea_threshold = quantile(steady_state_bioarea, 0.025)
bla = max(full_tube_table$id)

for (community in community_list){
  print(community)
  for (tube_id in unique(full_data$id[full_data$extinctions == T])){
    print(tube_id/bla*100)
    test_vector = full_data$bioarea[full_data$id == tube_id & full_data$time > 300 & full_data$community == community]
    temp = full_data$time[full_data$id == tube_id & full_data$time > 300 & full_data$community == community]
    return_time = temp[min(which(test_vector > bioarea_threshold))] - 300
    full_tube_table$return_time[full_tube_table$id == tube_id & full_tube_table$community == community] = return_time
  }
}

                                        
## Direct effects plots ########################################################
temp_data = full_data[(full_data$n_extinctions == 0 | full_data$extinctions == T) & full_data$time>350 & full_data$time < 450,]

  # color scale
colors = c("red", "blue", "green")
names(colors) = levels(full_data$correlation_label)
fillScale = scale_fill_manual(name = "Treatment",values = colors)
alpha_box = 0.3
  
  # Plot base
  gg_base_direct = ggplot(data = temp_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label)) +
    facet_grid(.~community) +
    fillScale
  
  # Simpson plot
  panel_simpson = gg_base_direct +
    geom_boxplot(aes(y = simpson), alpha = alpha_box, outlier.shape = NA) +
    labs(y = expression(paste(alpha, "-diversity (Simpson's index)")), x = "Extinctions")
  
  # beta div plot
  panel_beta = gg_base_direct +
    geom_boxplot(aes(y = BetaDiv), alpha = alpha_box, outlier.shape = NA) +
    labs(y = expression(paste(beta, "-diversity")), x = "Extinctions")
  
  ## bioarea plot
  panel_bioarea = gg_base_direct +
    geom_boxplot(aes(y = log(bioarea)), alpha = alpha_box, outlier.shape = NA) +
    labs(y = "log(Bioarea)", x = "Extinctions")
  
  ## return time:
  panel_recovery = ggplot(data = full_tube_table[full_tube_table$extinctions == T,], mapping = aes(x = as.factor(n_extinctions), fill = correlation_label)) +
    facet_grid(.~community) +
    fillScale +
    geom_boxplot(aes(y = return_time), alpha = alpha_box, outlier.shape = NA) +
    labs(y = "Recovery time", x = "Extinctions") +
    ylim(0, 190)
  
  
cairo_pdf(file = "./out/figures/direct_effects.pdf", 11.7, 8.3)
ggarrange(panel_bioarea, panel_recovery, panel_simpson, panel_beta,
            labels = "auto",
            common.legend = TRUE, legend = "bottom")
dev.off()

## indirect effects

temp_data = full_data[full_data$extinctions == F & full_data$time>350,]

# Plot base
gg_base_direct = ggplot(data = temp_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label)) +
  facet_grid(.~community) +
  fillScale

# Simpson plot
panel_simpson = gg_base_direct +
  geom_boxplot(aes(y = simpson), alpha = alpha_box, outlier.shape = NA) +
  labs(y = expression(paste(alpha, "-diversity (Simpson's index)")), x = "Extinctions")



## bioarea plot
panel_bioarea = gg_base_direct +
  geom_boxplot(aes(y = log(bioarea)), alpha = alpha_box, outlier.shape = NA) +
  labs(y = "log(Bioarea)", x = "Exinctions")


cairo_pdf(file = "./out/figures/indirect_effects.pdf", 11.7*1.1/2, 8.3*1.1)
ggarrange(panel_bioarea, panel_simpson,
          labels = "auto",
          common.legend = TRUE, legend = "bottom", nrow = 2)
dev.off()

## moving windows of alpha
tmin_list = seq(300, 570, by = 10)
tmax_list = seq(330, 600, by = 10)
N = length(tmin_list)

plot = list()
for (k in 1:N){
  tmin = tmin_list[k]
  tmax = tmax_list[k]
  
  temp_data = full_data[(full_data$n_extinctions == 0 | full_data$extinctions == T) & full_data$time > tmin & full_data$time <= tmax,]
  
  plot[[k]] =  ggplot(data = temp_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label)) +
    facet_grid(.~community) +
    fillScale +
    geom_boxplot(aes(y = simpson), alpha = alpha_box, outlier.shape = NA) +
    labs(y = "", x = "", title = paste("]", tmin, ",", tmax, "]", sep = ""))
}
  
cairo_pdf(file = "./out/figures/moving_window_alpha1.pdf", 11.7*5, 8.3*5)
ggarrange(plotlist = plot, common.legend = T)
dev.off()

##  
tmin_list = seq(300, 570, by = 30)
tmax_list = seq(330, 600, by = 30)
N = length(tmin_list)


plot = list()
for (k in 1:N){
  tmin = tmin_list[k]
  tmax = tmax_list[k]
  
  temp_data = full_data[(full_data$n_extinctions == 0 | full_data$extinctions == T) & full_data$time > tmin & full_data$time <= tmax,]
  
  plot[[k]] =  ggplot(data = temp_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label)) +
    facet_grid(.~community) +
    fillScale +
    geom_boxplot(aes(y = simpson), alpha = alpha_box, outlier.shape = NA) +
    labs(y = "", x = "", title = paste("]", tmin, ",", tmax, "]", sep = "")) +
    theme(strip.text.x = element_text(size = 8, angle = 90))
}

cairo_pdf(file = "./out/figures/moving_window_alpha2.pdf", 11.7, 8.3)
ggarrange(plotlist = plot, common.legend = T, nrow = 2, ncol = 5)
dev.off()



  



