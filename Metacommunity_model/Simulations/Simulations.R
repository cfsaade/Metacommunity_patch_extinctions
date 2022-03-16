## Simulations.R
##
## Execute the code below to reproduce the simulations (and reproduce the figures using Plot_simulations.R)
##
## !!!!!
## This code is written to run in parallel on 10 cores and takes a while to run. Please ensure that your computer has
## enough core or change the parameter mc.core (l. 209)
## !!!!



rm(list = ls())
source("./misc/DesignFunctions.R")
source("./misc/DynamicsFunctions.R")
source("./misc/Parameters.R")

library(ggplot2)
library(ggstance)
library(lme4)
library(qpcR)
library(MASS)
library(MuMIn)
library(xtable)
library(adespatial)
library(parallel)
options(na.action = na.fail)


dir.create("../out/simulations", recursive = T)

n_spec = 3

simpson.index = function(community){
  if (max(community) == 0){
    return(1)
  }
  N = sum(community)
  D = 0
  for (n in community){
    D = D + (n/N)^2
  }
  return(1/D)
}

get.initial.state = generate.initial.state

bioarea_function = function(community){
  ## a function for calculating total bioarea given a community (ble, col, tet)
  ## values are the average bioarea of each species relative to that of Tet
  return(9.62*community[1] + 2.08*community[2] + community[3])
}


tube_table <- read.csv("./misc/SimulationsTubeTable.csv")
landscape_table <- read.csv("./misc/SimulationsLandscapeTable.csv")

simul = function(scenario = "fit"){
  ## dynamics simulation
  parameters = generate.parameters(n_spec, scenario)
  
  landscape = 1
  density_data = tube_table[tube_table$landscape == landscape,]
  density_data$correlation = landscape_table$cor[1]
  density_data$n_extinctions = landscape_table$extinctions[1]
  
  initial_species_vect = tube_table$initial.species[tube_table$landscape == landscape]
  initial.state = get.initial.state(initial_species_vect)
  extinction.map = which(tube_table$extinctions[tube_table$landscape == landscape])
  
  out = dynamics(initial.state, parameters, extinction.map, Thalf = 300, Tmax = 600, by = 10)
  #write.csv(out, file=paste("out_simulation_", landscape, ".csv", sep = ""))
  
  density_data$time = out[1, 1]
  density_matrix = matrix(out[1, -1], nrow = 16)
  colnames(density_matrix) = letters[1:n_spec]
  density_data = cbind(density_data, density_matrix)
  
  for (t in 2:nrow(out)){
    temp_density_data = tube_table[tube_table$landscape == landscape,]
    temp_density_data$correlation = landscape_table$cor[1]
    temp_density_data$n_extinctions = landscape_table$extinctions[1]
    
    temp_density_data$time = out[t, 1]
    density_matrix = matrix(out[t, -1], nrow = 16)
    colnames(density_matrix) = letters[1:n_spec]
    
    temp_density_data = cbind(temp_density_data, density_matrix)
    
    density_data = rbind(density_data, temp_density_data)
  }
  for (landscape in 2:max(landscape_table$id)){
    print(paste(scenario, landscape))
    if (landscape %% 5 == 1){
      parameters = generate.parameters(n_spec, scenario)
    }
    new_density_data = tube_table[tube_table$landscape == landscape,]
    new_density_data$correlation = landscape_table$cor[landscape]
    new_density_data$n_extinctions = landscape_table$extinctions[landscape]
    
    initial_species_vect = tube_table$initial.species[tube_table$landscape == landscape]
    initial.state = get.initial.state(initial_species_vect)
    extinction.map = which(tube_table$extinctions[tube_table$landscape == landscape])
    
    out = dynamics(initial.state, parameters, extinction.map, Thalf = 300, Tmax = 600, by = 10)
    
    new_density_data$time = out[1, 1]
    density_matrix = matrix(out[1, -1], nrow = 16)
    colnames(density_matrix) = letters[1:n_spec]
    new_density_data = cbind(new_density_data, density_matrix)
    
    for (t in 2:nrow(out)){
      temp_density_data = tube_table[tube_table$landscape == landscape,]
      temp_density_data$correlation = landscape_table$cor[landscape]
      temp_density_data$n_extinctions = landscape_table$extinctions[landscape]
      
      temp_density_data$time = out[t, 1]
      density_matrix = matrix(out[t, -1], nrow = 16)
      colnames(density_matrix) = letters[1:n_spec]
      
      temp_density_data = cbind(temp_density_data, density_matrix)
      
      new_density_data = rbind(new_density_data, temp_density_data)
    }
    density_data = rbind(density_data, new_density_data)
  }
  
  ## analysis
  
  density_data$community = scenario
  
  n = nrow(density_data)
  
  
  
  ## Table with informations on the experimental design (treatments, extinction positions...)
  tube_table = read.csv("./misc/SimulationsTubeTable.csv")
  landscape_table = read.csv("./misc/SimulationsLandscapeTable.csv")
  
  density_data$correlation_label = rep("", n)
  
  density_data$correlation_label[!density_data$correlation] = "Dispersed extinctions"
  density_data$correlation_label[density_data$correlation] = "Clumped extinctions"
  density_data$n_extinctions_label = paste(density_data$n_extinctions, "extinctions", sep = " ")
  density_data$n_extinctions_factor = as.factor(density_data$n_extinctions)
  
  density_data$extinction_label = rep("", n)
  density_data$extinction_label[density_data$extinctions] = "Extinct patches"
  density_data$extinction_label[!density_data$extinctions] = "Non-extinct patches"
  
  
  tube_table$n_extinctions = landscape_table$extinctions[tube_table$landscape]
  tube_table$correlation = landscape_table$cor[tube_table$landscape]
  
  
  tube_table$correlation_label[!tube_table$correlation] = "Dispersed extinctions"
  tube_table$correlation_label[tube_table$correlation] = "Clumped extinctions"
  tube_table$n_extinctions_label = paste(tube_table$n_extinctions, "extinctions", sep = " ")
  
  tube_table$extinction_label = rep("", nrow(tube_table))
  tube_table$extinction_label[tube_table$extinctions] = "Extinct patches"
  tube_table$extinction_label[!tube_table$extinctions] = "Non-extinct patches"
  
  
  density_data$BetaDiv = 0  ## beta diversity
  density_data$ReplDiv = 0  ## beta diversity due to species replacement
  density_data$RichDiv = 0  ## beta diversity due to species richness differences
  
  density_data$landscape_richness = rep(0, n)
  
  
  for (t in unique(density_data$time)){
    for (k in unique(density_data$landscape)){
      test = density_data$time == t & density_data$landscape == k
      density_data[test, letters[1:n_spec]] = density_data[test, letters[1:n_spec]] * (density_data[test, letters[1:n_spec]] > 0)
      density_data$landscape_richness[test] = sum(apply(density_data[test, letters[1:n_spec]], 2, sum) > 0)
    }
  }
  
  
  
  for (landscape in unique(density_data$landscape)){
    if (landscape%%10 == 1){
      print(paste(scenario, "analysis", landscape))
    }
    for (time in unique(density_data$time)){
      beta = beta.div.comp(density_data[density_data$landscape == landscape & density_data$time == time,
                                        letters[1:3]], quant = T)$part
      density_data$BetaDiv[density_data$landscape == landscape & density_data$time == time] = beta["BDtotal"]
      density_data$ReplDiv[density_data$landscape == landscape & density_data$time == time] = beta["Repl"]
      density_data$RichDiv[density_data$landscape == landscape & density_data$time == time] = beta["RichDif"]
    }
  }
  
  density_data$bioarea = 0
  density_data$simpson = 0
  for (k in 1:n){
    if (k%%1000 == 1){
      print(paste(scenario, floor(k/n*100), '%'))
    }
    density_data$bioarea[k] = bioarea_function(unlist(density_data[k, letters[1:n_spec]]))
    density_data$simpson[k] = simpson.index(density_data[k, letters[1:n_spec]])
  }
  
  write.csv(density_data, file= paste("../out/simulations/density_data_", scenario, ".csv", sep = ""))
}

scenarios_list = c("fit", "trade_off", "randomized", "no_interaction")

mclapply(scenarios_list, simul, mc.cores = 10, mc.preschedule = F, mc.silent = F)  ## applying simul to all scenarios

full_data = read.csv("../out/simulations/density_data_fit.csv")
for (scenario in scenarios_list[2:4]){
  temp = read.csv(paste("../out/simulations/density_data_", scenario, ".csv", sep = ""), header = T)
  full_data = rbind(full_data, temp)
}

write.csv(full_data, "../out/simulations/full_data.csv")





















