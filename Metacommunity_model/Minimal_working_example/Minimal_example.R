## Minimal_working_example.R

# provides minimal working examples on how to run the metacommunity model

rm(list = ls())

## sources ##################################################################################################################
source("./DesignFunctions.R")
source("./Parameters.R")
source("./DynamicsFunctions.R")


# Metacommunity without extinctions with interaction parameters fitted from empirical data ("emp." scenario) #################

  # generating a random initial state:
  initial.state1 = generate.initial.state()

  # loading parameters 
  parameters1 = generate.parameters(scenario = "emp.", migration_rate = 1/100)
  
  # specifying an extinction vector (in this case an empty one)
  extinction.vector1 = c()
  
  example1 = dynamics(initial.state = initial.state1,
                      parameters = parameters1,
                      extinctions = extinction.vector1,
                      Thalf = 500,  #time at which the extinctions happen (if extinction.vector1 was not empty)
                      Tmax = 1000,  #endpoint of the simulation
                      by = 20)      #time step for output
  
  #plotting the dynamics in patch 1
  plot(example1[,1], example1[, 32+2], type = "l", col = "green")
  points(example1[,1], example1[, 16+2], type = "l", col = "red")
  points(example1[,1], example1[, 2], type = "l", col = "blue")
  
  #plotting the dynamics in patch 9
  plot(example1[,1], example1[, 32+10], type = "l", col = "green")
  points(example1[,1], example1[, 16+10], type = "l", col = "red")
  points(example1[,1], example1[, 10], type = "l", col = "blue")
  

# Metacommunity with 8 clumped extinctions ("comp.-col." scenario) #################
  
  
  # generating an initial state from an existing species map
  initial.state2 = generate.initial.state(c("Tet", "Ble", "Col", "Col",
                                            "Col", "Tet", "Ble", "Ble",
                                            "Tet", "Col", "Tet", "Ble",
                                            "Ble", "Tet", "Tet", "Col"))
  
  # loading parameters 
  parameters2 = generate.parameters(scenario = "comp.-col.", migration_rate = 1/100)
  
  # specifying an extinction vector (half extinct, half non-extinct, clumped in space)
  extinction.vector2 = c(1:8)
  
  example2 = dynamics(initial.state = initial.state2,
                      parameters = parameters2,
                      extinctions = extinction.vector2,
                      Thalf = 500,  #time at which the extinctions happen (if extinction.vector1 was not empty)
                      Tmax = 1000,  #endpoint of the simulation
                      by = 5)      #time step for output
  
  #plotting the dynamics in patch 1
  plot(example2[,1], example2[, 32+2], type = "l", col = "green")
  points(example2[,1], example2[, 16+2], type = "l", col = "red")
  points(example2[,1], example2[, 2], type = "l", col = "blue")
  
  #plotting the dynamics in patch 9
  plot(x = example2[,1], y = example2[, 32+10], type = "l", col = "green")
  points(example2[,1], example2[, 16+10], type = "l", col = "red")
  points(example2[,1], example2[, 10], type = "l", col = "blue")
  
  