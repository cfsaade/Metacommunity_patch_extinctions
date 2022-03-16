
library(deSolve)


## derivative #################################################################################################
# Returns the derivative of a 3 species system (X, Y, Z) in size*size patches
#
# Parameters
#   t : time, unused but necessary for ode function
#   state : the poplulation densities in each of the 16 patches as a 48 value vector (1:16 = X, 17:32 = Y, 33:48 = Z)
#
#
# X, Y, and Z filter fed on an abiotic ressource (Beverton-Holt dynamics)
# Z eats Y and X (linear feeding rate)


derivative <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    
    # vector of derivatives (output)
    dState = matrix(rep(0, n_patches*n_spec), nrow = n_patches)
    
    # rewriting the species densities as a matrix    
    temp_state = matrix(state, nrow = n_patches)
    
    # computing local dynamics (competitive Lotka-volterra) and emmigration
    for (patch in 1:n_patches){
      dState[patch,] = temp_state[patch,] %*% diag(r) -
        (temp_state[patch,] %*% competition_matrix) %*% (diag(temp_state[patch,])) -
        (temp_state[patch,] %*% diag(migration_vector))
    }
    
    # computing immigration
    for(spec in 1:n_spec){
      dState[, spec] = dState[, spec] +
        neighbour_matrix %*% temp_state[, spec] * migration_vector[spec]
    }
    
    # returning the derivative
    list(dState)
  })
}

#### Dynamics ####################################################################################
## Returns the dynamics of a metacommunity undergoing extinction events after Thalf and recovering until Tmax
# value : out, a list of 3 elements, each the dynamics of 1 specie with time in row
# and patche in columns

dynamics <- function(initial.state, parameters, extinctions, Thalf = 500, Tmax = 1000, by = 20){
  n_patches = parameters$n_patches
  n_spec = parameters$n_spec
  times <- c(seq(0, Thalf, by = by), Thalf + 0.1)
  Thalf_index = length(times)
  ode.out <- ode(y = initial.state, times = times, func = derivative, parms = parameters, method = "lsoda", rtol = 10**-9, atol = 10**-9)
  for (spec_count in 0:(n_spec - 1)){
    ode.out[Thalf_index, extinctions + 1 + spec_count*n_patches] = 0
  }
  times = seq(Thalf + 0.1, Tmax + 0.1, by = by)
  ode.out = rbind(ode.out[1:Thalf_index,], ode(y = ode.out[Thalf_index,2:(n_spec*n_patches+1)],
                                               times = times, func = derivative, parms = parameters, method = "lsoda", rtol = 10**-9, atol = 10**-9))
  
  ode.out[,2:(n_spec*n_patches+1)] = ode.out[,2:(n_spec*n_patches+1)] * as.numeric(ode.out[,2:(n_spec*n_patches+1)] >= 0.01)
  return(ode.out)
}

#### dynamics_cyclic ####################################################################################################
# Returns the dynamics of a metacommunity undergoing cyclic extinctions events
# arguments :
#   Tassembly : time before first extinction
#   Text : time between each extinction
#   Tend : time after last extinction
#   n_extinctions : number of extinctions cycles
# value : out, a list of 3 elements, each the dynamics of 1 specie with time in row
# and patche in columns

dynamics_cyclic <- function(initial.state, parameters, n_extinct_patches, treatment = "autocor", Tassembly = 150, Text = 150, Tend = 150, n_extinctions = 10, by = 20){
  n_patches = parameters$n_patches
  n_spec = parameters$n_spec
  
  times <- c(seq(0, Tassembly, by = by))
  ode.out <- ode(y = initial.state, times = times, func = derivative, parms = parameters, method = "lsoda", rtol = 10**-9, atol = 10**-9)
  
  extinctions = draw.extinction.set(n_extinct_patches, treatment = treatment)
  for (spec_count in 0:(n_spec - 1)){
    ode.out[nrow(ode.out), extinctions + 1 + spec_count*n_patches] = 0
  }
  
  if (n_extinctions > 1){
    for (cycle in 1:n_extinctions-1){
      Tstart = ode.out[nrow(ode.out), 1]
      times = seq(Tstart, Tstart + Text, by = by)
      ode.out = rbind(ode.out[-nrow(ode.out),], ode(y = ode.out[nrow(ode.out),2:(n_spec*n_patches+1)],
                                                    times = times, func = derivative, parms = parameters, method = "lsoda", rtol = 10**-9, atol = 10**-9))
      extinctions = draw.extinction.set(n_extinct_patches, treatment = treatment)
      for (spec_count in 0:(n_spec - 1)){
        ode.out[nrow(ode.out), extinctions + 1 + spec_count*n_patches] = 0
      }
    }
  }
  
  Tstart = ode.out[nrow(ode.out), 1]
  times = seq(Tstart, Tstart + Tend, by = by)
  ode.out = rbind(ode.out[-nrow(ode.out),], ode(y = ode.out[nrow(ode.out),2:(n_spec*n_patches+1)],
                                                times = times, func = derivative, parms = parameters, , method = "lsoda", rtol = 10**-9, atol = 10**-9))
  
  
  
  ode.out[,2:(n_spec*n_patches+1)] = ode.out[,2:(n_spec*n_patches+1)] * as.numeric(ode.out[,2:(n_spec*n_patches+1)] >= 0.01)
  return(ode.out)
}


## plot ---------------------------------------------------------------- only working for 4*4 landscapes
plot.dynamics = function(out){
  time = out[,1]
  par(mfrow = c(4,4))
  species_name = c("Tet", "Col", "Ble")
  for (spe in 0:2){
    ymax = max(out[, 1 + 16*spe + 1:16])
    for (patch in 1:16){
      plot(time, out[,1 + 16*spe + patch], type = "l", ylab = species_name[spe + 1], main = paste("Patch", patch), ylim = c(0, ymax))
    }
  }
}

## generate.initial.state
#
# Generates an initial state in the right format for "dynamics"

generate.initial.state = function(species.map = c()){
  initial.state = rep(0, 48)
  if (length(species.map) != 16){
    species.map = sample(
      c(rep("Tet", 5), rep("Col", 5), rep("Ble", 5), sample(c("Tet", "Col", "Ble"), 1)),
      16, replace = F)
  }
  for (k in 1:16){
    if (species.map[k] == "Ble"){
      initial.state[k] = 50
    }
    if (species.map[k] == "Col"){
      initial.state[16 + k] = 500
    }
    if (species.map[k] == "Tet"){
      initial.state[32 + k] = 1000
    }
  }
  return(initial.state)
}

