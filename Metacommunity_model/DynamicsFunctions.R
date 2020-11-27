##  DynamicsFunctions.R


# Scripts related to simulating competitive lotka-volterra dynamics in a metacommunity


## Libraries #############################################################################################################

  # ode solver
  library(deSolve)

## derivative ############################################################################################################
  # Returns the derivative of a 3 species system (X, Y, Z) in a metacommunity of size*size patches (R equivalent of Eq. 1-4)
  #
  # Parameters
  #   t : time, unused but necessary for ode function
  #   state : the poplulation densities in each of the 16 patches as a 48 value vector (1:16 = X, 17:32 = Y, 33:48 = Z)
  
  
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
# and patches in columns

dynamics <- function(initial.state, parameters, extinctions, Thalf = 500, Tmax = 1000, by = 20){
  n_patches = parameters$n_patches
  n_spec = parameters$n_spec
  times <- c(seq(0, Thalf, by = by), Thalf + 0.1)
  Thalf_index = length(times)
  ode.out <- ode(y = initial.state, times = times, func = derivative, parms = parameters, method = "ode45")
  for (spec_count in 0:(n_spec - 1)){
    ode.out[Thalf_index, extinctions + 1 + spec_count*n_patches] = 0
  }
  times = seq(Thalf + 0.1, Tmax + 0.1, by = by)
  ode.out = rbind(ode.out[1:Thalf_index,], ode(y = ode.out[Thalf_index,2:(n_spec*n_patches+1)],
                                               times = times, func = derivative, parms = parameters, , method = "ode45"))
  
  ode.out[,2:(n_spec*n_patches+1)] = ode.out[,2:(n_spec*n_patches+1)] * as.numeric(ode.out[,2:(n_spec*n_patches+1)] >= 0.01)
  return(ode.out)
}


## generate.initial.state ##########################################################################################################
#
# A utility function providing an initial state in the right format for the "dynamics" functions
#   either from a species.map (a vector of length 16 with abbreviated species names ("Tet", "Ble", "Col") as elements)
#   or randomly if species.map() is not provided

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

