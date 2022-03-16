## parameters 


## Creating the neighbourhood matrix ############
neigh = generate.neigh.matrix(4, T)

generate.parameters = function(n_spec = 3, scenario = "neutral", migration_rate = 1/100){
  if (scenario == "fit"){ ### values fitted from the experiment (all except migration)
    r_vect = c(3.693081e-02, 5.381216e-02, 1.778028e-01)
    r_vect = r_vect*(r_vect > 0)
    
    migration_vect = rep(migration_rate, n_spec)
    
    comp_vect = c(3.413686e-04, 8.149289e-04, 1.885881e-03,
                  1.498016e-05, 2.326073e-05, 1.335655e-04,
                  9.374581e-07, 5.388102e-06, 3.661246e-05)
    comp_vect = comp_vect*(comp_vect > 0)
    
    comp_matrix = matrix(comp_vect, n_spec, byrow = T)
    
    param = list(size = 4,                      # number of rows = number of columns
                 n_patches = 16,                 # total number of patches (= size*size)
                 n_spec = n_spec,
                 
                 migration_vector = migration_vect,
                 competition_matrix = comp_matrix,
                 r = r_vect,
                 neighbour_matrix = neigh)
    return(param)
  }
  
  if (scenario == "trade_off"){ ### values fitted from the experiment (all except migration)
    r_vect = c(3.693081e-02, 5.381216e-02, 1.778028e-01)
    r_vect = r_vect*(r_vect > 0)
    
    migration_vect = c(migration_rate/5, migration_rate, 2*migration_rate)
    migration_vect = migration_vect*(migration_vect > 0)
    
    comp_vect = c(3.413686e-04, 8.149289e-04, 1.885881e-03,
                  1.498016e-05, 2.326073e-05, 1.335655e-04,
                  9.374581e-07, 5.388102e-06, 3.661246e-05)
    comp_vect = comp_vect*(comp_vect > 0)
    
    comp_matrix = matrix(comp_vect, n_spec, byrow = T)
    
    param = list(size = 4,                      # number of rows = number of columns
                 n_patches = 16,                 # total number of patches (= size*size)
                 n_spec = n_spec,
                 
                 migration_vector = migration_vect,
                 competition_matrix = comp_matrix,
                 r = r_vect,
                 neighbour_matrix = neigh)
    return(param)
  }
  
  if (scenario == "no_interaction"){
    r_vect = c(3.693081e-02, 5.381216e-02, 1.778028e-01)
    r_vect = r_vect*(r_vect > 0)
    
    migration_vect = rep(migration_rate, n_spec)
    
    comp_vect = c(3.413686e-04, 0,            0,
                  0,            2.326073e-05, 0,
                  0,            0,            3.661246e-05)
    comp_vect = comp_vect*(comp_vect > 0)
    
    comp_matrix = matrix(comp_vect, n_spec, byrow = T)
    
    param = list(size = 4,                      # number of rows = number of columns
                 n_patches = 16,                 # total number of patches (= size*size)
                 n_spec = n_spec,
                 
                 migration_vector = migration_vect,
                 competition_matrix = comp_matrix,
                 r = r_vect,
                 neighbour_matrix = neigh)
    return(param)
  }
  
  if (scenario == "randomized"){
    r_vect = c(3.693081e-02, 5.381216e-02, 1.778028e-01)
    r_vect = r_vect*(r_vect > 0)
    
    migration_vect = rep(migration_rate, n_spec)
    
    comp_vect = c(3.413686e-04, 8.149289e-04, 1.885881e-03,
                  1.498016e-05, 2.326073e-05, 1.335655e-04,
                  9.374581e-07, 5.388102e-06, 3.661246e-05)
    
    diag_index = sample(c(1, 2, 3))
    
    off_diag = comp_vect[c(2, 3, 6)]
    off_diag = off_diag[diag_index]
    comp_vect[c(2, 3, 6)] = off_diag
    
    off_diag = comp_vect[c(4, 7, 8)]
    off_diag = off_diag[diag_index]
    comp_vect[c(4, 7, 8)] = off_diag
    
    
    comp_matrix = matrix(comp_vect, n_spec, byrow = T)
    
    param = list(size = 4,                      # number of rows = number of columns
                 n_patches = 16,                 # total number of patches (= size*size)
                 n_spec = n_spec,
                 
                 migration_vector = migration_vect,
                 competition_matrix = comp_matrix,
                 r = r_vect,
                 neighbour_matrix = neigh)
    return(param)
  }
  
  if (scenario == "shuffled"){
    ## shuffle all the interspecific interaction terms
    r_vect = c(3.693081e-02, 5.381216e-02, 1.778028e-01)
    r_vect = r_vect*(r_vect > 0)
    
    migration_vect = rep(migration_rate, n_spec)
    
    comp_vect = c(3.413686e-04, 8.149289e-04, 1.885881e-03,
                  1.498016e-05, 2.326073e-05, 1.335655e-04,
                  9.374581e-07, 5.388102e-06, 3.661246e-05)
    
    diag_index = sample(c(1, 2, 3))
    
    off_diag = comp_vect[c(2, 3, 6, 4, 7, 8)]
    off_diag = sample(off_diag)
    comp_vect[c(2, 3, 6, 4, 7, 8)] = off_diag
    
    
    comp_matrix = matrix(comp_vect, n_spec, byrow = T)
    
    param = list(size = 4,                      # number of rows = number of columns
                 n_patches = 16,                 # total number of patches (= size*size)
                 n_spec = n_spec,
                 
                 migration_vector = migration_vect,
                 competition_matrix = comp_matrix,
                 r = r_vect,
                 neighbour_matrix = neigh)
    return(param)
  }
  
}
