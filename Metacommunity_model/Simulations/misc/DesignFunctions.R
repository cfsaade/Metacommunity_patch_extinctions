library(raster)
library(landscapemetrics)
library(pracma)

# position2rowcol ------------------------------  done
# Transforms the position (int in 1:size^2) to
# row and column numbers in a size*size grid

position2rowcol = function(position, size = 4){
  row = ceiling(position/size)
  column = (position-1)%%size + 1
  return(c(row, column))
}


# nlink -----------------------------------------  done
# defininit la fonction nlink(position, size, border = T)
# prend en argument la position d'un point dans
# une grille de size*size (position dans 1:size*size,
# numerotation de gauche a droite et de haut en
# bas)

nlink = function(position, size = 4, border = T){
  if (size == 4 & border == T){
    if (position %in% c(1, 4, 13, 16)){
      return(2)
    }
    else if (position %in% c(6, 7, 10, 11)){
      return(4)
    }
    else return(3)
  }
  if (border == F){
    return(4)
  }
  if (size != 4 & border == T){
    pos.rowcol = position2rowcol(position, size = size)
    return(4 - (pos.rowcol[1] %in% c(1, size)) - (pos.rowcol[2] %in% c(1, size)))
  }
}

# nlink.set ------------------------------------  done
# returns the number of links of a points set
# (the sum of nlink for each point of the set)
# in a size*size grid

nlink.set = function(points.set, size = 4, border = T){
  return(sum(sapply(points.set, nlink, size = size, border = border)))
}



# grid.distance --------------------------------  done
# computes the distance (number of links) 
# between two points in a size*size grid
#
# border indicates if there is a border (TRUE)
# or if the grid is a torus (FALSE)

grid.distance = function(pos1, pos2, size = 4, border = T){
  pos1.rowcol = position2rowcol(pos1, size = size)
  pos2.rowcol = position2rowcol(pos2, size = size)
  
  dist.vector = abs(pos1.rowcol - pos2.rowcol)
  if (border == F & dist.vector[1] == size - 1){
    dist.vector[1] = 1
  }
  if (border == F & dist.vector[2] == size - 1){
    dist.vector[2] = 1
  }
  return(sum(dist.vector))
}


# overdispersion -------------------------------- done
# computes the mean distance between points
# in a set of points (measure of dispertion)

overdispersion = function(points.set){
  dist = 0
  n = length(points.set)
  for (k in 1:n){
    dist = dist + sum(sapply(X = points.set[-k], FUN = grid.distance, pos2 = points.set[k]))
  }
  return(dist/(n*(n-1)))
}

# plot.set -------------------------------------- not working on anything else than a 4*4 grid yet
# plots a set of points
plot.set = function(points.set, main = ''){
  coords = sapply(points.set, position2rowcol)
  plot(rep(1:4, 4), c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), pch = 19, ylim = c(4.5,0.5), xlim = c(0.5,4.5),
       xlab = '', ylab = '', cex = 2, main = main)
  points(coords[2,], coords[1,], pch = 19, col = "red", cex = 2)
}

# draw.random.set ------------------------------------  obsolete, not used anymore
# returns a random set drawn among the sets
# with the same connectivity as the target.set
# and the highest dispersion (lowest pladj)
draw.random.set = function(target.set){
  target.nlink = nlink.set(target.set)
  candidates = list()   # list with random sets of the same nlink as target
  n = 1
  for (k in 1:10000){
    newset = sample(1:16, length(target.set), F)
    if (nlink.set(newset) == target.nlink){
      candidates[[n]] = newset
      n = n + 1
    }
  }
  
  candidates = lapply(candidates, sort)
  
  candidates = unique(candidates)
  
  candidates.disp = sapply(candidates, percent.like.adjacencies)
  
  candidates = candidates[candidates.disp == max(candidates.disp)]
  return(sample(candidates, 1)[[1]])
}

# distance2nonextinct -----------------------------------------------  done
# returns a list of distances between each patch of the set
# and the nearest non-extinct patch
#
distance2nonextinct = function(points.set, size = 4, border = T){
  dist.list = c()
  for (k in points.set){
    dist.list = c(dist.list, min(sapply(c(1:(size*size))[-points.set], FUN =  grid.distance, pos2 = k, size = size, border = border)))
  }
  return(dist.list)
}


# connectivity2nonexctinct ------------------------------------------  done
# returns the number of adjacent non-extinct patches
# for each patch in the set
#
connectivity2nonextinct = function(points.set, size = 4, border = T){
  connectivity.list = c()
  for (k in points.set){
    connectivity.list = c(connectivity.list, sum((sapply(c(1:(size*size))[-points.set], FUN =  grid.distance, pos2 = k, size = size, border = border) == 1)))
  }
  return(connectivity.list)
}

# percent.like.adjacencies ---------------------------------------------------- done
percent.like.adjacencies = function(points.set, size = 4){
  a = rep(0, size*size)
  a[points.set] = 1
  landscape = raster(matrix(a, nrow = size))
  return(lsm_c_pladj(landscape)$value[2])
}

#### draw.extinction.set ##############################################################################################
# draws a random set of patches with average connectivity == 3
# clumped set if treatment == autocor (default)
# dispersed set if treatment != autocor

draw.extinction.set = function(n_extinctions, treatment = "autocor"){
  
  candidates = list()  ## list of possible sets to draw from
  
  if(n_extinctions == 0){  ## no extinction -> empty set
    return(c())
  }
  
  matrix.set = combs(1:16, n_extinctions) ## getting all the subsets of "n_extinctions" patches
  candidates = list()
  for (i in 1:length(matrix.set[,1])){
    candidates[[i]] = matrix.set[i,]
  }
  
  candidates.nlink = sapply(candidates, nlink.set)
  candidates = candidates[candidates.nlink == 3*n_extinctions]  ## filtering all sets with average connectivity != 3
  
  
  
  candidates = lapply(candidates, sort)  ## removing duplicates
  candidates = unique(candidates)
  
  
  candidates.pladj = sapply(candidates, percent.like.adjacencies) ## pladj is a measure of spatial autocorrelation
  
  if (treatment == 'autocor'){
    candidates = candidates[candidates.pladj == max(candidates.pladj)]
  }
  if (treatment != 'autocor'){
    candidates = candidates[candidates.pladj == min(candidates.pladj)]
  }
  return(sample(candidates, 1)[[1]])
}


## generate.neigh.matrix --------------------------------------- done
#
# generates the neighbour matrix of a size*size grid
# columns normalized to sum to 1

generate.neigh.matrix = function(size = 4, border = T){
  npatches = size*size
  neigh = matrix(rep(0, npatches**2), nrow = npatches)
  for (k in 1:npatches){
    for (l in 1:npatches){
      if (grid.distance(k, l, size = size, border = border) == 1){
        neigh[k, l] = 1
        neigh[l, k] = 1
      }
    }
  }
  for (k in 1:npatches){
    neigh[,k] = neigh[,k]/sum(neigh[,k])
  }
  return(neigh)
}

