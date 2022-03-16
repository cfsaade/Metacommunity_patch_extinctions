### Full design
## creates a landscape table and a tube table that contain the extinction maps and initial species for
## extinctions ranging from 2 to 15 and 2 correlation regimes
rm(list = ls())
source("DesignFunctions.R")

n_rep = 100
n_landscapes = 5*n_rep
n_patches = n_landscapes * 16
## table of landscapes properties
landscape.tab = data.frame(id = 1:n_landscapes)                                          # Numerotation (150 landscapes)

landscape.tab$extinctions = rep(c(0, 4, 4, 8, 8), n_rep)                                 # Number of patches going extinct in
landscape.tab$replicate = rep(1:n_rep, each = 5)                                                 # each lanscape

landscape.tab$cor = rep(c(F, F, T, F, T), n_rep)                     # Color of the extinctions :
# T = autocorrelation
# F = no autocorrelation



## Table of the patches
tube.tab = data.frame(id = c(1:n_patches))  # numerotation du tube parmi tous les tubes (utile pour les videos)

tube.tab$landscape = ceiling(tube.tab$id/16) # number of the landscape (1 to 150)

tube.tab$position = rep(1:16, n_landscapes)     # la position dans le patch, 1 en haut a gauche puis 
# numerotation de gauche a droite et de haut en bas

tube.tab$row = ceiling(tube.tab$position/4)
tube.tab$column = (tube.tab$position-1)%%4 + 1

## Plans des espèces initiales (sur une base de 4 esp. , a refaire apres choix des especes)
# 1 plan pour tous les landscapes (replicats reels) ? 1 plan par landscape ?
tube.tab$initial.species = rep(0, n_patches)
for (i in 1:n_landscapes){
  tube.tab$initial.species[tube.tab$landscape == i] = sample(x = c(rep(c("Tet", "Col", "Ble"), 5), sample(c("Tet", "Col", "Ble"), 1)),
                                                             size = 16, replace = F)
}


## Plans des extinctions

tube.tab$extinctions = rep(F, n_patches)
tube.tab$distance2nonextinct = rep(0, n_patches)
tube.tab$connectivity2nonextinct = rep(NA, n_patches)


landscape.tab$maxdistance2nonextinct = rep(0, n_landscapes)
landscape.tab$averagedistance2nonextinct = rep(0, n_landscapes)
landscape.tab$pladj = rep(0, n_landscapes)
landscape.tab$averageconnectivity2nonextinct = rep(0,n_landscapes)

for (landscape in 1:n_landscapes){
  print(landscape)
  cor = ifelse(landscape.tab$cor[landscape.tab$id == landscape], "autocor", "anticor")
  n_ext = landscape.tab$extinctions[landscape.tab$id == landscape]
  
  set = draw.extinction.set(n_ext, cor)
  distance = distance2nonextinct(set)
  connectivity = connectivity2nonextinct(set)
  
  tube.tab$extinctions[tube.tab$landscape == landscape][set] = T
  tube.tab$distance2nonextinct[tube.tab$landscape == landscape][set] = distance
  tube.tab$connectivity2nonextinct[tube.tab$landscape == landscape][set] = connectivity
  
  landscape.tab$maxdistance2nonextinct[landscape] = max(distance)
  landscape.tab$averagedistance2nonextinct[landscape] = mean(distance)
  landscape.tab$pladj[landscape] = percent.like.adjacencies(set)
  landscape.tab$averageconnectivity2nonextinct[landscape] = mean(connectivity2nonextinct(set))
}


## saving the tables in csv
write.csv(tube.tab, file = "SimulationsTubeTable.csv")
write.csv(landscape.tab, file = "SimulationsLandscapeTable.csv")


