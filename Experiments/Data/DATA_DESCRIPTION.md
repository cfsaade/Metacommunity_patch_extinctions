
# density_data.csv
A table with the densities of each species over time in the experimental metacommunities, as well as treatments informations. Values in columns:

- "file": the name of the video sample used for species count and identification
- "date": date at which the sample was taken (yyyy-mm-dd)
- "time": the time of day at which the sample was taken (hour)
- "magnification": magnification setting of the microscope during video acquisition
- "landscape": a number (1-15) indicating from which landscape the sample was taken
- "tube_id": a number (1-240) indicating from which tube the sample was taken
- "position": a number (1-16) indicating the position of the tube in its landscape (numbered from left to right and top to bottom)
- "measure_point": the number of the measurement (t0 to t11)
- "initial_species": the species that was initially seeded into the tube ("Ble", "Tet" or "Col")
- "extinctions": indicates if the tube is an extinct (TRUE) or non-extinct (FALSE) patch
- "n_extinctions": the number of extinctions in this tube's landscape (0, 4, 8)
- "correlation": the spatial clumping of extinctions in this tube's landscape (TRUE: clumped extinctions, FALSE: dispersed extinctions/no extinctions)
- "bioarea": the bioarea per volume in the video sample (micrometer^2/mL, proxi of biomass in the tube)
- "Tet": the density (individuals/mL) of Tetrahymena thermophila
- "Col": the density (individuals/mL) of Colpidium sp.
- "Ble": the density (individuals/mL) of Blepharisma sp.
- "Ble_shadow": density of particles identified as the shadow of a Blepharisma sp. individual (indiv/mL)
- "Unassigned": particles not assigned to any species because of low indentification confidence (indiv/mL)
- "local_simpson": alpha-diversity (measured with the Simpson's index) inside the tube
- "landscape_simpson": Simpson's index computed at the scale of the landscape
- "hours": the time passed since the begining of the experiment (in hours)


# LandscapeTable.csv
A table with landscape-level informations on the experimental design. Each line corresponds to one landcape. Values in columns:

- "id": a unique number identifying each landscape (1-15)
- "extinctions": the number of extinctions in the landscape
- "cor": spatial treatment in said landscape (TRUE: clumped extinctions, FALSE: dispersed extinctions/no extinction)
- "maxdistance2nonextinct": the maximal distance between an extinct patch and the closest non-extinct patch
- "averagedistance2nonextinct": the average distance between extinct patches and the closest non-extinct patch
- "pladj": "percentage of like adjacencies", a measure of the spatial clumping of extinctions
- "averageconnectivity2nonextinct": the average connectivity (number of links) between extinct patches and non-extinct patches.

# TubeTable.csv
A table with tube-level informations on the experimental design. Each line corresponds to one tube. Values in columns:

- "id": a unique number identifying each tube (1-240)
- "landscape": the landscape in which this tube is (1-15)
- "position": a number (1-16) indicating the position of the tube in its landscape (numbered from left to right and top to bottom)
- "row"/"column": the position of the tube in its landscape expressed as row and column coordinates
- "initial_species": the species that was initially seeded into the tube ("Ble", "Tet" or "Col")
- "extinctions": indicates if the tube is an extinct (TRUE) or non-extinct (FALSE) patch
- "distance2nonextinct": distance to the nearest non-extinct patch
- "connectivity2nonextinct": Number of ajdacent non-extinct patches
- "n_extinctions": the number of extinctions in the tube's landscape
- "correlation": spatial treatment in the tube's landscape (TRUE: clumped extinctions, FALSE: dispersed extinctions/no extinction)
- "treatment_label": a numeric label (0-4) indicating the different combination of rate and clumping treatments (useful for subsetting the data)
- "correlation_label": a character string label for the spatial clumping treatement (same info as correlation) to make figures labels
- "n_extinctions_label": a character string label for the extinction rate treatement (same info as "n_extinctions") to make figures labels
- "extinction_label": a character string label indicating if the tube is an extinct patches (useful for figures labels)
- "treatment_label_bis": a character string label indicating he different combination of rate and clumping treatments (useful for figure labels)
- "return_point": the measurement point at which extinct patches come back to 95% of the pre-extinction bioarea
- "return_time": the time extinct patch takes to come back to 95% of the pre-extinction bioarea (in hours)
- "distance2extinct": the distance to the closest extinct patch



