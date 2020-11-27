Code and data for
"Spatial distribution of local patch extinctions drives recovery dynamics in metacommunities"
Camille Saade, Sonia Kéfi, Claire Gougat-Barbera, Benjamin Rosenbaum & Emanuel A. Fronhofer


Abstract:

Human activities lead more and more to the disturbance of plant and animal communities with local extinctions as a consequence. While these negative effects are clearly visible at a local scale it is less clear how such local patch extinctions affect regional processes, such as metacommunity dynamics and the distribution of diversity in space. Since such local extinctions may not be isolated events in space, it is crucial to investigate these questions in a spatially explicit framework, taking into account the clumping of local patch extinctions. Here, we use experimental microcosms and numerical simulations to explore the relationship between local patch extinctions and metacommunity dynamics. More specifically, we investigate the effects of patch extinction rate and spatial clumping of extinctions in a full factorial design. Experimentally, we found that local patch extinctions increased inter-patch ($\beta$-) diversity by creating differences between extinct and non-extinct patches and at the same time increased local ($\alpha$-) diversity by allowing inferior competitors to persist. Most importantly, recolonization dynamics depended more strongly on the spatial distribution of patch extinctions than on the extinction rate per se. Clumped local patch extinctions reduced mixing between extinct and non-extinct patches which led to slower recovery, lower $\alpha$-diversity in non-extinct patches and higher $\beta$-diversity. Results from a metacommunity model reproduced the experimental observations best when the model included a competition-colonization trade-off, giving a hint at the underlying mechanisms. Our results highlight that local patch extinctions can increase the diversity within and between communities, that the strength of these effects depends on the spatial distribution of extinctions and that the effects of local patch extinctions can spread regionally throughout a landscape. These findings are highly relevant for conservation and management of spatially structured communities under global change.


Folders:

- "Experiments": Data from the metacommunity experiment (species density over time) and code to conduct the statistical analyses (Tab 1, 2, S2 and S3) and plot the associated figures (Fig. 1 and 3). Data in "Experiments/Data". Run "Experiments/Script_fig1.R" and "Experiments/Script_fig3.R" to replicate the statistical analysis and plot the figures 1 and 3.

- "Lotka-Volterra_fit": Data from the single patch cultures (species density over time) and R code to fit competitive Lotka-Volterra equations to the data (used to parameterize the metacommunity model). Single patch time series are in "Lotka-Volterra_fit/Data". Run "Lotka-Volterra_fit/Model_fit.R" to conduct the fit.

- "Metacommunity_model": Code to run the metacommunity model, with a minimal working example (run "Metacommunity_model/Minimal_example.R").
