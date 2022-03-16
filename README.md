Code and data for
"Spatial distribution of local patch extinctions drives recovery dynamics in metacommunities"
Camille Saade, Sonia Kéfi, Claire Gougat-Barbera, Benjamin Rosenbaum & Emanuel A. Fronhofer


Abstract:

Human activities put ecosystems under increasing pressure, often resulting in local extinctions. However, it is unclear how local extinctions affect regional processes, such as the distribution of diversity in space, especially if extinctions show spatial patterns, such as being clustered. Therefore, it is crucial to investigate extinctions and their consequences in a spatially explicit framework. Using highly controlled microcosm experiments and theoretical models, we here ask how the number and spatial autocorrelation of extinctions interactively affect metacommunity dynamics. We found that local patch extinctions increased local ($\alpha$-) and inter-patch ($\beta$-) diversity by delaying the exclusion of inferior competitors. Importantly, recolonization dynamics depended more strongly on the spatial distribution than on the number of patch extinctions: clustered local patch extinctions resulted in slower recovery, lower $\alpha$-diversity and higher $\beta$-diversity. Our results highlight that the spatial distribution of perturbations should be taken into account when studying and managing spatially structured communities.

Folders:
	- "Experiments": Data from the metacommunity experiment (species density over time) and code to conduct the statistical analyses (Tab 1, 2, S2 and S3) and plot the associated figures (Fig. 2 and 4). Data in "Experiments/Data". Run "Experiments/Script_fig2.R" and "Experiments/Script_fig4.R" to replicate the statistical analysis and plot the figures 2 and 4.

	- "Lotka-Volterra_fit": Data from the single patch cultures (species density over time) and R code to fit competitive Lotka-Volterra equations to the data (used to parameterize the metacommunity model). Single patch time series are in "Lotka-Volterra_fit/Data". Run "Lotka-Volterra_fit/Model_fit.R" to conduct the fit.

	- "Metacommunity_model": Code to run the metacommunity model, with a minimal working example (folder "Metacommunity_model/Minimal_working_example"), the code to replicate the main simulation of the text (folder "Metacommunity_model/Simulations"), and code to replicate 10% of these simulations with another integrator and low error tolerance (solver 'lsoda' with rtol = atol = 10^-9 instead of solver 'ode45') to check that the results are free of integration error (folder "Metacommunity_model/Simulations_low_error_tolerance").
