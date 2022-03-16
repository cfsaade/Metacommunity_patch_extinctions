## Script_fig4.R

## This Script is used to make the mixed models for local variables in non-extinct patches (bioarea and alpha diversity) and associated figure (Fig.3)
## in order to check whether local patch extinctions effect can spread to non-extinct patches.

## Clearing the work environment #########################################################################################  
  rm(list = ls())  


## libraries and options ###############################################################################################

  # beta diversity
  library(adespatial)
  
  # plots
  library(ggplot2)
  library(ggpubr)
  
  # data normalization
  library(bestNormalize)
  
  # mixed models and model selection
  library(lme4)
  library(MASS)
  library(MuMIn)
  
  # ggplot theme
  theme_set(theme_classic() + theme(axis.title=element_text(size=8)))

  # necessary for some models to run
  options(na.action = na.fail)
  
  # setting the number of samples for predictions of the mixed models
  n_pred = 3000
  


#### Data importation ###########################################################################################################

  ## Table with species densities over time in each patch
  density_data = read.csv("./Data/density_data.csv")
  n = nrow(density_data)
  
  ## Table with information on the experimental design (treatments, extinction positions...)
  tube_table = read.csv("./Data/TubeTable.csv")
  landscape_table = read.csv("./Data/LandscapeTable.csv")
  
  
  ## adding explicit treatment labels to the tables
  density_data$correlation_label = rep("", n)
  
  density_data$correlation_label[!density_data$correlation] = "Dispersed extinctions"
  density_data$correlation_label[density_data$correlation] = "Clustered extinctions"
  density_data$correlation_label[density_data$n_extinctions == 0] = "No extinctions"
  density_data$n_extinctions_label = paste(density_data$n_extinctions, "extinctions", sep = " ")
  density_data$n_extinctions_factor = density_data$n_extinctions == 8
  density_data$n_extinctions_label[density_data$n_extinctions == 0] = "0 extinction"
  
  density_data$extinction_label = rep("", n)
  density_data$extinction_label[density_data$extinction] = "Extinct patches"
  density_data$extinction_label[!density_data$extinction] = "Non-extinct patches"
  
  
  tube_table$n_extinctions = landscape_table$extinctions[tube_table$landscape]
  tube_table$correlation = landscape_table$cor[tube_table$landscape]
  
  tube_table$treatment_label = c(rep(0, 48), rep(1, 48), rep(2, 48), rep(3, 48), rep(4, 48))
  
  tube_table$correlation_label[!tube_table$correlation] = "Dispersed extinctions"
  tube_table$correlation_label[tube_table$correlation] = "Clustered extinctions"
  tube_table$n_extinctions_label = paste(tube_table$n_extinctions, "extinctions", sep = " ")
  tube_table$n_extinctions_label[tube_table$n_extinctions == 0] = "0 extinction"
  
  
  
  tube_table$extinction_label = rep("", 24)
  tube_table$extinction_label[tube_table$extinctions] = "Extinct patches"
  tube_table$extinction_label[!tube_table$extinctions] = "Non-extinct patches"

  # adding the distance to closest extinct patch to density data
  density_data$distance2extinct = 0
  for (id in 1:240){
    density_data$distance2extinct[density_data$tube_id == id] = tube_table$distance2extinct[tube_table$id == id]
  }
  

#### Beta diversity computation #############################################################################################################

  density_data$BetaDiv = 0  ## beta diversity
  density_data$ReplDiv = 0  ## beta diversity due to species replacement
  density_data$RichDiv = 0  ## beta diversity due to species richness differences
  
  for (landscape in unique(density_data$landscape)){
    for (mp in unique(density_data$measure_point)[-8]){
      print(mp)
      beta = beta.div.comp(density_data[density_data$landscape == landscape & density_data$measure_point == mp,
                                        c("Ble", "Tet", "Col")], quant = T)$part
      density_data$BetaDiv[density_data$landscape == landscape & density_data$measure_point == mp] = beta["BDtotal"]
      density_data$ReplDiv[density_data$landscape == landscape & density_data$measure_point == mp] = beta["Repl"]
      density_data$RichDiv[density_data$landscape == landscape & density_data$measure_point == mp] = beta["RichDif"]
    }
  }



## Subsetting the data #################################################################################################
  end_data = density_data[density_data$measure_point %in% c("t11", "t12") & (density_data$distance2extinct == 1 | density_data$n_extinctions == 0),]
  post_ext_data = density_data[density_data$measure_point %in% c("t7", "t8") & (density_data$distance2extinct == 1 | density_data$n_extinctions == 0),]


## Alpha-diversity full model ############################################################################################

  # looking for the best normalization method and normalizing data
  bestNormalize(end_data$local_simpson)
  norm_simpson = orderNorm(end_data$local_simpson)
  end_data$norm_simpson = norm_simpson$x.t
  
  # declaring the full mixed model
  mod_simpson = lmer(norm_simpson ~ as.factor(n_extinctions)*correlation_label + (1 | landscape/measure_point),
                     data = end_data,
                     control = lmerControl(optimizer = "bobyqa"),
                     REML = F)

## Bioarea full model ################################################################################################
  ## same process as above
  bestNormalize(post_ext_data$bioarea, k = 20, r = 20)
  norm_bioarea = orderNorm(post_ext_data$bioarea)
  post_ext_data$norm_bioarea = norm_bioarea$x.t
  
  mod_bioarea = lmer(norm_bioarea ~ as.factor(n_extinctions)*correlation_label + (1 | landscape/measure_point),
                     data = post_ext_data,
                     control = lmerControl(optimizer = "bobyqa"),
                     REML = F)

## Plots for model diagnostic ################################################################################################
  par(mfrow = c(2, 2))
  
  ## qqplots
  qqnorm(resid(mod_simpson), main = "qqplot Simpson")
  qqline(resid(mod_simpson))
  qqnorm(resid(mod_bioarea), main = "qqplot bioarea")
  qqline(resid(mod_bioarea))
  
  ## fitted vs residuals
  plot(x = predict(mod_simpson), y = resid(mod_simpson),
       xlab = "predicted Simpson's index", ylab = "residuals")
  lines(x = c(-1, 1), y = c(0,0))
  plot(x = predict(mod_bioarea), y = resid(mod_bioarea),
       xlab = "predicted bioarea", ylab = "residuals")
  lines(x = c(-1, 1), y = c(0,0))

    
## predictions alpha ############################################################################################################
  # In this sections we run model predictions for alpha-diversity (average of models, weighted by their AICc)
  
  # running model comparison
  mod_simpson_dredge = dredge(mod_simpson)

  # making a list of models in the order given by model dredge (model comparison)
  mod_simpson_list = list()
  mod_simpson_list[[1]] = lmer(norm_simpson ~ correlation_label * as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = end_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_simpson_list[[2]] = lmer(norm_simpson ~ correlation_label +
                                 (1 | landscape/measure_point),
                               data = end_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_simpson_list[[3]] = lmer(norm_simpson ~ correlation_label + as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = end_data,
                               control = lmerControl(optimizer = "bobyqa"))

  ## getting the models weights and coefficients 
  mod_simpson_weights = mod_simpson_dredge$weight
  mod_simpson_coefficients = lapply(mod_simpson_list, coeffs)
  mod_simpson_vcov = lapply(mod_simpson_list, vcov)
  
  ## drawing parameters from the distribution of each model
  mod_simpson_coefficients_density = list()
  for (k in 1:3){
    mod_simpson_coefficients_density[[k]] = mvrnorm(n = as.integer(n_pred * mod_simpson_weights[[k]]),
                                                    mod_simpson_coefficients[[k]], mod_simpson_vcov[[k]])
  }

  predict_simpson = function(coeff, data, model){
    ext4 = data$n_extinctions == 4
    ext8 = data$n_extinctions == 8
    corDisp = data$correlation_label == "Dispersed extinctions"
    corNoExt = data$correlation_label == "No extinctions"
    if (model == 1){
      return(coeff[1] + corDisp*coeff[2] + corNoExt*coeff[3] + ext4*coeff[4] + corDisp*ext4*coeff[5])
    }
    if (model == 2){
      return(coeff[1] + corDisp*coeff[2] + corNoExt*coeff[3])
    }
    if (model == 3){
      return(coeff[1] + corDisp*coeff[2] + corNoExt*coeff[3] + ext4*coeff[4])
    }
  }

## predictions bioarea #########################################################################################
  #same process as above
  mod_bioarea_dredge = dredge(mod_bioarea)
  
  mod_bioarea_list = list()
  
  mod_bioarea_list[[1]] = lmer(norm_bioarea ~ (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[2]] = lmer(norm_bioarea ~ as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[3]] = lmer(norm_bioarea ~ as.factor(n_extinctions) * correlation_label +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[4]] = lmer(norm_bioarea ~ correlation_label +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[5]] = lmer(norm_bioarea ~ as.factor(n_extinctions) + correlation_label +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  
  
  
  
  
  
  
  mod_bioarea_weights = mod_bioarea_dredge$weight
  mod_bioarea_coefficients = lapply(mod_bioarea_list, coeffs)
  mod_bioarea_vcov = lapply(mod_bioarea_list, vcov)
  ## drawing parameters from the distribution of each model
  mod_bioarea_coefficients_density = list()
  for (k in 1:5){
    mod_bioarea_coefficients_density[[k]] = mvrnorm(n = as.integer(n_pred * mod_bioarea_weights[k]),
                                                    mod_bioarea_coefficients[[k]], mod_bioarea_vcov[[k]])
  }
  
  predict_bioarea = function(coeff, data, model){
    ext4 = data$n_extinctions == 4
    ext8 = data$n_extinctions == 8
    corDisp = data$correlation_label == "Dispersed extinctions"
    corNoExt = data$correlation_label == "No extinctions"
    if (model == 2){ #extinction number model
      return(coeff[1] + ext4*coeff[2] + ext8*coeff[3])
    }
    if (model == 5){ # clumping + n_ext
      return(coeff[1] + ext4*coeff[2] + ext8*coeff[3] + corDisp*coeff[4])
    }
    if (model == 1){ # intercept model
      return(coeff[1])
    }
    if (model == 3){ # full model
      return(coeff[1] + ext4*coeff[2] + ext8*coeff[3] + corDisp*coeff[4] + ext4*corDisp*coeff[5])
    }
    if (model == 4){ # clumping model
      return(coeff[1] + corDisp*coeff[2] + corNoExt*coeff[3])
    }
  }
  

## Making the predictions for all variables ####################################################################################
  
  ## table for saving the predictions
  predictions_table = data.frame(n_extinctions = c(0, 4, 4, 8, 8), correlation = c(F, F, T, F, T),
                                 correlation_label = c("No extinctions", "Clustered extinctions", "Dispersed extinctions",
                                                       "Clustered extinctions", "Dispersed extinctions"))
  predictions_table$norm_simpson = NA
  predictions_table$norm_bioarea = NA
  
  new_data = predictions_table
  
  for (k in 1:(n_pred-1)){
    predictions_table = rbind(predictions_table, new_data)
  }
  
  ## predicting alpha-diversity
  n = 1
  for (model in 1:3){
    coef_list = mod_simpson_coefficients_density[[model]]
    for (i in 1:nrow(coef_list)){
      simpson = predict_simpson(coef_list[i,], new_data, model)
      predictions_table$norm_simpson[n:(n+4)] = simpson
      n = n + 5
    }
  }
  
  ## predicting bioarea 
  n = 1
  for (model in 1:5){
    coef_list = mod_bioarea_coefficients_density[[model]]
    for (i in 1:nrow(coef_list)){
      bioarea = predict_bioarea(coef_list[i,], new_data, model)
      predictions_table$norm_bioarea[n:(n+4)] = bioarea
      n = n + 5
    }
  }
  
  ## denormalizing the predictions
  predictions_table$simpson = predict(norm_simpson, predictions_table$norm_simpson, inverse = T)
  predictions_table$bioarea = predict(norm_bioarea, predictions_table$norm_bioarea, inverse = T)

## summary of predictions ################################################################################
  # making a summary of the predictions with only the median and 95% CI for all values
  summary_table = data.frame(n_extinctions = c(0, 4, 4, 8, 8), correlation = c(F, F, T, F, T),
                             correlation_label = c("No extinctions", "Clustered extinctions", "Dispersed extinctions",
                                                   "Clustered extinctions", "Dispersed extinctions"))
  
  summary_table$simpson = NA
  summary_table$simpson_min = NA
  summary_table$simpson_max = NA
  
  summary_table$bioarea = NA
  summary_table$bioarea_min = NA
  summary_table$bioarea_max = NA
  
  for (n in c(0, 4, 8)){
    for (cor in c(T, F)){
      test_pred = (predictions_table$n_extinctions == n & predictions_table$correlation == cor)
      test_summary = (summary_table$n_extinctions == n & summary_table$correlation == cor)
      
      summary_table$simpson[test_summary] = median(predictions_table$simpson[test_summary], na.rm = T)
      summary_table$simpson_min[test_summary] = quantile(predictions_table$simpson[test_summary], 0.025, na.rm = T)
      summary_table$simpson_max[test_summary] = quantile(predictions_table$simpson[test_summary], 0.975, na.rm = T)
      
      summary_table$bioarea[test_summary] = median(predictions_table$bioarea[test_summary], na.rm = T)
      summary_table$bioarea_min[test_summary] = quantile(predictions_table$bioarea[test_summary], 0.025, na.rm = T)
      summary_table$bioarea_max[test_summary] = quantile(predictions_table$bioarea[test_summary], 0.975, na.rm = T)
    }
  }
  
  summary_table$simpson_center = (summary_table$simpson_max + summary_table$simpson_min)/2
  summary_table$simpson_height = summary_table$simpson_max - summary_table$simpson_min
  
  
  summary_table$bioarea_center = (log(summary_table$bioarea_max) + log(summary_table$bioarea_min))/2
  summary_table$bioarea_height = log(summary_table$bioarea_max) - log(summary_table$bioarea_min)



## plot predictions  ###############################################################################################
  ## Making Fig. 4
  
  ## color scale
  colors = c("red", "blue", "green")
  names(colors) = levels(post_ext_data$correlation_label)
  colScale = scale_colour_manual(name = "Treatment",values = colors)
  fillScale = scale_fill_manual(name = "Treatment",values = colors)
  
  # preparing the plots
  g = ggplot(data = end_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label,
                                            color = correlation_label, group = correlation_label)) +
    colScale + fillScale
  
  h = ggplot(data = post_ext_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label,
                                                 color = correlation_label, group = correlation_label)) +
    colScale + fillScale

  
  # various graphical parameters
  alpha_box = 0.3 ## boxes transparency
  alpha_points = 0.5
  box_width = 0.4
  dodge_width = 0.75


  panel_b = g +
    geom_jitter(aes(y = local_simpson), position=position_jitterdodge(jitter.width = 0.1), alpha = alpha_points, size = 0.5) +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = simpson_center,
                                        width = box_width, height = simpson_height),
              position=position_dodge(width = dodge_width), alpha = alpha_box, color = NA) +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = simpson,
                                        width = box_width, height = 0.005),
              position=position_dodge(width = dodge_width), alpha = 1, color = NA) +
    labs(y = expression(paste(alpha, "-diversity (Simpson's index)")),
         x = "Extinctions",
         fill = "Treatment",
         color = "Treatment") + guides(color=guide_legend(override.aes=list(fill=NA))) +
    theme(legend.title = element_text(size = 7), legend.text=element_text(size = 7))
  
  panel_a = h +
    geom_jitter(aes(y=log(bioarea)), position = position_jitterdodge(jitter.width = 0.1), alpha = alpha_points, size = 0.5) +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = bioarea_center,
                                        width = box_width, height = bioarea_height),
              position=position_dodge(width = dodge_width), alpha = alpha_box, color = NA) +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = log(bioarea),
                                        width = box_width, height = 0.01),
              position=position_dodge(width = dodge_width), alpha = 1, color = NA) +
    labs(y = expression(paste("Biomass (log(",mu,m^2,"/", mL, "))", sep = "")),
         x = "Extinctions",
         fill = "Treatment",
         color = "Treatment") + guides(color=guide_legend(override.aes=list(fill=NA))) +
    theme(legend.title = element_text(size = 7), legend.text=element_text(size = 7))


  fig4 = ggarrange(panel_a, panel_b,
          labels = "auto",
          common.legend = TRUE, legend = "bottom", nrow = 2) +
    theme(legend.title = element_text(size = 7), legend.text=element_text(size = 7))

  
  ggsave("./fig4.pdf",
         fig4,
         width = 110,
         height = 169,
         unit = "mm")
  
  
