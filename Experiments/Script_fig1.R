## Script_fig1.R
  
## This Script is used to make the mixed models for local variables in extinct patches (bioarea, recovery time and alpha diversity) and for beta-diversity,
## as well as the associated figures (all panels of Fig. 1)

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
  theme_set(theme_classic(base_size = 14))
  
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
  
  
## recovery time computation #############################################################################################################
      ## distribution of bioarea before perturbation
  steady_state_bioarea = density_data$bioarea[density_data$measure_point %in% c("t1", "t2", "t3", "t4", "t5", "t6")]
  
  
  ## setting the threshold for a return time as the 2.5% percentile
  bioarea_threshold = quantile(steady_state_bioarea, 0.025)
  tube_table$return_point = rep(0, 240)
  tube_table$return_time = rep(0, 240)
  for (tube_id in unique(density_data$tube_id[density_data$extinction == T])){
    test_vector = density_data$bioarea[density_data$tube_id == tube_id & density_data$hours >= 386]
    return_point = min(which(test_vector > bioarea_threshold))
    tube_table$return_point[tube_table$id == tube_id] = return_point
    tube_table$return_time[tube_id] = mean(density_data$hours[density_data$measure_point == paste("t", return_point + 6, sep = "")]) - 344
  }
  tube_table$return_time[tube_table$return_point == Inf] = 260
  
  tube_table$return_point[tube_table$return_point == Inf] = 6
  
  
## Beta diversity computation #############################################################################################################
  
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
  
  
  

## Subsetting the data ##################################################################################################
  ## we use only the data from after the extinctions in extinct patches for models
  
  post_ext_data = density_data[density_data$hours > 345 & density_data$n_extinctions != 0 & density_data$extinction == T,]
  post_ext_landscape_data = density_data[density_data$hours > 345 & density_data$n_extinctions != 0 & density_data$position == 1,]
  
  # all post-extinctions data for the figures
  fig_data = density_data[density_data$hours > 345 & (density_data$n_extinctions == 0 | density_data$extinction == T),]
  fig_data$correlation_label[fig_data$n_extinctions == 0] = "No extinction"
  
  fig_data_landscape = density_data[density_data$hours > 345 & density_data$position == 1,]
  fig_data_landscape$correlation_label[fig_data_landscape$n_extinctions == 0] = "No extinction"

  

## Alpha-diversity full model ############################################################################################
  ## looking for best normalization method
  bestNormalize(post_ext_data$bioarea)
  
  ## normalizing data
  norm_simpson = orderNorm(post_ext_data$local_simpson)
  post_ext_data$norm_simpson = norm_simpson$x.t
  
  ## mixed model declaration
  mod_simpson = lmer(norm_simpson ~ as.factor(n_extinctions)*as.factor(correlation) + (1 | landscape/measure_point),
                     data = post_ext_data,
                     control = lmerControl(optimizer = "bobyqa"),
                     REML = F)

## Beta-diversity full model #################################################################################################
  ## same process as above
  
  beta_data = na.omit(post_ext_landscape_data)
  bestNormalize(beta_data$BetaDiv) ## no transformation necessary
  mod_beta = lmer(BetaDiv~ as.factor(n_extinctions)*as.factor(correlation) + (1 | measure_point),
                  data = beta_data,
                  control = lmerControl(optimizer = "bobyqa"),
                  REML = F)
  
## Bioarea full model #######################################################################################################
  ## same process as above
  
  bestNormalize(post_ext_data$bioarea)
  norm_bioarea = orderNorm(post_ext_data$bioarea)
  post_ext_data$norm_bioarea = norm_bioarea$x.t
  
  mod_bioarea = lmer(norm_bioarea ~ as.factor(n_extinctions)*as.factor(correlation) + (1 | landscape/measure_point),
                     data = post_ext_data,
                     control = lmerControl(optimizer = "bobyqa"),
                     REML = F)

  
  
## Return Time full model ###############################################################################################
  ## same process as above
  
  return_time_data = tube_table[tube_table$extinctions & tube_table$return_point != Inf,]
  
  bestNormalize(return_time_data$return_time)
  norm_return_time = orderNorm(return_time_data$return_time)
  return_time_data$norm_return_time = norm_return_time$x.t
  
  mod_return_time = lmer(norm_return_time ~ as.factor(n_extinctions)*as.factor(correlation) + (1 | landscape),
                         data = return_time_data,
                         control = lmerControl(optimizer = "bobyqa"),
                         REML = F)

  
## Plots for model diagnostic ################################################################################################
  par(mfrow = c(2, 2))
  
  ## qqplots
    qqnorm(resid(mod_simpson), main = "qqplot Simpson")
    qqline(resid(mod_simpson))
    qqnorm(resid(mod_beta), main = "qqplot beta")
    qqline(resid(mod_beta))
    qqnorm(resid(mod_bioarea), main = "qqplot bioarea")
    qqline(resid(mod_bioarea))
    qqnorm(resid(mod_return_time), main = "qqplot recovery time")
    qqline(resid(mod_return_time))
  
  ## fitted vs residuals
    plot(x = predict(mod_simpson), y = resid(mod_simpson),
         xlab = "predicted Simpson's index", ylab = "residuals")
    lines(x = c(0, 10), y = c(0,0))
    plot(x = predict(mod_beta), y = resid(mod_beta),
         xlab = "predicted beta diversity", ylab = "residuals")
    lines(x = c(0, 10), y = c(0,0))
    plot(x = predict(mod_bioarea), y = resid(mod_bioarea),
         xlab = "predicted bioarea", ylab = "residuals")
    lines(x = c(0, 12000), y = c(0,0))
    plot(x = predict(mod_return_time), y = resid(mod_return_time),
         xlab = "predicted recovery time", ylab = "residuals")
    lines(x = c(0, 300), y = c(0,0))
  
  
## predictions alpha ############################################################################################################
  # In this sections we run model predictions for alpha-diversity (average of models, weighted by their AICc)
    
  ## running model comparison
  mod_simpson_dredge = dredge(mod_simpson)
  
  # making a list of models in the order given by model dredge (model comparison)
  mod_simpson_list = list()
  mod_simpson_list[[1]] = lmer(norm_simpson ~ as.factor(correlation) +
                                 (1 | landscape/measure_point),
                     data = post_ext_data,
                     control = lmerControl(optimizer = "bobyqa"))
  
  mod_simpson_list[[2]] = lmer(norm_simpson ~ as.factor(correlation) + as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_simpson_list[[3]] = lmer(norm_simpson ~ as.factor(correlation)*as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_simpson_list[[4]] = lmer(norm_simpson ~ (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_simpson_list[[5]] = lmer(norm_simpson ~ as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  ## getting the models weights and coefficients 
  mod_simpson_weights = mod_simpson_dredge$weight
  mod_simpson_coefficients = lapply(mod_simpson_list, coeffs)
  mod_simpson_vcov = lapply(mod_simpson_list, vcov)
  
  ## drawing parameters from the distribution of each model
  mod_simpson_coefficients_density = list()
  for (k in 1:5){
    mod_simpson_coefficients_density[[k]] = mvrnorm(n = as.integer(n_pred * mod_simpson_weights[k]),
                                             mod_simpson_coefficients[[k]], mod_simpson_vcov[[k]])
  }
  
  predict_simpson = function(coeff, data, model){
    ext = data$n_extinctions == 8
    cor = data$correlation
    if (model == 1){
      return(coeff[1] + cor*coeff[2])
    }
    if (model == 2){
      return(coeff[1] + cor*coeff[2] + ext*coeff[3])
    }
    if (model == 3){
      return(coeff[1] + cor*coeff[2] + ext*coeff[3] + (ext*cor)*coeff[4])
    }
    if (model == 4){
      dummy = rep(1, nrow(data))
      return(dummy*coeff[1])
    }
    if (model == 5){
      return(coeff[1] + ext*coeff[2])
    }
  }
  
## predict beta ###########################################################################################
  ## same process as above
  mod_beta_dredge = dredge(mod_beta)
  
  mod_beta_list = list()
  mod_beta_list[[1]] = lmer(BetaDiv ~ as.factor(correlation)*as.factor(n_extinctions) +
                              (1 | measure_point),
                            data = beta_data,
                            control = lmerControl(optimizer = "bobyqa"))
  
  mod_beta_weights = mod_beta_dredge$weight
  mod_beta_coefficients = lapply(mod_beta_list, coeffs)
  mod_beta_vcov = lapply(mod_beta_list, vcov)
  ## drawing parameters from the distribution of each model
  mod_beta_coefficients_density = list()
  for (k in 1){
    mod_beta_coefficients_density[[k]] = mvrnorm(n = n_pred, mod_beta_coefficients[[k]], mod_beta_vcov[[k]])
  }
  
  predict_beta = function(coeff, data, model){
    ext = data$n_extinctions == 8
    cor = data$correlation
      return(coeff[1] + cor*coeff[2] + ext*coeff[3] + (ext*cor)*coeff[4])
  }
  
## predict bioarea ########################################################################
  # same process as above
  mod_bioarea_dredge = dredge(mod_bioarea)
  
  mod_bioarea_list = list()
  mod_bioarea_list[[1]] = lmer(norm_bioarea ~ as.factor(correlation) +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[2]] = lmer(norm_bioarea ~ (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[3]] = lmer(norm_bioarea ~ as.factor(correlation) + as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[4]] = lmer(norm_bioarea ~ as.factor(n_extinctions) +
                                 (1 | landscape/measure_point),
                               data = post_ext_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_bioarea_list[[5]] = lmer(norm_bioarea ~ as.factor(correlation)*as.factor(n_extinctions) +
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
    ext = data$n_extinctions == 8
    cor = data$correlation
    if (model == 1){
      return(coeff[1] + cor*coeff[2])
    }
    if (model == 3){
      return(coeff[1] + cor*coeff[2] + ext*coeff[3])
    }
    if (model == 5){
      return(coeff[1] + cor*coeff[2] + ext*coeff[3] + (ext*cor)*coeff[4])
    }
    if (model == 2){
      dummy = rep(1, nrow(data))
      return(dummy*coeff[1])
    }
    if (model == 4){
      return(coeff[1] + ext*coeff[2])
    }
  }
  
## predict return time ##########################################################################
  # same process as above
  mod_return_time_dredge = dredge(mod_return_time)
  
  mod_return_time_list = list()
  mod_return_time_list[[1]] = lmer(norm_return_time ~ as.factor(correlation) +
                                 (1 | landscape),
                               data = return_time_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_return_time_list[[2]] = lmer(norm_return_time ~ (1 | landscape),
                               data = return_time_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_return_time_list[[3]] = lmer(norm_return_time ~ as.factor(correlation) + as.factor(n_extinctions) +
                                 (1 | landscape),
                               data = return_time_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_return_time_list[[4]] = lmer(norm_return_time ~ as.factor(n_extinctions) +
                                 (1 | landscape),
                               data = return_time_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  mod_return_time_list[[5]] = lmer(norm_return_time ~ as.factor(correlation)*as.factor(n_extinctions) +
                                 (1 | landscape),
                               data = return_time_data,
                               control = lmerControl(optimizer = "bobyqa"))
  
  
  
  
  
  mod_return_time_weights = mod_return_time_dredge$weight
  mod_return_time_coefficients = lapply(mod_return_time_list, coeffs)
  mod_return_time_vcov = lapply(mod_return_time_list, vcov)
  ## drawing parameters from the distribution of each model
  mod_return_time_coefficients_density = list()
  for (k in 1:5){
    mod_return_time_coefficients_density[[k]] = mvrnorm(n = as.integer(n_pred * mod_return_time_weights[k]),
                                                    mod_return_time_coefficients[[k]], mod_return_time_vcov[[k]])
  }
  
  predict_return_time = function(coeff, data, model){
    ext = data$n_extinctions == 8
    cor = data$correlation
    if (model == 1){
      return(coeff[1] + cor*coeff[2])
    }
    if (model == 3){
      return(coeff[1] + cor*coeff[2] + ext*coeff[3])
    }
    if (model == 5){
      return(coeff[1] + cor*coeff[2] + ext*coeff[3] + (ext*cor)*coeff[4])
    }
    if (model == 2){
      dummy = rep(1, nrow(data))
      return(dummy*coeff[1])
    }
    if (model == 4){
      return(coeff[1] + ext*coeff[2])
    }
  }
  
  
  
## Making the predictions for all variables ####################################################################################
  
  ## table for saving the predictions
  predictions_table = data.frame(n_extinctions = c(4, 4, 8, 8), correlation = c(F, T, F, T))
  predictions_table$norm_simpson = NA
  predictions_table$beta = NA
  predictions_table$norm_bioarea = NA
  predictions_table$norm_return_time = NA
  
  new_data = predictions_table
  
  for (k in 1:(n_pred-1)){
    predictions_table = rbind(predictions_table, new_data)
  }
  
  ## predicting alpha-diversity
  n = 1
  for (model in 1:5){
    coef_list = mod_simpson_coefficients_density[[model]]
    for (i in 1:nrow(coef_list)){
      simpson = predict_simpson(coef_list[i,], new_data, model)
      predictions_table$norm_simpson[n:(n+3)] = simpson
      n = n +4
    }
  }
  
  ## predicting beta-diversity
  n = 1
  for (model in 1){
    coef_list = mod_beta_coefficients_density[[model]]
    for (i in 1:nrow(coef_list)){
      beta = predict_beta(coef_list[i,], new_data, model)
      predictions_table$beta[n:(n+3)] = beta
      n = n +4
    }
  }
  
  ## predicting bioarea
  n = 1
  for (model in 1:5){
    coef_list = mod_bioarea_coefficients_density[[model]]
    for (i in 1:nrow(coef_list)){
      bioarea = predict_bioarea(coef_list[i,], new_data, model)
      predictions_table$norm_bioarea[n:(n+3)] = bioarea
      n = n +4
    }
  }
  
  ## predicting recovery time
  n = 1
  for (model in 1:5){
    coef_list = mod_return_time_coefficients_density[[model]]
    for (i in 1:nrow(coef_list)){
      return_time = predict_return_time(coef_list[i,], new_data, model)
      predictions_table$norm_return_time[n:(n+3)] = return_time
      n = n +4
    }
  }
  
  ## making labels for figure
  predictions_table$correlation_label = rep("", n_pred*4)
  predictions_table$correlation_label[!predictions_table$correlation] = "Dispersed extinctions"
  predictions_table$correlation_label[predictions_table$correlation] = "Clustered extinctions"
  
  
  ## Denormalizing the predictions
  predictions_table$simpson = predict(norm_simpson, predictions_table$norm_simpson, inverse = T)
  predictions_table$bioarea = predict(norm_bioarea, predictions_table$norm_bioarea, inverse = T)
  predictions_table$return_time = predict(norm_return_time, predictions_table$norm_return_time, inverse = T)
  
## summary of predictions ################################################################################
  # making a summary of the predictions with only the median and 95% CI for all values
  
  summary_table = data.frame(n_extinctions = c(4, 4, 8, 8), correlation = c(F, T, F, T))
  
  summary_table$simpson = NA
  summary_table$simpson_min = NA
  summary_table$simpson_max = NA
  
  summary_table$beta = NA
  summary_table$beta_min = NA
  summary_table$beta_max = NA
  
  summary_table$bioarea = NA
  summary_table$bioarea_min = NA
  summary_table$bioarea_max = NA
  
  summary_table$return_time = NA
  summary_table$return_time_min = NA
  summary_table$return_time_max = NA
  
  for (n in c(4, 8)){
    for (cor in c(T, F)){
      test_pred = (predictions_table$n_extinctions == n & predictions_table$correlation == cor)
      test_summary = (summary_table$n_extinctions == n & summary_table$correlation == cor)
      
      summary_table$simpson[test_summary] = median(predictions_table$simpson[test_summary], na.rm = T)
      summary_table$simpson_min[test_summary] = quantile(predictions_table$simpson[test_summary], 0.025, na.rm = T)
      summary_table$simpson_max[test_summary] = quantile(predictions_table$simpson[test_summary], 0.975, na.rm = T)
      
      summary_table$beta[test_summary] = median(predictions_table$beta[test_summary], na.rm = T)
      summary_table$beta_min[test_summary] = quantile(predictions_table$beta[test_summary], 0.025, na.rm = T)
      summary_table$beta_max[test_summary] = quantile(predictions_table$beta[test_summary], 0.975, na.rm = T)
      
      summary_table$bioarea[test_summary] = median(predictions_table$bioarea[test_summary], na.rm = T)
      summary_table$bioarea_min[test_summary] = quantile(predictions_table$bioarea[test_summary], 0.025, na.rm = T)
      summary_table$bioarea_max[test_summary] = quantile(predictions_table$bioarea[test_summary], 0.975, na.rm = T)
      
      summary_table$return_time[test_summary] = median(predictions_table$return_time[test_summary], na.rm = T)
      summary_table$return_time_min[test_summary] = quantile(predictions_table$return_time[test_summary], 0.025, na.rm = T)
      summary_table$return_time_max[test_summary] = quantile(predictions_table$return_time[test_summary], 0.975, na.rm = T)
    }
  }
  
  summary_table$correlation_label = rep("", 4)
  summary_table$correlation_label[!summary_table$correlation] = "Dispersed extinctions"
  summary_table$correlation_label[summary_table$correlation] = "Clustered extinctions"
  
  summary_table$simpson_center = (summary_table$simpson_max + summary_table$simpson_min)/2
  summary_table$simpson_height = summary_table$simpson_max - summary_table$simpson_min
  
  summary_table$beta_center = (summary_table$beta_max + summary_table$beta_min)/2
  summary_table$beta_height = summary_table$beta_max - summary_table$beta_min
  
  summary_table$bioarea_center = (log(summary_table$bioarea_max) + log(summary_table$bioarea_min))/2
  summary_table$bioarea_height = log(summary_table$bioarea_max) - log(summary_table$bioarea_min)
  
  summary_table$return_time_center = (summary_table$return_time_max + summary_table$return_time_min)/2
  summary_table$return_time_height = summary_table$return_time_max - summary_table$return_time_min
  
  
  
## plot predictions ##############################################################################################
  ## making Fig. 1
  
  
  ## setting the color scale 
  #colors = c("red", "blue")
  colors = c("red", "blue", "green")
  names(colors) = levels(post_ext_data$correlation_label)
  colScale = scale_colour_manual(name = "Treatment",values = colors)
  fillScale = scale_fill_manual(name = "Treatment",values = colors)
  
  
  ## preparing the plots
  #g = ggplot(data = post_ext_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label,
  #                                               color = correlation_label, group = correlation_label))
  
  #h = ggplot(data = post_ext_landscape_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label,
  #                                                         color = correlation_label, group = correlation_label))
  
  g = ggplot(data = fig_data, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label,
                                          color = correlation_label, group = correlation_label))

  h = ggplot(data = fig_data_landscape, mapping = aes(x = as.factor(n_extinctions), fill = correlation_label,
                                                    color = correlation_label, group = correlation_label))




  ## graphical parameters
  alpha_box = 0.3 ## boxes transparency
  alpha_points = 0.5 ## points transparency
  box_width = 0.4 ## width of boxes
  dodge_width = 0.75 ## horizontal offset
  
  panel_a = g +
    geom_jitter(aes(y = local_simpson), position=position_jitterdodge(jitter.width = 0.1), alpha = alpha_points) +
    colScale + 
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = simpson_center,
                                        width = box_width, height = simpson_height),
              position=position_dodge(width = dodge_width), alpha = alpha_box, color = NA) +
    fillScale +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = simpson,
                                        width = box_width, height = 0.01),
              position=position_dodge(width = dodge_width), alpha = 1, color = NA) +
    labs(y = expression(paste(alpha, "-diversity (Simpson's index)")),
         x = "Extinctions",
         fill = "Treatment",
         color = "Treatment") + guides(color=guide_legend(override.aes=list(fill=NA)))
  
  panel_b = h +
    geom_jitter(aes(y=BetaDiv), position = position_jitterdodge(jitter.width = 0.1), alpha = alpha_points) +
    colScale + 
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = beta_center,
                                        width = box_width, height = beta_height),
              position=position_dodge(width = dodge_width), alpha = alpha_box, color = NA) +
    fillScale +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = beta,
                                        width = box_width, height = 0.001),
              position=position_dodge(width = dodge_width), alpha = 1, color = NA) +
    labs(y = expression(paste(beta, "-diversity")),
         x = "Extinctions",
         fill = "Treatment",
         color = "Treatment") + guides(color=guide_legend(override.aes=list(fill=NA)))
  
  panel_c = g +
    geom_jitter(aes(y=log(bioarea)), position = position_jitterdodge(jitter.width = 0.1), alpha = alpha_points) +
    colScale + 
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = bioarea_center,
                                        width = box_width, height = bioarea_height),
              position=position_dodge(width = dodge_width), alpha = alpha_box, color = NA) +
    fillScale +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = log(bioarea),
                                        width = box_width, height = 0.05),
              position=position_dodge(width = dodge_width), alpha = 1, color = NA) +
    labs(y = expression(paste("Bioarea per volume (log(",mu,m^2,"/", mL, "))", sep = "")),
         x = "Extinctions",
         fill = "Treatment",
         color = "Treatment") + guides(color=guide_legend(override.aes=list(fill=NA)))
  
  panel_d = ggplot(data = tube_table[tube_table$extinctions == T,], mapping = aes(x = as.factor(n_extinctions), fill = correlation_label,
                                                                                  color = correlation_label, group = correlation_label)) +
    geom_jitter(aes(y=return_time), position = position_jitterdodge(jitter.width = 0.1), alpha = alpha_points) +
    colScale + 
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = return_time_center,
                                        width = box_width, height = return_time_height),
              position=position_dodge(width = dodge_width), alpha = alpha_box, color = NA) +
    fillScale +
    geom_tile(data = summary_table, aes(x = as.factor(n_extinctions), y = return_time,
                                        width = box_width, height = 1),
              position=position_dodge(width = dodge_width), alpha = 1, color = NA) +
    labs(y = "Recovery time (h)",
         x = "Extinctions",
         fill = "Treatment",
         color = "Treatment") + guides(color=guide_legend(override.aes=list(fill=NA)))
  
  cairo_pdf(file = "./fig1.pdf", 11.7, 8.3)
  ggarrange(panel_a, panel_b, panel_c, panel_d,
            labels = "auto",
            common.legend = TRUE, legend = "bottom")
  dev.off()