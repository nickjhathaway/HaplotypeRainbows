
#' @title Prep a dataframe for clustering functions, this sorts the alleles by their within sample fraction rather than the population rank  
#' The input data is processed so that the information needed for the haplotype rainbow plotting functions
#'
#' @param inputData the input data 
#' @param sampleCol the name of the sample column 
#' @param targetCol the name of the target column
#' @param popUIDCol the name of the identifier column 
#' @param relAbundCol the name of the relative abundance column  
#' @param minPopSize the minimum population size of a target, if a target has less than unique alleles than the target is dropped 
#' @param colorOuput the number of colors to spread across 
#' @param barHeight the height of the final full bars per sample, controls whether they touch of not, e.g. barHeight==1 will create touching bars
#'
#' @returns the prepped data frame for the plotting function 
#' @export
#'
prepForRainbowArrangedByFrac <-function(inputData, sampleCol = s_Sample, targetCol= p_name, popUIDCol = h_popUID, relAbundCol = c_AveragedFrac, minPopSize = 1, colorOuput = 11, barHeight = 0.80){
  inputData_filt = inputData %>% 
    group_by({{sampleCol}})  %>% 
    mutate(targetNumber = length(unique({{targetCol}}))) %>% 
    group_by() %>% 
    mutate("{{sampleCol}}" := as.character({{sampleCol}})) %>% 
    group_by({{sampleCol}}) %>% 
    group_by() %>% 
    mutate("{{sampleCol}}" := factor({{sampleCol}})) %>% 
    group_by({{sampleCol}}) %>% 
    arrange({{popUIDCol}}) %>%
    group_by({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    summarise("{{relAbundCol}}" := sum({{relAbundCol}})) %>% 
    group_by({{sampleCol}}, {{targetCol}})
  inputData_filt = inputData_filt%>%
    group_by({{sampleCol}}, {{targetCol}}) %>% 
    mutate(totalAbund = sum({{relAbundCol}})) %>% 
    mutate("{{relAbundCol}}" :={{relAbundCol}}/totalAbund)%>% 
    arrange({{relAbundCol}}) 
  
  inputData_filt = inputData_filt %>% 
    group_by({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    mutate(s_COI = length(unique({{popUIDCol}})))
  
  
  inputData_filt = inputData_filt %>% 
    group_by() %>% 
    group_by({{sampleCol}}, {{targetCol}}) %>% 
    mutate(relAbundCol_mod = {{relAbundCol}} * barHeight) %>% 
    mutate(fracCumSum = cumsum({{relAbundCol}}) - {{relAbundCol}}) %>% 
    mutate(fracModCumSum = cumsum(relAbundCol_mod) - relAbundCol_mod) %>% 
    mutate(fakeFrac = 1/unique(s_COI))  %>% 
    mutate(fakeFracMod = fakeFrac * barHeight)  %>% 
    mutate(fakeFracCumSum = cumsum(fakeFrac) - fakeFrac) %>% 
    mutate(fakeFracModCumSum = cumsum(fakeFracMod) - fakeFracMod)
  
  inputData_filt_popName = inputData_filt %>% 
    select({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    unique() %>% 
    group_by({{targetCol}}, {{popUIDCol}}) %>% 
    summarise(samp_n = n()) %>% 
    arrange({{targetCol}}, desc(samp_n)) %>% 
    group_by({{targetCol}}) %>% 
    mutate(popid = row_number())  %>% 
    group_by({{targetCol}}) %>% 
    mutate(maxPopid = max(popid))
  
  
  
  inputData_filt = inputData_filt %>% 
    left_join(inputData_filt_popName)
  
  
  inputData_filt_tarFilt = inputData_filt %>% 
    filter(maxPopid >= minPopSize) %>% 
    group_by()
  
  inputData_filt_tarFilt = inputData_filt_tarFilt%>%
    group_by() %>% 
    mutate("{{targetCol}}" := factor({{targetCol}}))
  
  
  colorsOutput = colorOuput;
  #colorsOutput =  length(unique(inputData_filt_tarFilt${{targetCol}}))
  targetNumber = 0
  targetToHue = tibble()
  tempTarCol = inputData_filt_tarFilt %>% select({{targetCol}})
  
  for(tarname in levels(tempTarCol[[1]] ) ) {
    targetToHueForTarget = tibble("{{targetCol}}" := tarname, hueMod = (targetNumber%% colorsOutput) + 1)
    targetToHue = targetToHue %>% 
      bind_rows(targetToHueForTarget)
    targetNumber = targetNumber + 1;
  }
  inputData_filt_tarFilt = inputData_filt_tarFilt %>% 
    group_by({{targetCol}}) %>% 
    mutate(popidFrac = (popid-1)/(maxPopid))
  tempTarCol = inputData_filt_tarFilt %>% select({{targetCol}})
  
  targetToHue = targetToHue %>% 
    mutate("{{targetCol}}" := factor({{targetCol}}, levels = levels(tempTarCol[[1]])))
  inputData_filt_tarFilt = inputData_filt_tarFilt %>%
    left_join(targetToHue) %>%
    mutate(popidPerc = 100 * popidFrac) %>% 
    mutate(popidFracRegColor = round(abs((popidPerc + (hueMod/colorsOutput) *100) %% 200 -0.0001 ) %% 100) ) %>% 
    mutate(popidPercLog = log((popidFrac * 99) + 1 , base = 100) * 100 ) %>% 
    mutate(popidFracLogColor = round(abs((popidPercLog + (hueMod/colorsOutput) *100) %% 200 -0.0001 ) %% 100) ) 
  return(inputData_filt_tarFilt)
}


#' @title Prep a dataframe for clustering functions, this sorts the alleles by their population rank  
#' The input data is processed so that the information needed for the haplotype rainbow plotting functions
#'
#' @param inputData the input data 
#' @param sampleCol the name of the sample column 
#' @param targetCol the name of the target column
#' @param popUIDCol the name of the identifier column 
#' @param relAbundCol the name of the relative abundance column  
#' @param minPopSize the minimum population size of a target, if a target has less than unique alleles than the target is dropped 
#' @param colorOuput the number of colors to spread across 
#' @param barHeight the height of the final full bars per sample, controls whether they touch of not, e.g. barHeight==1 will create touching bars
#'
#' @returns the prepped data frame for the plotting function 
#' @export
prepForRainbow <-function(inputData, sampleCol = s_Sample, targetCol= p_name, popUIDCol = h_popUID, relAbundCol = c_AveragedFrac, minPopSize = 1, colorOuput = 11, barHeight = 0.80){
  inputData_filt = inputData %>% 
    group_by({{sampleCol}})  %>% 
    mutate(targetNumber = length(unique({{targetCol}}))) %>% 
    group_by() %>% 
    mutate("{{sampleCol}}" := as.character({{sampleCol}})) %>% 
    group_by({{sampleCol}}) %>% 
    group_by() %>% 
    mutate("{{sampleCol}}" := factor({{sampleCol}})) %>% 
    group_by({{sampleCol}}) %>% 
    arrange({{popUIDCol}}) %>% 
    group_by({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    summarise("{{relAbundCol}}" := sum({{relAbundCol}})) %>% 
    group_by({{sampleCol}}, {{targetCol}})
  inputData_filt = inputData_filt%>%
    group_by({{sampleCol}}, {{targetCol}}) %>% 
    mutate(totalAbund = sum({{relAbundCol}})) %>% 
    mutate("{{relAbundCol}}" :={{relAbundCol}}/totalAbund)
  
  inputData_filt = inputData_filt %>% 
    group_by({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    mutate(s_COI = length(unique({{popUIDCol}})))
  
  
  inputData_filt = inputData_filt %>% 
    group_by() %>% 
    group_by({{sampleCol}}, {{targetCol}}) %>% 
    mutate(relAbundCol_mod = {{relAbundCol}} * barHeight) %>% 
    mutate(fracCumSum = cumsum({{relAbundCol}}) - {{relAbundCol}}) %>% 
    mutate(fracModCumSum = cumsum(relAbundCol_mod) - relAbundCol_mod) %>% 
    mutate(fakeFrac = 1/unique(s_COI))  %>% 
    mutate(fakeFracMod = fakeFrac * barHeight)  %>% 
    mutate(fakeFracCumSum = cumsum(fakeFrac) - fakeFrac) %>% 
    mutate(fakeFracModCumSum = cumsum(fakeFracMod) - fakeFracMod)
  
  inputData_filt_popName = inputData_filt %>% 
    select({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    unique() %>% 
    group_by({{targetCol}}, {{popUIDCol}}) %>% 
    summarise(samp_n = n()) %>% 
    arrange({{targetCol}}, desc(samp_n)) %>% 
    group_by({{targetCol}}) %>% 
    mutate(popid = row_number())  %>% 
    group_by({{targetCol}}) %>% 
    mutate(maxPopid = max(popid))
  
  
  
  inputData_filt = inputData_filt %>% 
    left_join(inputData_filt_popName)
  
  
  inputData_filt_tarFilt = inputData_filt %>% 
    filter(maxPopid >= minPopSize) %>% 
    group_by()
  
  inputData_filt_tarFilt = inputData_filt_tarFilt%>%
    group_by() %>% 
    mutate("{{targetCol}}" := factor({{targetCol}}))
  
  
  colorsOutput = colorOuput;
  #colorsOutput =  length(unique(inputData_filt_tarFilt${{targetCol}}))
  targetNumber = 0
  targetToHue = tibble()
  tempTarCol = inputData_filt_tarFilt %>% select({{targetCol}})
  
  for(tarname in levels(tempTarCol[[1]] ) ) {
    targetToHueForTarget = tibble("{{targetCol}}" := tarname, hueMod = (targetNumber%% colorsOutput) + 1)
    targetToHue = targetToHue %>% 
      bind_rows(targetToHueForTarget)
    targetNumber = targetNumber + 1;
  }
  inputData_filt_tarFilt = inputData_filt_tarFilt %>% 
    group_by({{targetCol}}) %>% 
    mutate(popidFrac = (popid-1)/(maxPopid))
  tempTarCol = inputData_filt_tarFilt %>% select({{targetCol}})
  
  targetToHue = targetToHue %>% 
    mutate("{{targetCol}}" := factor({{targetCol}}, levels = levels(tempTarCol[[1]])))
  inputData_filt_tarFilt = inputData_filt_tarFilt %>%
    left_join(targetToHue) %>%
    mutate(popidPerc = 100 * popidFrac) %>% 
    mutate(popidFracRegColor = round(abs((popidPerc + (hueMod/colorsOutput) *100) %% 200 -0.0001 ) %% 100) ) %>% 
    mutate(popidPercLog = log((popidFrac * 99) + 1 , base = 100) * 100 ) %>% 
    mutate(popidFracLogColor = round(abs((popidPercLog + (hueMod/colorsOutput) *100) %% 200 -0.0001 ) %% 100) ) 
  return(inputData_filt_tarFilt)
}



#' @title Prep a dataframe for clustering functions, this preps for coloring by shading rather than by a rainbow of colors   
#' The input data is processed so that the information needed for the haplotype rainbow plotting functions
#'
#' @param inputData the input data 
#' @param sampleCol the name of the sample column 
#' @param targetCol the name of the target column
#' @param popUIDCol the name of the identifier column 
#' @param relAbundCol the name of the relative abundance column  
#' @param minPopSize the minimum population size of a target, if a target has less than unique alleles than the target is dropped 
#' @param baseColors the colors to create shades from 
#' @param barHeight the height of the final full bars per sample, controls whether they touch of not, e.g. barHeight==1 will create touching bars

#'
#' @returns the prepped data frame for the plotting function 
#' @export
prepForRainbowShade <-function(inputData, sampleCol = s_Sample, targetCol= p_name, popUIDCol = h_popUID, relAbundCol = c_AveragedFrac, minPopSize = 3, baseColors = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33'), barHeight = 0.80){
  inputData_filt = inputData %>% 
    group_by({{sampleCol}})  %>% 
    mutate(targetNumber = length(unique({{targetCol}}))) %>% 
    group_by() %>% 
    mutate("{{sampleCol}}" := as.character({{sampleCol}})) %>% 
    group_by({{sampleCol}}) %>% 
    group_by() %>% 
    mutate("{{sampleCol}}" := factor({{sampleCol}})) %>% 
    group_by({{sampleCol}}) %>% 
    arrange({{popUIDCol}}) %>% 
    group_by({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    summarise("{{relAbundCol}}" := sum({{relAbundCol}})) %>% 
    group_by({{sampleCol}}, {{targetCol}})
  inputData_filt = inputData_filt%>%
    group_by({{sampleCol}}, {{targetCol}}) %>% 
    mutate(totalAbund = sum({{relAbundCol}})) %>% 
    mutate("{{relAbundCol}}" :={{relAbundCol}}/totalAbund)
  
  inputData_filt = inputData_filt %>% 
    group_by({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    mutate(s_COI = length(unique({{popUIDCol}})))
  
  
  inputData_filt = inputData_filt %>% 
    group_by() %>% 
    group_by({{sampleCol}}, {{targetCol}}) %>% 
    mutate(relAbundCol_mod = {{relAbundCol}} * barHeight) %>% 
    mutate(fracCumSum = cumsum({{relAbundCol}}) - {{relAbundCol}}) %>% 
    mutate(fracModCumSum = cumsum(relAbundCol_mod) - relAbundCol_mod) %>% 
    mutate(fakeFrac = 1/unique(s_COI))  %>% 
    mutate(fakeFracMod = fakeFrac * barHeight)  %>% 
    mutate(fakeFracCumSum = cumsum(fakeFrac) - fakeFrac) %>% 
    mutate(fakeFracModCumSum = cumsum(fakeFracMod) - fakeFracMod)
  
  inputData_filt_popName = inputData_filt %>% 
    select({{sampleCol}}, {{targetCol}}, {{popUIDCol}}) %>% 
    unique() %>% 
    group_by({{targetCol}}, {{popUIDCol}}) %>% 
    summarise(samp_n = n()) %>% 
    arrange({{targetCol}}, desc(samp_n)) %>% 
    group_by({{targetCol}}) %>% 
    mutate(popid = row_number())  %>% 
    group_by({{targetCol}}) %>% 
    mutate(maxPopid = max(popid))
  
  
  
  inputData_filt = inputData_filt %>% 
    left_join(inputData_filt_popName)
  
  
  inputData_filt_tarFilt = inputData_filt %>% 
    filter(maxPopid >= minPopSize) %>% 
    group_by()
  
  inputData_filt_tarFilt = inputData_filt_tarFilt%>%
    group_by() %>% 
    mutate("{{targetCol}}" := factor({{targetCol}}))

  
  popCountsWithColors = inputData_filt %>% 
    group_by({{targetCol}}, {{popUIDCol}}) %>% 
    count() %>% 
    group_by({{targetCol}}) %>% 
    mutate(total = sum(n)) %>% 
    arrange({{targetCol}}, n) %>% 
    mutate(freq = n/total) %>% 
    group_by({{targetCol}}) %>% 
    mutate(p_uniqHaps = n()) %>% 
    mutate(h_id = row_number()) %>% 
    mutate(h_id_freq = h_id/p_uniqHaps)  %>% 
    mutate(h_id_freq_mod = h_id_freq * 0.75 + 0.25) %>% 
    mutate(cumFreq = cumsum(freq)) %>%
    mutate(modCumFreq = cumFreq * 0.75 + 0.25) %>%  
    mutate("{{targetCol}}" := factor({{targetCol}})) %>% 
    mutate(p_color_ID = (as.numeric({{targetCol}}) %% length(baseColors)) + 1 ) %>% 
    mutate(p_hue = (p_color_ID/length(baseColors)) ) %>% 
    mutate(p_baseColor = baseColors[p_color_ID]) %>% 
    mutate(h_color = hsv(p_hue, alpha = cumFreq)) %>% 
    mutate(h_color_mod = alpha(p_baseColor, alpha = modCumFreq)) %>% 
    mutate(h_color_byFreq = alpha(p_baseColor, alpha = h_id_freq))%>% 
    mutate(h_color_byFreq_mod = alpha(p_baseColor, alpha = h_id_freq_mod))
  

  inputData_filt_tarFilt = inputData_filt_tarFilt %>% 
    group_by({{targetCol}}) %>% 
    mutate(popidFrac = (popid-1)/(maxPopid))

  inputData_filt_tarFilt = inputData_filt_tarFilt %>%
    left_join(popCountsWithColors) 
  return(inputData_filt_tarFilt)
}


#' @title Create a ggplot object for the plotting by a rainblow of colors    
#'
#' @param prepData the prepped data, the columns used should be the same as the prepped 
#' @param sampleCol the name of the sample column 
#' @param targetCol the name of the target column
#' @param popUIDCol the name of the identifier column 
#' @param relAbundCol the name of the relative abundance column  
#' @param colorCol the color by which to color 
#' @param colors the colors to spread across, the number supplied should be the same number given to the prep functions 
#'
#' @returns a ggplot2 object that can be further modified 
#' @export
#'
genRainbowHapPlotObjActualFracLogColor <-function(prepData, sampleCol = s_Sample, targetCol= p_name, popUIDCol = h_popUID, relAbundCol = c_AveragedFrac, colorCol = popidFracLogColor, colors = RColorBrewer::brewer.pal(11, "Spectral")){
  sofonias_theme = theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
    theme(axis.line.x = element_line(color="black", size = 0.3),axis.line.y =
            element_line(color="black", size = 0.3))+
    theme(text=element_text(size=12, family="Helvetica"))+
    theme(axis.text.y = element_text(size=12))+
    theme(axis.text.x = element_text(size=12)) +
    theme(legend.position = "bottom") + 
    theme(plot.title = element_text(hjust = 0.5))
  sampleNamesDf = prepData %>% 
    group_by() %>% 
    select({{sampleCol}}) %>% 
    unique() %>% 
    arrange({{sampleCol}})
  min_sampleLevel = min(as.numeric(sampleNamesDf[[1]]))
  max_sampleLevel = max(as.numeric(sampleNamesDf[[1]]))
  allLevelsLabels = as.character(sampleNamesDf[[1]])
  names(allLevelsLabels) = as.numeric(sampleNamesDf[[1]])
  plotLabels = c()
  for(lev in min_sampleLevel:max_sampleLevel){
    if(lev %in% names(allLevelsLabels)){
      plotLabels = c(plotLabels, allLevelsLabels[as.character(lev)])
    }else{
      plotLabels = c(plotLabels, "")
    }
  }
    return (ggplot(prepData) + 
              geom_rect(aes(xmin = as.numeric({{targetCol}}) -0.5,
                            xmax = as.numeric({{targetCol}}) +0.5,
                            ymin = as.numeric({{sampleCol}}) + fracModCumSum - 0.5,
                            ymax = as.numeric({{sampleCol}}) + fracModCumSum + relAbundCol_mod - 0.5, 
                            fill = {{colorCol}},
                            "{{sampleCol}}" = {{sampleCol}},
                            "{{popUIDCol}}" = {{popUIDCol}},
                            "{{targetCol}}"= {{targetCol}},
                            "{{relAbundCol}}" = {{relAbundCol}}
              ), 
              color = "black") + 
              scale_fill_gradientn(colours = colors) + 
              scale_y_continuous(breaks = min_sampleLevel:max_sampleLevel, labels = plotLabels ) + 
              sofonias_theme + 
              theme(axis.text.x = element_blank()) + 
              guides(fill = "none"))
} 



#' @title Create a ggplot object for the plotting by a generating shades of input colors     
#'
#' @param prepData the prepped data, the columns used should be the same as the prepped 
#' @param sampleCol the name of the sample column 
#' @param targetCol the name of the target column
#' @param popUIDCol the name of the identifier column 
#' @param relAbundCol the name of the relative abundance column  
#' @param colorCol the color by which to color 
#'
#' @returns a ggplot2 object that can be further modified 
#' @export
#'
genRainbowHapPlotObjShade <-function(prepData, sampleCol = s_Sample, targetCol= p_name, popUIDCol = h_popUID, relAbundCol = c_AveragedFrac, colorCol = h_color_byFreq_mod){
  sofonias_theme = theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
    theme(axis.line.x = element_line(color="black", size = 0.3),axis.line.y =
            element_line(color="black", size = 0.3))+
    theme(text=element_text(size=12, family="Helvetica"))+
    theme(axis.text.y = element_text(size=12))+
    theme(axis.text.x = element_text(size=12)) +
    theme(legend.position = "bottom") + 
    theme(plot.title = element_text(hjust = 0.5))
  sampleNamesDf = prepData %>% 
    group_by() %>% 
    select({{sampleCol}}) %>% 
    unique() %>% 
    arrange({{sampleCol}})
  min_sampleLevel = min(as.numeric(sampleNamesDf[[1]]))
  max_sampleLevel = max(as.numeric(sampleNamesDf[[1]]))
  allLevelsLabels = as.character(sampleNamesDf[[1]])
  names(allLevelsLabels) = as.numeric(sampleNamesDf[[1]])
  plotLabels = c()
  for(lev in min_sampleLevel:max_sampleLevel){
    if(lev %in% names(allLevelsLabels)){
      plotLabels = c(plotLabels, allLevelsLabels[as.character(lev)])
    }else{
      plotLabels = c(plotLabels, "")
    }
  }
  return (ggplot(prepData) + 
            geom_rect(aes(xmin = as.numeric({{targetCol}}) -0.5,
                          xmax = as.numeric({{targetCol}}) +0.5,
                          ymin = as.numeric({{sampleCol}}) + fracModCumSum - 0.5,
                          ymax = as.numeric({{sampleCol}}) + fracModCumSum + relAbundCol_mod - 0.5, 
                          fill = {{colorCol}},
                          "{{sampleCol}}" = {{sampleCol}},
                          "{{popUIDCol}}" = {{popUIDCol}},
                          "{{targetCol}}"= {{targetCol}},
                          "{{relAbundCol}}" = {{relAbundCol}}
            ), 
            
            color = "black") + 
            
            scale_fill_identity() + 
            scale_y_continuous(breaks = min_sampleLevel:max_sampleLevel, labels = plotLabels ) + 
            sofonias_theme + 
            theme(axis.text.x = element_blank()) + 
            guides(fill = "none"))
} 

#' @title Create a ggplot object for the plotting by a rainblow of colors    
#'
#' @param prepData the prepped data, the columns used should be the same as the prepped 
#' @param sampleCol the name of the sample column 
#' @param targetCol the name of the target column
#' @param popUIDCol the name of the identifier column 
#' @param relAbundCol the name of the relative abundance column  
#' @param colorCol the color by which to color 
#' @param colors the colors to spread across, the number supplied should be the same number given to the prep functions 
#'
#' @returns a ggplot2 object that can be further modified 
#' @export
#'
genRainbowHapPlotObj <-function(prepData, sampleCol = s_Sample, targetCol= p_name, popUIDCol = h_popUID, relAbundCol = c_AveragedFrac, colorCol = popidFracLogColor, colors = RColorBrewer::brewer.pal(11, "Spectral")){
  genRainbowHapPlotObjActualFracLogColor(prepData, {{sampleCol}}, {{targetCol}}, {{popUIDCol}}, {{relAbundCol}}, {{colorCol}}, colors)
}


#' @title Add naming for the x-axis targets and reshape the output plot to have less empty space around plot 
#'
#' @param rainbow_plot The rainbow plot object created the gen rainbow plot functions 
#' @param prepped_allele_data the "preped" table created by the prep functions 
#' @param sampleCol the sample name column to use, this does not need to be the same as what was used for gen rainbow plot but can be a new name but must be a factor with same number of levels
#' @param targetCol the character name of the target column to use for the naming of the columns, this does not need to be the same as what was used for gen rainbow plot but can be a new name but must be a factor with same number of levels
#'
#' @returns the ggplot object 
#' @export
#'
add_haplotype_rainbow_axis_theming <-function(rainbow_plot, prepped_allele_data,
                                              sampleCol = s_Sample,
                                              targetCol = p_name){
  target_levels = levels(prepped_allele_data %>% pull({{targetCol}}))
  sample_levels = levels(prepped_allele_data %>% pull({{sampleCol}}))
  return(rainbow_plot + 
           scale_x_continuous(
             labels = target_levels, 
             breaks = 1:n_distinct(target_levels), 
             expand = c(0,0)
           ) + 
           scale_y_continuous(
             labels = sample_levels, 
             breaks = 1:n_distinct(sample_levels), 
             expand = c(0,0)
           ) + 
           theme(axis.text.y = element_text(family = "mono"), 
                 axis.text.x = element_text(family = "mono", angle = -90, hjust = 0)))
}



#' @title resort_prepped_samples_by_clustering 
#' Reset the sample factor levels to be in order that similar samples based on haplotype sharing
#'
#' @param prepped_allele_data the prepped data table 
#' @param sampleCol the name of the sample column 
#' @param targetCol the name of the target column
#' @param popUIDCol the name of the identifier column 
#' @param relAbundCol the name of the relative abundance column  
#' @param target_sample_coverage_freq_cut_off the sample coverage per target cut off, targets must have at least this freq of sample coverage to be used in the clustering
#' @param by_major_allele whether or not to cluster just by the major allele 
#'
#' @returns the prepped data frame with the resorted sample levels 
#' @export
#'
resort_prepped_samples_by_clustering <- function(
    prepped_allele_data,
    sampleCol = s_Sample,
    targetCol = p_name,
    popUIDCol = h_popUID,
    relAbundCol = c_AveragedFrac,
    target_sample_coverage_freq_cut_off = 0.80,
    by_major_allele = FALSE
) {
  
  stopifnot(is.data.frame(prepped_allele_data))
  
  # total distinct samples (used for coverage frequency)
  n_samples_total <- dplyr::n_distinct(dplyr::pull(prepped_allele_data, {{ sampleCol }}))
  
  # targets with sufficient sample coverage
  targets_keep <- prepped_allele_data %>%
    ungroup() %>% 
    group_by({{targetCol}}) %>% 
    dplyr::summarise(
      sample_count = dplyr::n_distinct({{ sampleCol }})
    ) %>%
    dplyr::mutate(sample_freq = sample_count / n_samples_total) %>%
    dplyr::filter(sample_freq >= target_sample_coverage_freq_cut_off) %>%
    dplyr::pull({{ targetCol }}) %>% 
    unique()
  if(length(targets_keep) == 0){
    stop("no targets above sample coverage cut off")
  }
  dat <- prepped_allele_data %>%
    filter({{targetCol}} %in% targets_keep)
  
  # optionally keep only major allele per (sample, target)
  if (by_major_allele) {
    dat <- dat %>%
      dplyr::group_by({{ sampleCol }}, {{ targetCol }}) %>%
      dplyr::slice_max(order_by = {{ relAbundCol }}, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
  }
  
  # build sample x popUID presence/absence matrix
  dat_sp <- dat %>%
    mutate(new_identifier = paste0({{targetCol}}, "-", {{popUIDCol}})) %>% 
    ungroup() %>% 
    dplyr::select({{ sampleCol }},  new_identifier) %>%
    unique() %>% 
    dplyr::mutate(marker = 1L) %>%
    tidyr::pivot_wider(
      names_from  = new_identifier,
      values_from = marker,
      values_fill = 0L
    )
  
  # matrix for clustering
  dat_sp_mat <- as.matrix(dat_sp[,2:ncol(dat_sp)])
  rownames(dat_sp_mat) <- dat_sp %>% 
    dplyr::pull({{ sampleCol }})
  dat_sp_mat_hc <- stats::hclust(stats::dist(dat_sp_mat), method = "ward.D2")
  sample_levels <- rownames(dat_sp_mat)[dat_sp_mat_hc$order]
  
  # add any samples that were filtered off to the end 
  any_missing_samples = prepped_allele_data %>%
    filter(!({{ sampleCol }} %in% sample_levels)) %>% 
    ungroup() %>% 
    dplyr::pull({{ sampleCol }})
  
  sample_levels = c(as.character(sample_levels), as.character(any_missing_samples))
  
  prepped_allele_data %>%
    dplyr::mutate({{ sampleCol }} := factor({{ sampleCol }}, levels = sample_levels))
}
