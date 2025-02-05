######### Conformal prediction framework as described in Siebers et al. 2025########
###### Author: Maud Siebers #######
###### Last modified: 5-02-2025 #######


#### example with Esthwaite Water UK with Polymer AC

###### packages #######
library(sf)
library(raster)
library(RPostgres)
library(dplyr)
library(lubridate)
library(stars)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(Metrics)
library(magrittr)
library(RColorBrewer)
library(ggplot2)

######## Read in file ########
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

#### read in MSI estimate file and add a year column
Chl_MSI<- read.csv('Chl_MSI_EsthwaiteWater.csv') %>% 
  mutate(Year = year(as.Date(Date_MSI, format = '%d/%m/%Y')))


###### Observed available data 
Observed <- read.csv('Observed_chl_EsthwaiteWater.csv') %>%
  rename(Observed = Chlorophyll) %>%
  mutate(
    Date = as.Date(Date, format = '%d/%m/%Y'),
    Observed = as.numeric(Observed),
    doy = yday(Date),
    year = year(Date)
  )


#### non-conformity score formula 
non_conformity_scores <- function(df_conform, df_test,df_Observed) {
  y_conf <- df_conform$MSI_chl
  y_conf_pred <- df_conform$Matchup
  alpha_conf <- abs(y_conf - y_conf_pred)
  
  alpha_test <- c()
  
  for (i in 1:nrow(df_test)) {
    y_test <- df_test$MSI_chl[i]
    y_test_pred <- df_test$Matchup[i]
    
    
    # Matchup and Non Matchup points in test data set are calculated differently
    if (df_test$Source[i] == 'Matchup') {
      alpha_test <- c(alpha_test, abs(y_test - y_test_pred))
    } else {
      obs_list <- df_Observed[c('Observed', 'doy', 'year')]
      differ <- abs(df_test$doy[i] - obs_list$doy)
      lowest_indices <- order(differ)[1:3]
      
      matching_multiplier <- 2
      values <- c()
      
      for (z_obs in lowest_indices) {
        if (df_Observed$year[z_obs] == df_test$Year[i]) {
          y_adjusted <- df_Observed$Observed[z_obs]
          values <- c(values, rep(y_adjusted, matching_multiplier))
        } else {
          values <- c(values, df_Observed$Observed[z_obs])
        }
      }
      
      y_adjusted <- mean(values)
      alpha_test <- c(alpha_test, abs(y_test - y_adjusted))
    }
  }
  
  return(list(alpha_conf = alpha_conf, alpha_test = alpha_test))
}



###################################################
### loop through conformal prediction ####


# Extract years
years_list <- c(unique(Chl_MSI$Year))
Total <- data.frame(iteration = NULL,Date = NULL, Observed = NULL,
                    
                    Source = NULL,MSI_chl = NULL,conformal = NULL,
                    upper = NULL, lower = NULL, uncertainty = NULL, 
                    uncertainty_diff = NULL)

test_year = 2017

for (test_year in years_list){
  ## test data = 1 year (test_year) + all no observation data

    test_data <- Chl_MSI %>%
      filter(Year == test_year | is.na(Date_observed)) %>%
      mutate(Source = ifelse(!is.na(Date_observed),'Matchup', 'No Matchup'))

  ### conformity data  - all years - test_year 
  conformity_data <- Chl_MSI %>%
    na.omit() %>% 
    filter(Year != test_year)
  
  # #### calculation of conformity scores
  test_data$doy <- yday(test_data$Date_MSI)   # day of year calculation
  conformity_data$doy <- yday(conformity_data$Date_observed) # day of year calculation

  results <- non_conformity_scores(conformity_data, test_data, Observed)
  
  
  final <- data.frame(Date = test_data$Date_MSI, Observed = test_data$Matchup, OWT =test_data$OWT, Source = test_data$Source,
                      iteration = test_year)
  # 
  
  conf_alpha <- conformity_data
  conf_alpha$alpha <- results$alpha_conf
  conf_alpha_sorted <- arrange(conf_alpha, alpha)
  epsilon <- 0.10
  qhat <- quantile(conf_alpha_sorted$alpha, 1 - epsilon, type =1)
  for (i in 1:length(conf_alpha_sorted$alpha)){
    if (conf_alpha_sorted$alpha[[i]] <= qhat){
      conf_alpha_sorted$conformal[[i]] <- 'Yes'
    }else{
      conf_alpha_sorted$conformal[[i]] <- 'No'
    }
  }
  
  test_alpha <- test_data
  test_alpha$alpha <- results$alpha_test
  
  for (i in 1:length(test_data$MSI_chl)){
    
    Y_test_pred <- test_data$MSI_chl[[i]]
    final$MSI_chl[[i]] <- Y_test_pred
    
    alpha_test <- results$alpha_test[[i]]
    final$alpha_test[[i]] <- alpha_test
    alpha_conf <- conf_alpha_sorted$alpha
    ### determine which alpha to use 
    pos <- findInterval(alpha_test, alpha_conf)
    final$position[[i]]<-  as.numeric(pos)
    final$quantile[[i]] <- as.numeric(pos/(length(conf_alpha_sorted$alpha)))
    final$conformal[[i]] <- ifelse(alpha_test <= qhat, 'Yes','No')
    # use quantiles from calibration estimates as baseline intervals
    p_value <- sum(alpha_conf >= alpha_test)/(length(alpha_conf)+1)
    final$p[[i]] <- as.numeric(p_value)
    
    epsilon <- 0.05
    
    lower_quantile <- quantile(alpha_conf, epsilon)
    
    upper_quantile <- quantile(alpha_conf, 1 - epsilon)
    
    
    median_abs_deviation <- median(alpha_conf) 
    
    scale_factor <- alpha_test / median_abs_deviation
    
    adjusted_lower_quantile <- lower_quantile * scale_factor
    
    adjusted_upper_quantile <- upper_quantile * scale_factor
    
    int2 <- (adjusted_upper_quantile - adjusted_lower_quantile)/2
    interval <- c(max(0,final$MSI_chl[[i]] - int2), final$MSI_chl[[i]] + int2)
    
    final$upper[[i]] <- as.numeric(interval[[2]])
    final$lower[[i]] <- as.numeric(interval[[1]])
  }
  
  final$uncertainty <- as.numeric(final$quantile) * 100
  final$conformal <- factor(as.character(final$conformal), levels = c("Yes", "No"))
  ### calculate the uncertainty
  final$width <- (as.numeric(final$upper) - as.numeric(final$lower))
  final$uncertainty_diff <- ((final$width/2)/as.numeric(final$MSI_chl))*100
  
  
  
  
  
  
  final <- final %>% dplyr::select(iteration,Date, Observed,Source,MSI_chl,conformal,upper, lower, uncertainty, uncertainty_diff,alpha_test)
  
  Total <- rbind(Total, final)
  
}


Total$conformal <- factor(as.character(Total$conformal), levels = c("Yes", "No"))

ggplot(data = Total, aes(x = as.Date(Date, format = '%d/%m/%Y'), y = as.numeric(MSI_chl))) +
  geom_point(aes(col = uncertainty, pch = conformal), size = 3) +
  geom_errorbar(
    aes(ymin = as.numeric(lower), ymax = as.numeric(upper)),
    width = 0.2  # Adjust width of the error bars as needed
    , col = 'red') +
  geom_line(data = Observed, aes(x = Date, y = Observed), col = 'black') +
  xlab("Date") +
  ylab("Predicted Chlorophyll a conc. (ug/L)") +
  xlim(as.Date('2016-01-01'), as.Date('2022-12-31')) +
  theme(panel.background = element_rect(fill = "darkgrey"))+
  facet_grid(iteration~conformal)+
  ylim(0,100)+
  scale_color_gradientn(colours = brewer.pal(9, 'OrRd'), limits = c(0, 100), na.value = "grey")+
  ggtitle('Conformal prediciton')+
  labs(color = "Conformity", pch = 'Accepted values')



conformal_counts <- Total %>%
  group_by(Date, Observed, Source, MSI_chl, conformal,alpha_test) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(
    names_from = conformal,
    values_from = count,
    values_fill = list(count = 0))  

conformal_counts <- conformal_counts %>% mutate(Conformal_percentage = (Yes/(Yes+No)*100))

Total_filter <- Total %>%
  group_by(Date, Observed,Source,MSI_chl, alpha_test) %>%
  summarise(
    upper = mean(as.numeric(upper), na.rm = TRUE),
    lower = mean(as.numeric(lower), na.rm = TRUE),
    uncertainty = mean(as.numeric(uncertainty), na.rm = TRUE),
    uncertainty_diff = mean(as.numeric(uncertainty_diff), na.rm = TRUE),
    .groups = 'drop')

Total_filter <- Total_filter %>%
  left_join(conformal_counts[c("Date", "Observed", "Source", "MSI_chl",'Conformal_percentage')], 
            by = c("Date", "Observed", "Source", "MSI_chl")) 

Total_filter<- Total_filter%>% 
  mutate(conformal = ifelse(Conformal_percentage > 90, 'Yes', 'No'))




######## plotting of points with their uncertainty metrics 

Total_filter$conformal_2 <- ifelse(Total_filter$conformal == 'Yes', 'Accepted values', 'Rejected values')
ggplot(data = Total_filter, aes(x = as.Date(Date, format = '%d/%m/%Y'), y = as.numeric(MSI_chl))) +
  # Observed line and points with linetype mapping for the legend
  geom_line(data = Observed, aes(x = Date, y = Observed, linetype = "In situ chl-a"), color = "black") +
  geom_point(data = Observed, aes(x = Date, y = Observed), shape = 16, color = "black", show.legend = TRUE) +  # Use shape for In situ chl-a
  
  # Points for uncertainty with color gradient
  geom_point(aes(color = uncertainty, shape = conformal_2), size = 3) +
  
  # Error bars
  geom_errorbar(
    aes(ymin = as.numeric(lower), ymax = as.numeric(upper)),
    width = 0.2, col = 'red') +
  
  # Labels and limits
  xlab("Date") +
  xlim(as.Date('2016-01-01'), as.Date('2022-12-31')) +
  ylab(expression('MSI chl-a estimate (mg m'^-3*')')) +
  coord_cartesian(ylim = c(0, 150)) +
  
  
  # Color gradient for uncertainty
  scale_color_gradientn(
    name = 'Non-conformity (%)', 
    colours = c("yellow", "orange", 'red'), 
    limits = c(0, 100), 
    na.value = "grey"
  ) +
  
  # Customize legend for 'In situ chl-a'
  scale_linetype_manual(name = element_blank(), values = c("In situ chl-a" = "solid")) +  # Line for 'In situ chl-a'
  scale_shape_manual(name = element_blank(), 
                     values = c(  # Point for 'In situ chl-a'
                       'Accepted values' = 16, 
                       'Rejected values' = 8)) + 
  ggtitle('MSI chl-a with uncertainty %')+
  # Theme
  theme_bw() 


########### Write a CSV file with the conformal outputs 
# writing_file <- Total_filter
# writing_file$Prediction <- as.numeric(writing_file$MSI_chl)
# writing_file$alpha_test <- as.numeric(writing_file$alpha_test)
# 
# write_csv(writing_file, paste0('Conformal_output_iterative.csv'))



