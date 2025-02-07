######### Conformal prediction framework as described in Siebers et al. 2025########
###### Author: Maud Siebers #######
###### Last modified: 5-02-2025 #######


#### example with Estwaithe Water UK with Polymer AC

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
Chl_MSI<- read.csv('Chl_MSI_ExampleLake.csv') %>%   ### amend to your own MSI estimate files
  mutate(Year = year(as.Date(Date, format = '%d/%m/%Y')),
         Date = as.Date(Date, format = '%d/%m/%Y'))


###### Observed available data 
Observed <- read.csv('Observed_chl.csv') %>%  #### amend to your own in situ chl
  rename(Observed = Chlorophyll) %>%
  mutate(
    Date = as.Date(Date, format = '%d/%m/%Y'),
    Observed = as.numeric(Observed),
    doy = yday(Date),
    Year = year(Date)
  )


#### non-conformity score formula 
non_conformity_scores <- function(df_conform, df_test,df_Observed) {
  y_conf <- df_conform$MSI_chl 
  y_conf_pred <- df_conform$Observed
  alpha_conf <- abs(y_conf - y_conf_pred)
  
  alpha_test <- c()
  
  for (i in 1:nrow(df_test)) {
    y_test <- df_test$MSI_chl[i]
    y_test_pred <- df_test$Matchup[i]
    
    
    # Matchup and Non Matchup points in test data set are calculated differently
    if (df_test$Source[i] == 'Matchup') {
      alpha_test <- c(alpha_test, abs(y_test - y_test_pred))
    } else {
      obs_list <- df_Observed[c('Observed', 'doy', 'Year')]
      differ <- abs(df_test$doy[i] - obs_list$doy)
      lowest_indices <- order(differ)[1:3]
      
      matching_multiplier <- 2
      values <- c()
      
      for (z_obs in lowest_indices) {
        if (df_Observed$Year[z_obs] == df_test$Year[i]) {
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
### Loop through different years of test data ###


# Extract unique years in the dataset to loop through
years_list <- c(unique(Chl_MSI$Year))  

# Set up final dataframe 
Conformal_output <- data.frame(iteration = NULL,Date = NULL, Observed = NULL,
                    Source = NULL,MSI_chl = NULL,conformal = NULL,
                    upper = NULL, lower = NULL, uncertainty = NULL)


### iteration loop####


for (test_year in years_list){
  ## test data = 1 year (test_year) + all no observation data
    Chl_MSI_merged <- merge(Observed, Chl_MSI, by = c("Date", 'Year'), all = TRUE)
    test_data <- Chl_MSI_merged %>%
      filter(Year == test_year | is.na(Observed)) %>%
      mutate(Source = ifelse(!is.na(Observed),'Matchup', 'No Matchup'))

  ### conformity data  - all years - test_year 
  conformity_data <- Chl_MSI_merged %>%
    na.omit() %>% 
    filter(Year != test_year)
  
  # #### calculation of conformity scores
  test_data$doy <- yday(test_data$Date)   # day of year calculation
  conformity_data$doy <- yday(conformity_data$Date) # day of year calculation

  
  results <- non_conformity_scores(conformity_data, test_data, Observed) # output of non-conformity score calculation
  
  ### Temporary file to keep data stored throughout loop
  temp_file <- data.frame(Date = test_data$Date, Observed = test_data$Observed, OWT =test_data$OWT, Source = test_data$Source,
                      iteration = test_year)
  
  ## extraction of conformity scores from 'results' file and combining with conformtity file
  conf_alpha <- conformity_data
  conf_alpha$alpha <- results$alpha_conf
  
  # sort non-conformity scores (conformity dataset)
  conf_alpha_sorted <- arrange(conf_alpha, alpha)
  epsilon <- 0.05 #set threshold 
  
  #determine which points in conformity set fall outside of confidence interval
  qhat <- quantile(conf_alpha_sorted$alpha, 1 - epsilon, type =1)
  for (i in 1:length(conf_alpha_sorted$alpha)){
    if (conf_alpha_sorted$alpha[[i]] <= qhat){
      conf_alpha_sorted$conformal[[i]] <- 'Yes'
    }else{
      conf_alpha_sorted$conformal[[i]] <- 'No'
    }
  }
 
  ## extraction of test scores from 'results' file and combining with test file
  test_alpha <- test_data
  test_alpha$alpha <- results$alpha_test
  
  #relating the alpha conformity to the alpha test
  
  for (i in 1:length(test_data$MSI_chl)){
    
    Y_test_pred <- test_data$MSI_chl[[i]]
    temp_file$MSI_chl[[i]] <- Y_test_pred
    
    alpha_test <- results$alpha_test[[i]]
    temp_file$alpha_test[[i]] <- alpha_test
    alpha_conf <- conf_alpha_sorted$alpha
    ### determine the position of the alpha_test in sorted alpha_confs
    pos <- findInterval(alpha_test, alpha_conf)
    temp_file$position[[i]]<-  as.numeric(pos)
    #calculate its exact quatile based on the position
    temp_file$quantile[[i]] <- as.numeric(pos/(length(conf_alpha_sorted$alpha)))
    #determin if this is conformal or not 
    temp_file$conformal[[i]] <- ifelse(alpha_test <= qhat, 'Yes','No')
    # use quantiles from calibration estimates as baseline intervals
    p_value <- sum(alpha_conf >= alpha_test)/(length(alpha_conf)+1)
    temp_file$p[[i]] <- as.numeric(p_value)
    
    #upper and lower quantile 
    lower_quantile <- quantile(alpha_conf, epsilon)
    
    upper_quantile <- quantile(alpha_conf, 1 - epsilon)
    
    ### scale the prediction intervals
    median_abs_deviation <- median(alpha_conf) 
    scale_factor <- alpha_test / median_abs_deviation
    adjusted_lower_quantile <- lower_quantile * scale_factor
    adjusted_upper_quantile <- upper_quantile * scale_factor
    
    #calculate intervals
    int2 <- (adjusted_upper_quantile - adjusted_lower_quantile)/2
    interval <- c(max(0,temp_file$MSI_chl[[i]] - int2), temp_file$MSI_chl[[i]] + int2)
    
    temp_file$upper[[i]] <- as.numeric(interval[[2]])
    temp_file$lower[[i]] <- as.numeric(interval[[1]])
  }
  
  #non-conformtiy scores to uncertainty percentages
  temp_file$uncertainty <- as.numeric(temp_file$quantile) * 100
  
  #For easy of reading put a Yes No classification on 
  temp_file$conformal <- factor(as.character(temp_file$conformal), levels = c("Yes", "No"))

  
  
  temp_file <- temp_file %>% dplyr::select(iteration,Date, Observed,Source,MSI_chl,conformal,upper, lower, uncertainty, alpha_test)
  
  Conformal_output <- rbind(Conformal_output, temp_file)
  
}


### plot the different iterations seperately 
ggplot(data = Conformal_output, aes(x = as.Date(Date, format = '%d/%m/%Y'), y = as.numeric(MSI_chl))) +
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
  ggtitle('Conformal prediciton per iteration')+
  labs(color = "Conformity", pch = 'Accepted values')

######## plotting of points with their uncertainty metrics 

Conformal_output$conformal_2 <- ifelse(Conformal_ouput$conformal == 'Yes', 'Accepted values', 'Rejected values')
ggplot(data = Conformal_ouput, aes(x = as.Date(Date, format = '%d/%m/%Y'), y = as.numeric(MSI_chl))) +
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
writing_file <- transform(Conformal_output, 
                          MSI_chl = as.numeric(MSI_chl), 
                          alpha_test = as.numeric(alpha_test), 
                          upper = as.numeric(upper), 
                        lower = as.numeric(lower))


write_csv(writing_file, paste0('Conformal_output_iterative.csv'))



