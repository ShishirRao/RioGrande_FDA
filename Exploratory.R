library(ggplot2)
library(dplyr)
library(tidyr)

setwd("F:/wild/PhD/ICON/Internship/RE__HEC-EFM_Results_for_FDA/")
trends<-read.table("Trends_flow_RGSM.05112020.csv", sep="," ,header=T)

#select only the response variables
fish_response = trends %>% select(1,2,3,4,6,8,9,11,12)

head(fish_response)
summary(fish_response)


#convert from wide to long format
fish_response = fish_response %>% gather(key = variable, value = "value",-num,-year,factor_key = TRUE)
fish_response$year = as.factor(fish_response$year)

fish_response_plot = ggplot(fish_response, aes(x = variable,y=value))+
  geom_boxplot()+facet_wrap(variable~.,scales = "free")+
  geom_jitter(aes(color = year),shape=16, position=position_jitter(0.2))+
  theme_bw()+ggtitle("Fish response")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave("fish_response_variability.jpg", fish_response_plot, device = "jpg",
       path = "F:/wild/PhD/ICON/Internship/RE__HEC-EFM_Results_for_FDA/",
       scale = 2, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)

##########################################################################

#select only the flow variables
flow_variables = trends %>% select(1,2,14:29)

#convert to long format
flow_variables = flow_variables %>% gather(key = variable, value = "value",-num,-year,factor_key = TRUE)
flow_variables$year = as.factor(flow_variables$year)

flow_variability = ggplot(flow_variables, aes(x = variable,y=value))+
  geom_boxplot()+facet_wrap(variable~.,scales = "free")+
  geom_jitter(aes(color = year),shape=16, position=position_jitter(0.2))+
  theme_bw()+ggtitle("Flow variability")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave("flow_variability.jpg", flow_variability, device = "jpg",
       path = "F:/wild/PhD/ICON/Internship/RE__HEC-EFM_Results_for_FDA/",
       scale = 2, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)


##########################################################################

#select only the EFM variables
flow_variables = trends %>% select(1,2,14:32)

#convert to long format
flow_variables = flow_variables %>% gather(key = variable, value = "value",-num,-year,factor_key = TRUE)
flow_variables$year = as.factor(flow_variables$year)

flow_variability = ggplot(flow_variables, aes(x = variable,y=value))+
  geom_boxplot()+facet_wrap(variable~.,scales = "free")+
  geom_jitter(aes(color = year),shape=16, position=position_jitter(0.2))+
  theme_bw()+ggtitle("Flow variability")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave("flow_variability.jpg", flow_variability, device = "jpg",
       path = "F:/wild/PhD/ICON/Internship/RE__HEC-EFM_Results_for_FDA/",
       scale = 2, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)


##########################################################################

#select only the EFM variables
EFM_variables = trends %>% select(1,2,33:74)

#convert to long format
EFM_variables = EFM_variables %>% gather(key = variable, value = "value",-num,-year,factor_key = TRUE)
EFM_variables$year = as.factor(EFM_variables$year)

EFM_variability = ggplot(EFM_variables, aes(x = variable,y=value))+
  geom_boxplot()+facet_wrap(variable~.,scales = "free")+
  geom_jitter(aes(color = year),shape=16, position=position_jitter(0.2))+
  theme_bw()+ggtitle("EFM variability")+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave("EFM_variability.jpg", EFM_variability, device = "jpg",
       path = "F:/wild/PhD/ICON/Internship/RE__HEC-EFM_Results_for_FDA/",
       scale = 2, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)
