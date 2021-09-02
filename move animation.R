##### Path Animation #####

library(anipaths)
library(ggmap)
library(tidyverse)

source("./Sera_Functions.R")

# Load dataset and convert to latlong
output <- read.csv('./movdata/Sera_datasets_10Mar21_HMMclassified_20210427.csv') 

# remove missing relocs
df <- output %>%
  drop_na(x,y) %>%
  droplevels()

output.latlong <- AddLatLong(df)

# get background map
# https://console.cloud.google.com/google/maps-apis/overview
ggmap::register_google(key = "AIzaSyAyZnyz0E5teo9SbaKvyvoxkV1sAz-ty10", write = TRUE)
background <- ggmap::get_googlemap(center = c(37.800, 1.097),
                                   zoom = 12,
                                   maptype = "satellite")


##### Animate 
# create animation dataset - cohort
ele <- filter(output.latlong, cohort %in% c(1,2,3)) %>% # filter here
  mutate(date = as.POSIXct(date, tz = "Africa/Nairobi")) %>%
  filter(as.Date(date) >= as.Date("2019-05-02")) %>%
  mutate(cohort = as.factor(cohort))


# color pallets sourced at https://colorbrewer2.org/
color.pal <- c('#a50f14', '#de2d26', '#fb6a4a', #red
               '#005a32', '#238b45', #green  
               '#08519c', '#3182bd', '#6baed6') #blue

# animation function
animate_paths(paths = ele,
              delta.t = "day", # daily movement
              coord = c("location.long", "location.lat"),
              Time.name = "date",
              ID.name = "MovDataID",
              pt.cex = 2,
              pt.colors = color.pal,
              background = background,
              tail.wd = .5,
              tail.length = 1,
              method = "mp4")


##### Animate all cohorts + residents
# create animation dataset - cohort
ele <- output.latlong %>%
  mutate(date = as.POSIXct(date, tz = "Africa/Nairobi")) %>%
  filter(as.Date(date) >= as.Date("2019-05-02")) %>%
  mutate(cohort = as.factor(cohort))


color.pal <- c('#a50f14', '#de2d26', '#fb6a4a', #red
               '#005a32', '#238b45', #green  
               '#08519c', '#3182bd', '#6baed6') #blue
               '#000000', '#252525', '#525252', '#737373', '#969696') # black
              
# animation function
animate_paths(paths = ele,
              delta.t = "day", # daily movement
              coord = c("location.long", "location.lat"),
              Time.name = "date",
              ID.name = "MovDataID",
              pt.cex = 2,
              pt.colors = color.pal,
              background = background,
              tail.wd = .5,
              tail.length = 1,
              method = "mp4")


