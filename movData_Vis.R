## HMM Data Vis Code
library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(mapview)
library(ggmap)
library(anipaths)

# DESCRIPTION
# 1. Build trajectories and plot post-release data
# 2. Animation code
# 3. Export shapefiles of trajectories (polylines) and GPS relocations (points)

# load dataset 
output <- read.csv('./movdata/Sera_datasets_10Mar21_HMMclassified_20210524.csv') 

# remove missing relocs
df <- output %>%
  drop_na(x,y) %>%
  droplevels()

##### 1. Plot Trajectories #####
# Have to split dataframe and build trajectories for each individual
# Final output converted to an sf object and plotted in mapview

# split dataframe and create ID column
df.split <- split(df, df$id)
df.split <- lapply(df.split, function(x) x %>% mutate(row.id = seq.int(nrow(.))))

# prep spatial data as a list
sf.split <- lapply(df.split, st_as_sf, coords = c("x", "y"), crs = 32637)
coords <- lapply(sf.split, st_coordinates)
begin.coord <- lapply(coords, function(x) as.matrix(x[1:nrow(x)-1,]))
end.coord <- lapply(coords, function(x) as.matrix(x[2:nrow(x),]))

# for loop creates lines in inner loop, and builds the trajectories in outer loop
sl.list <- vector("list", length(df.split))
sldf.list <- vector("list", length(df.split))

system.time({
  
for(i in 1:length(df.split)){
  l <- vector("list", nrow(begin.coord[[i]])) # adaptive list for lines 
  #create lines
  for (j in seq_along(l)) {
    l[[j]] <- Lines(list(Line(rbind(begin.coord[[i]][j, ], end.coord[[i]][j,]))), as.character(j))
  }
  # create trjaectories
  sl.list[[i]] <- SpatialLines(l, proj4string = CRS('+proj=utm +init=epsg:32637'))
  # remove the last row of the dataframe because it doesn't have a polyline
  sldf.list[[i]] <- SpatialLinesDataFrame(sl.list[[i]], data = head(df.split[[i]], -1), match.ID = FALSE)
}
  
})

# Create sf object
traj.sf <- lapply(sldf.list, st_as_sf) # convert to sf 
names(traj.sf) <- names(df.split)

# filter to first x months post-release 
traj.sf.filter <- lapply(traj.sf, filter, as.Date(date) >= as.Date(releaseDate) + 1 & 
                           # adjust time post release (n months * 30 days)
                           as.Date(date) <= (as.Date(releaseDate) + 6*30))

## Preview using mapview
library(mapview)
waterpoints <- st_read("./spatial data/Sera_PermanentWaterSources_May2021.shp")

# Lingwezi
mapview(traj.sf.filter$`1366`["viterbi"]) + waterpoints

# Pokot
mapview(traj.sf.filter$`1387`["viterbi"]) + waterpoints

# Nchurai
mapview(traj.sf.filter$`1520`["viterbi"]) + waterpoints


##### 3. Write to shapefiles #####

# set wd before exporting! 
setwd('./movdata/traj output/trajectories')


## write filtered trajectories to disk

# by individual
sapply(names(traj.sf.filter), 
       function (x) st_write(traj.sf.filter[[x]], paste(x, "traj_6month_2021-05-20", sep="_"), driver = "ESRI Shapefile"))

## write GPS relocations to disk
# set wd before exporting!

# create locs
df %>% group_by()
locs.sf <- st_as_sf(df, coords = c("x", "y"), crs = 32637) # for some reason doesn't work with ..36
locs.sf <- split(locs.sf, locs.sf$id)

# filter to 6 months post release
locs.sf.filter <- lapply(locs.sf, filter, as.Date(date) >= as.Date(releaseDate) + 1 & 
                           # adjust time post release (n months * 30 days)
                           as.Date(date) <= (as.Date(releaseDate) + 6*30))


# # write to disk
# set wd
setwd('..')
setwd('./relocations')

# full
#st_write(ele.sf, "Sera_datasets_10Mar21_HMMclassified_20210427", driver = "ESRI Shapefile")
# by individual
sapply(names(locs.sf.filter),
       function (x) st_write(locs.sf.filter[[x]], paste(x, "locs_6month_2021-05-20", sep="_"), driver = "ESRI Shapefile"))
# rds
saveRDS(ele.sf, "EleCollars_211119_locs_2019-11-30.rds")




