# Santos et al. (in preparation)
# Assessing spatial patterns of sand flies
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(spatialEco)
library(readr)
library(dplyr)
#library(tidyr)
library(raster)
#library(rgdal)


# loading data ------------------------------------------------------------

spdata <- read_csv("./data/sandfly_data_paracambi_analysis.csv")

CRS = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# sample raster for kernel output
mapbiomas <- raster("./data/raster/mapbiomas_col5_paracambi_2003.tif")


# calculate relative sand fly frequencies by point ------------------------

# sampling effort
sampl <- spdata %>%
  filter(point_id != "NA") %>%
  group_by(point_id) %>%
  summarize(effort = n_distinct(month))

spp_freq <- spdata %>%
  filter(point_id != "NA",
         sp_code != "negative") %>%
  group_by(point_id) %>%
  summarize(n = sum(total)) %>%
  left_join(sampl, by = "point_id") %>%
  mutate(freq = round(n/effort, 2),
         n_log = log(n),
         freq_log = log(freq)) 


# convert into spatial points ---------------------------------------------

# converting spdata to spatial points
points <- spdata %>%
  filter(point_id != "NA") %>%
  dplyr::select(point_id, lon, lat) %>%
  distinct() %>%
  arrange(point_id) %>%
  left_join(spp_freq, by = "point_id") %>%
  drop_na()

points <- SpatialPointsDataFrame(coords = data.frame(points$lon, points$lat),
                                      data = points,
                                      proj4string = CRS)


# estimating kernel density -----------------------------------------------

# using frequencies (log)
kernel_map_freq <- sp.kde(points, y = points$freq_log, bw = 0.02, newdata = mapbiomas)

# using raw counts (log)
kernel_map_n <- sp.kde(points, y = points$n_log, bw = 0.02, newdata = mapbiomas)

# saving output
writeRaster(kernel_map_freq, "./outputs/04_kernel_freq_log.tif")
writeRaster(kernel_map_n, "./outputs/04_kernel_n_log.tif")
