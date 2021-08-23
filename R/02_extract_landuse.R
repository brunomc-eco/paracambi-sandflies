# Santos et al. (in preparation)
# Extracting land use/cover data from mapbiomas in each sand fly sampling point
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(raster)
library(rgdal)
library(data.table)


# loading data ------------------------------------------------------------

# point data
spdata <- read_csv("./data/sandfly_data_paracambi_analysis.csv")

# mapbiomas rasters
mapbiomas <- list.files("./data/raster", pattern = "paracambi", full.names = TRUE) %>%
  stack()

# names of mapbiomas categories present in Paracambi
cod <- data.frame(cod = c(3, 15, 21, 25, 24),
                  name = c("forest", "pasture", "agripasture", "nonveg", "urban"))

# buffer size for data extraction, in meters
buf = 200

# wgs84 datum
CRS <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 

#mapbiomas_2002 <- raster("./data/raster/mapbiomas-brazil-collection-50-riodejaneiro-2002.tif")


# extracting land use values per point/year -------------------------------

# converting spdata to spatial points
points <- spdata %>%
  dplyr::select(point_id, ano, lon_campo, lat_campo) %>%
  distinct() %>%
  arrange(point_id) %>%
  drop_na() 

year <- sort(unique(points$ano))

landcov_vals <- list()
for(i in 1:length(year)){
  
  # select only points in this year
  points_year <- filter(points, ano == year[[i]])
  
  # convert to spatial points
  points_year <- SpatialPointsDataFrame(coords = data.frame(points_year$lon_campo, 
                                                            points_year$lat_campo),
                                   data = points_year,
                                   proj4string = CRS)
  
  # extract pixel values in buffers
  vals <- raster::extract(mapbiomas[[i]], points_year, buffer=buf)
  
  # calculate frequency of each category in buffers
  freq <- lapply(vals, table)
  
  tables <- list()
  for(j in 1:length(freq)){
    tables[[j]] <- data.frame(id = points_year$point_id[[j]],
                              unlist(freq[[j]]))
  }
  
  landcov_vals[[i]] <- rbindlist(tables) %>%
    #mutate(area_m2 = Freq*30) %>% # calculate area based on pixel size of 30m
    mutate(year = year[[i]])
  
}


# generate output table
landcov <- rbindlist(landcov_vals) %>%
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  rename(forest = "3",
         pasture = "15",
         agripasture = "21",
         nonveg = "25",
         urban = "24") %>%
  relocate(year, .after = id) %>%
  mutate(point_year_id = str_c(id, "_", year)) %>%
  relocate(point_year_id, .after = year) %>%
  replace_na(list(forest = 0,
                  pasture = 0,
                  agripasture = 0,
                  nonveg = 0,
                  urban = 0)) %>%
  mutate(sum = forest+pasture+agripasture+nonveg+urban,
         forest_f = forest/sum,
         pasture_f = pasture/sum,
         agripasture_f = agripasture/sum,
         nonveg_f = nonveg/sum,
         urban_f = urban/sum)


# save output table
write_csv(landcov, "./outputs/02_land_covers.csv")
