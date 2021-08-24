# Santos et al. (in preparation)
# Extracting land use/cover data from mapbiomas in each sand fly sampling point
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(tidyr)
library(raster)
library(rgdal)
library(data.table)


# loading data ------------------------------------------------------------

# point data
spdata <- read_csv("./data/sandfly_data_paracambi_analysis.csv")

# mapbiomas rasters
mapbiomas <- list.files("./data/raster", pattern = "paracambi", full.names = TRUE) %>%
  stack()

# names of mapbiomas categories in the study area
cod <- data.frame(cod = c(3, 15, 21, 25, 24, 33),
                  name = c("forest", "pasture", "agripasture", "nonveg", "urban", "water"))

# buffer size for data extraction, in meters
buf = 200

# wgs84 datum
CRS <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 


# extracting land use values per point/year -------------------------------

# getting coordinate data
points <- spdata %>%
  dplyr::select(point_id, year, point_year_id, lon, lat) %>%
  distinct() %>%
  arrange(point_id) %>%
  drop_na()

years <- sort(unique(points$year))

landcov_vals <- list()
for(i in 1:length(years)){
  
  # select only points in this year
  points_year <- filter(points, year == years[[i]])
  
  # convert to spatial points
  points_year <- SpatialPointsDataFrame(coords = data.frame(points_year$lon, 
                                                            points_year$lat),
                                   data = points_year,
                                   proj4string = CRS)
  
  # extract pixel values in buffers
  vals <- raster::extract(mapbiomas[[i]], points_year, buffer=buf)
  
  # calculate frequency of each category in buffers
  freq <- lapply(vals, table)
  
  tables <- list()
  for(j in 1:length(freq)){
    tables[[j]] <- data.frame(id = points_year$point_year_id[[j]],
                              unlist(freq[[j]])) %>%
      mutate(freq_percent = Freq/sum(Freq))
  }
  
  landcov_vals[[i]] <- rbindlist(tables)
  
}


# generate output table
landcov <- rbindlist(landcov_vals) %>%
  dplyr::select(!Freq) %>%
  pivot_wider(names_from = Var1, values_from = freq_percent) 

if(ncol(landcov) == 7){
  
  landcov <- landcov %>%
    rename(point_year_id = "id",
           forest = "3",
           pasture = "15",
           agripasture = "21",
           nonveg = "25",
           urban = "24",
           water = "33") %>%
    replace_na(list(forest = 0,
                    pasture = 0,
                    agripasture = 0,
                    nonveg = 0,
                    urban = 0,
                    water = 0))
  
} else {
  
  landcov <- landcov %>%
    rename(point_year_id = "id",
           forest = "3",
           pasture = "15",
           agripasture = "21",
           nonveg = "25",
           urban = "24") %>%
    replace_na(list(forest = 0,
                    pasture = 0,
                    agripasture = 0,
                    nonveg = 0,
                    urban = 0))
}

# save output table
write_csv(landcov, paste0("./outputs/02_land_covers_", buf, "m_buffer.csv"))
