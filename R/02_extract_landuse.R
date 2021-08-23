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
mapbiomas <- raster("./data/raster/mapbiomas-brazil-collection-50-riodejaneiro-2002.tif")
CRS = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # wgs84 datum

# Paracambi shapefile
paracambi <- readOGR("./data/shp/Paracambi_2020.shp")

# crop mapbiomas to Paracambi extent
mapbiomas_paracambi <- crop(mapbiomas, paracambi)

# names of mapbiomas categories
cod <- data.frame(cod = c(3, 15, 21, 25, 24),
                  name = c("forest", "pasture", "agripasture", "nonveg", "urban"))

# buffer size for data extraction, in meters
buf = 200


# extracting land use values per point/year -------------------------------

points <- spdata %>%
  dplyr::select(point_id, ano, lon_campo, lat_campo) %>%
  filter(ano == 2002) %>%
  distinct() %>%
  arrange(point_id) %>%
  drop_na()

# converting to spatial points
points <- SpatialPointsDataFrame(coords = data.frame(points$lon_campo, points$lat_campo),
                            data = points,
                            proj4string = CRS)

# extract land cover percentages 
cobvals <- raster::extract(mapbiomas_paracambi, points, buffer=buf)

#cobvals <- list()
#for(i in 1:length(points)){
#  cobvals[[i]] <- raster::extract(mapbiomas_paracambi, 
#                                  points[[i]],
#                                  buffer=buffer) #extrai os valores dos pixels dentro de cada buffer
#}

#freq <- list()
#for(i in 1:length(cobvals)){
#  freq[[i]] <- lapply(cobvals[[i]], table) #calcula a frequencia de cada tipo de cobertura dentro do #buffer
#}

# get frequencies of mapbiomas categories in each buffer
freq <- lapply(cobvals, table)

# generate tables
tables <- list()
for(i in 1:length(freq)){
  tables[[i]] <- data.frame(id = points$point_id[[i]],
                            unlist(freq[[i]]))
}

pixels <- rbindlist(tables) %>%
  mutate(area_m2 = Freq*30) # calculating area based on pixel size of 30m

landcov <- pixels %>%
  pivot_wider(names_from = Var1, values_from = Freq) %>%
  rename(point_year_id = id,
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

# save output table
write_csv(landcov, "./outputs/02_land_covers.csv")
