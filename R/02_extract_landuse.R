# Extract land use data from mapbiomas in each sand fly sampling point

library(readr)
library(dplyr)
library(tidyr)
library(raster)
library(rgdal)
library(data.table)

#load point data
dados <- read_csv("./data/flebs_paracambi_ginelza_completo_check.csv")

#load rasters
CRS = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
mapbiomas <- raster("./data/raster/mapbiomas-brazil-collection-50-riodejaneiro-2002.tif")
paracambi <- readOGR("./data/shp/paracambi.shp")
mapbiomas_paracambi <- crop(mapbiomas, paracambi)

cod <- data.frame(cod = c(3, 15, 21, 25, 24),
                  name = c("forest", "pasture", "agripasture", "nonveg", "urban"))

#### Part 1: Preparing point data 2002 ####

points <- dados %>%
  dplyr::select(point_id, ano, lon_campo, lat_campo) %>%
  filter(ano == 2002) %>%
  distinct() %>%
  arrange(point_id) %>%
  drop_na()

points <- SpatialPointsDataFrame(coords = data.frame(points$lon_campo, points$lat_campo),
                            data = points,
                            proj4string = CRS)

#### Part 2: extract land use values per point/year ####

#buffer size in meters
buffer = 200

#extract land cover percentages 
cobvals <- raster::extract(mapbiomas_paracambi, points, buffer=buffer)

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

freq <- lapply(cobvals, table)

#tables
tables <- list()
for(i in 1:length(freq)){
  tables[[i]] <- data.frame(id = points$point_id[[i]],
                            unlist(freq[[i]]))
}

pixels <- rbindlist(tables) %>%
  mutate(area_m2 = Freq*30)

g <- f %>%
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

#save output table
write_csv(g, "./data/outputs/02_land_covers.csv")


#figure with buffers

plot(mapbiomas_paracambi)
plot(points, add = TRUE)

buf_pts <- buffer(points, 200)

plot(buf_pts, add = TRUE)

x <- crop(mapbiomas_paracambi, buf_pts)
x <- mask(mapbiomas_paracambi, buf_pts)

writeRaster(x, "./data/outputs/02_landcov_in_buffers.tif")

x <- mask(mapbiomas_paracambi, paracambi)

writeRaster(x, "./data/outputs/02_landcov_in_paracambi.tif")
