#kernel?

library(spatialEco)
library(readr)
library(dplyr)
library(tidyr)
library(raster)
library(rgdal)

#spdata
dados <- read_csv("./data/flebs_paracambi_ginelza_completo_check.csv")
spp_abund <- dados %>%
  filter(dout_mest == "doutorado",
         especie != "negativa") %>%
  group_by(point_id, especie) %>%
  summarize(n = sum(total)) %>%
  pivot_wider(names_from = especie, 
              values_from = n) %>%
  replace(., is.na(.), 0) %>%
  rename("L.int" = "L. intermedia",
         "L.fis" = "L. fischeri",
         "L.mig" = "L. migonei",
         "L.pel" = "L. pelloni",
         "B.bru" = "B. brumpti",
         "L.whi" = "L. whitmani",
         "L.ant" = "L. antunesi",
         "L.cor" = "L. cortelezzii",
         "L.sor" = "L. sordelli",
         "L.sch" = "L. schreiberi",
         "B.nit" = "B. nitzulescui") %>%
  mutate(total = sum(c_across(L.int:L.cor)),
         log_total = log(total))

#add coords

points <- dados %>%
  filter(ano == 2002) %>%
  dplyr::select(point_id, ano, lon_campo, lat_campo) %>%
  distinct() %>%
  left_join(spp_abund, by = "point_id") %>%
  drop_na()

CRS = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

spatial <- SpatialPointsDataFrame(coords = data.frame(points$lon_campo, points$lat_campo),
                       data = spp_abund,
                       proj4string = CRS)

writeOGR(spatial, "./data/outputs/spp_abund_doutorado.shp", driver = "ESRI Shapefile", layer = "spp_abund")

#load rasters
mapbiomas <- raster("./data/raster/mapbiomas-brazil-collection-50-riodejaneiro-2002.tif")
paracambi <- readOGR("./data/shp/paracambi.shp")
mapbiomas_paracambi <- crop(mapbiomas, extent(paracambi)+0.05)


#check
plot(mapbiomas_paracambi)
plot(paracambi, add = TRUE)
plot(spatial, add = TRUE)

#kernel

weighted <- sp.kde(spatial, y = spatial$log_total, bw = 0.02, newdata = mapbiomas_paracambi)

#check
plot(weighted)
plot(paracambi, add = TRUE)
plot(spatial, add = TRUE)
writeRaster(weighted, "./results/kernel_0-02_log.tif")

#variogram
library(gstat)

lzn.vgm <- variogram(log(total)~1, data=a) # calculates sample variogram values 
lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph", 900, 1)) # fit model

plot(lzn.vgm)
