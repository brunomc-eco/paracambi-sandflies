# Santos et al. (in preparation)
# Calculating correlations between sand flies and land use/cover
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


# loading data ------------------------------------------------------------

spdata <- read_csv("./data/sandfly_data_paracambi_analysis.csv")

#landcov <- read_csv("./data/outputs/02_land_covers_mod.csv")
landcov <- read_csv("./outputs/02_land_covers.csv")


# calculate relative sand fly frequencies by point/year -------------------

# sampling effort in days
sampl <- spdata %>%
  group_by(point_year_id) %>%
  summarize(dias = n_distinct(mes))

spp_freq <- spdata  %>%
  group_by(point_year_id, especie) %>%
  summarize(n = sum(total)) %>%
  left_join(sampl, by = "point_year_id") %>%
  mutate(freq = round(n/dias, 2)) %>%
  filter(str_starts(point_year_id, "_") == FALSE,
         especie != "negativa",
         especie != "L. sp")


# join species and land cover data ----------------------------------------

spp_landcov <- spp_freq %>% 
  left_join(landcov, by = "point_year_id")


# calculate correlation ---------------------------------------------------

#spnames <- unique(spp_landcov$especie)
spnames <- c("L. intermedia", "L. migonei", "L. fischeri", "L. whitmani")
covs <- c("forest", "pasture", "agripasture", "nonveg", "urban")
covs <- c("forest_f", "pasture_f", "agripasture_f", "nonveg_f", "urban_f")

colnamesr <- c(str_c(covs, "_r"))
colnamesp <- c(str_c(covs, "_p"))

correl_test <- matrix(nrow=length(spnames), ncol=(length(covs)*2)+1)
colnames(correl_test) <- c("species", colnamesr, colnamesp)
  
for(i in 1:length(spnames)){
  
  # filter data for this species
  data_cor_sp <- spp_landcov %>%
    filter(especie == spnames[[i]])
  correl_test[i,1] <- spnames[[i]]
  
  # calculate correlations with each land category
  for(j in 1:5){
    data_cor <- data_cor_sp %>%
      dplyr::select(especie, point_year_id, freq, covs[[j]]) # mudar aqui n ou freq para var dep da cor
    data_cor <- as.data.frame(data_cor)
    t <- cor.test(data_cor[,3], data_cor[,4], method = "pearson")
    correl_test[i,colnamesr[j]] <- round(t$estimate, digits = 3)
    correl_test[i,colnamesp[j]] <- round(t$p.value, digits = 3)
  }
  
}

correl_test <- as_tibble(correl_test)

# save final cor table
write_csv(correl_test, "./outputs/03_correl_test.csv")
