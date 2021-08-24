# Santos et al. (in preparation)
# Calculating correlations between sand flies and land use/cover
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(stringr)


# loading data ------------------------------------------------------------

spdata <- read_csv("./data/sandfly_data_paracambi_analysis.csv")

landcov <- read_csv("./outputs/02_land_covers_200m_buffer.csv")

study_spnames <- c("Nyssomyia intermedia*", "Nyssomyia whitmani*",
                   "Migonemyia migonei*", "Pintomyia fischeri*")


# calculate relative sand fly frequencies by point/year -------------------

# sampling effort
sampl <- spdata %>%
  filter(str_starts(point_year_id, "_") == FALSE) %>%
  group_by(point_year_id) %>%
  summarize(effort = n_distinct(month))

spp_freq <- spdata %>%
  filter(str_starts(point_year_id, "_") == FALSE,
         sp_code != "negative",
         sp_code != "L_sp") %>%
  group_by(point_year_id, sp_name) %>%
  summarize(n = sum(total)) %>%
  left_join(sampl, by = "point_year_id") %>%
  mutate(freq = round(n/effort, 2)) 


# join species and land cover data ----------------------------------------

spp_landcov <- spp_freq %>% 
  left_join(landcov, by = "point_year_id")


# calculate correlation ---------------------------------------------------

covers <- names(landcov)[-1]

colnamesr <- c(str_c(covers, "_r"))
colnamesp <- c(str_c(covers, "_p"))

correl_test <- matrix(nrow=length(study_spnames), ncol=length(covers)*2+1)
colnames(correl_test) <- c("species", colnamesr, colnamesp)
  
for(i in 1:length(study_spnames)){
  
  # filter data for this species
  data_cor_sp <- spp_landcov %>%
    filter(sp_name == study_spnames[[i]])
  correl_test[i,1] <- study_spnames[[i]]
  
  # calculate correlations with each land category
  for(j in 1:length(covers)){
    data_cor <- data_cor_sp %>%
      dplyr::select(sp_name, point_year_id, freq, covers[[j]]) # mudar aqui n ou freq para var dep da cor
    data_cor <- as.data.frame(data_cor)
    t <- cor.test(data_cor[,3], data_cor[,4], method = "pearson")
    correl_test[i,colnamesr[j]] <- round(t$estimate, digits = 3)
    correl_test[i,colnamesp[j]] <- round(t$p.value, digits = 3)
  }
  
}

correl_test <- as_tibble(correl_test)

# save final cor table
write_csv(correl_test, "./outputs/03_correl_test.csv")
