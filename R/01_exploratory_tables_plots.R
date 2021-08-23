# Santos et al. (in preparation)
# Exploring full data table and summarizing initial results
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)


# loading data ------------------------------------------------------------

spdata <- read_csv("./data/sandfly_data_paracambi_analysis.csv", 
                  col_types = cols(mes = col_date(format = "%d/%m/%Y"),
                                   ano = col_date(format = "%Y"),
                                   dout_mest = col_factor(levels = c("mestrado", "doutorado")))) 

spnames <- c("B. brumpti" = "Brumptomyia brumpti", 
                "B. nitzulescui" = "Brumptomyia nitzulescui", 
                "B. sp" = "Brumptomyia sp.", 
                "L. antunesi" = "Nyssomyia antunesi", 
                "L. cortelezzii" = "Evandromyia cortelezzii", 
                "L. fischeri" = "Pintomyia fischeri*", 
                "L. hirsuta" = "Psychodopygus hirsutus hirsutus ", 
                "L. intermedia" = "Nyssomyia intermedia*", 
                "L. migonei" = "Migonemyia migonei*", 
                "L. monticola" = "Pintomyia monticola", 
                "L. pelloni" = "Psathyromyia pelloni", 
                "L. pessoai" = "Pintomyia pessoai", 
                "L. schreiberi" = "Micropygomyia schreiberi", 
                "L. shannoni" = "Psathyromyia shannoni", 
                "L. sordelli" = "Sciopemyia sordellii", 
                "L. sp" = "Unidentified", 
                "L. whitmani" = "Nyssomyia whitmani*")


# species and counts by period --------------------------------------------

tab <- spdata %>%
  filter(especie != "negativa") %>%
  group_by(especie, dout_mest) %>%
  summarise("F" = sum(femeas),
            "M" = sum(machos),
            "Total" = sum(total)) %>%
  pivot_wider(names_from = dout_mest, 
              values_from = c("F", "M", "Total")) %>%
  replace(is.na(.), 0) %>%
  mutate(Total = Total_mestrado + Total_doutorado) %>%
  mutate(Species = recode(especie, !!!spnames)) %>%
  relocate(Species, .before = especie) %>%
  arrange(Species) %>%
  adorn_totals()
  
write_csv(tab, "./outputs/01_tab_periods.csv")


# comparing abundances by period ------------------------------------------

study_spnames <- c("L. intermedia", "L. migonei", "L. fischeri", "L. whitmani")

abund <- spdata %>%
  filter(especie %in% study_spnames) %>%
  mutate(Period = recode(dout_mest, "mestrado" = "1993-1994", "doutorado" = "2001-2003"),
         Species = recode(especie, !!!spnames))

tiff(filename = "./outputs/01_plot_period.tif",
     width = 1400, height = 1300, units = "px",
     res = 300,
     pointsize = 12,
     compression = "lzw")
ggplot(abund, aes(x = Species, y = total, fill = Period)) +
  geom_boxplot() +
  facet_wrap(~Species, scale="free") +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x=element_blank()) +
  xlab("") + ylab("Number of sand flies")
dev.off()

t_test <- matrix(nrow=5, ncol=6)
colnames(t_test) <- c("total_mest", "total_dout", "n_obs", "t", "df", "p-value")
rownames(t_test) <- c(study_spnames, "total")
for(i in 1:4){
  a <- filter(abund, especie == study_spnames[i], dout_mest == "mestrado")
  d <- filter(abund, especie == study_spnames[i], dout_mest == "doutorado")
  t <- t.test(a$total, d$total)
  t_test[i,1] <- sum(a$total)
  t_test[i,2] <- sum(d$total)
  t_test[i,3] <- nrow(a) + nrow(d)
  t_test[i,4] <- round(t$statistic, digits = 3)
  t_test[i,5] <- round(t$parameter, digits = 3)
  t_test[i,6] <- round(t$p.value, digits = 3)
}
c <- filter(spdata, dout_mest == "mestrado")
e <- filter(spdata, dout_mest == "doutorado")
t <- t.test(c$total, e$total)
t_test[5,1] <- sum(c$total)
t_test[5,2] <- sum(e$total)
t_test[5,3] <- nrow(c) + nrow(e)
t_test[5,4] <- round(t$statistic, digits = 3)
t_test[5,5] <- round(t$parameter, digits = 3)
t_test[5,6] <- round(t$p.value, digits = 3)

t_test <- data.frame(especie = row.names(t_test), t_test)
write_csv(t_test, "./outputs/01_tab_t-test.csv")


# counts by sex -----------------------------------------------------------

tab1 <- spdata %>%
  filter(dout_mest == "doutorado",
         especie != "negativa") %>%
  group_by(especie) %>%
  summarize(femeas = sum(femeas),
            machos = sum(machos),
            total = sum(total)) %>%
  arrange(desc(total))

write_csv(tab1, "./outputs/01_tab_sex.csv")



# monthly abundance -------------------------------------------------------

spp_total_mes <- abund %>%
  filter(dout_mest == "doutorado") %>%
  dplyr::select(especie, mes, total) %>%
  filter(especie != "negativa") %>%
  drop_na() %>%
  group_by(especie, mes) %>%
  summarise(total = sum(total))

# all species in a single plot

levels(spp_total_mes$especie) <- study_spnames

tiff(filename = "./outputs/01_plot_seasonality.tif",
     width = 1400, height = 1800, units = "px",
     res = 300,
     pointsize = 12,
     compression = "lzw")
ggplot(spp_total_mes, aes(x = mes, y = total)) +
  geom_line() + xlab("") + ylab("Number of sand flies") + 
  geom_vline(xintercept = 16283, linetype="dashed", color = "gray") +
  facet_wrap(~ especie, nrow = 4) +
  theme(strip.text = element_text(face = "italic"))
dev.off()

# few consistent observations to use this seasonality plot


# counts by trap type -----------------------------------------------------

tab2 <- spdata %>%
  filter(#dout_mest == "doutorado",
         especie != "negativa") %>%
  group_by(especie, tipo) %>%
  summarise(total = sum(total)) %>%
  drop_na() %>%
  pivot_wider(names_from = tipo, 
              values_from = total) %>%
  replace_na(list(luminosa = 0, manual = 0)) %>%
  mutate(total = luminosa + manual) %>%
  arrange(desc(total)) %>%
  adorn_totals()

write_csv(tab2, "./outputs/01_tab_trap_type.csv")


plot2data <- spdata %>%
  filter(especie != "negativa",
         especie != "L. sp",
         especie != "B. sp",
         tipo != "NA") %>%
  group_by(especie, tipo) %>%
  summarise(total = sum(total)) %>%
  mutate(percent = total/sum(total),
         Species = recode(especie, !!!spnames),
         Method = recode(tipo, "luminosa" = "Light trap", "manual" = "Manual capture")) %>%
  drop_na()


tiff(filename = "./outputs/01_plot_trap_type.tif",
     width = 2000, height = 1400, units = "px",
     res = 300,
     pointsize = 10,
     compression = "lzw")
ggplot(plot2data, aes(x=reorder(Species, desc(Species)), y=percent, fill=Method)) +
  geom_bar(stat="identity", width = .7, lwd=0.1) +
  theme_minimal()+
  geom_text(aes(label=paste0(sprintf("%.0f", percent*100),"%")),
            position=position_stack(vjust=0.5), colour="white", lwd = 2.5) +
  coord_flip() +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(y="", x="")
dev.off()


# counts by urban/rural zone ----------------------------------------------

tab3 <- spdata %>%
  filter(#dout_mest == "doutorado",
         especie != "negativa") %>%
  group_by(especie, area) %>%
  summarise(total = sum(total)) %>%
  drop_na() %>%
  pivot_wider(names_from = area, 
              values_from = total) %>%
  replace_na(list(semiurbano = 0, rural = 0)) %>%
  mutate(total = semiurbano + rural) %>%
  arrange(desc(total))

write_csv(tab3, "./outputs/01_tab_zones.csv")


plot3data <- spdata %>%
  filter(especie != "negativa",
         especie != "L. sp",
         especie != "B. sp",
         area != "NA",
         especie != "L. antunesi") %>%
  group_by(especie, area) %>%
  summarise(total = sum(total)) %>%
  mutate(percent = total/sum(total),
         Species = recode(especie, !!!spnames),
         Zone = recode(area, "rural" = "Rural", "semiurbano" = "Semi-urban")) %>%
  drop_na()

tiff(filename = "./outputs/01_plot_zones.tif",
     width = 2000, height = 1400, units = "px",
     res = 300,
     pointsize = 10,
     compression = "lzw")
ggplot(plot3data, aes(x=reorder(Species, desc(Species)), y=percent, fill=Zone)) +
  geom_bar(stat="identity", width = .7, lwd=0.1) +
  theme_minimal()+
  geom_text(aes(label=paste0(sprintf("%.0f", percent*100),"%")),
            position=position_stack(vjust=0.5), colour="white", lwd = 2.5) +
  coord_flip() +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(y="", x="")
dev.off()


# counts by distance to forest border -------------------------------------

tab4 <- spdata %>%
  filter(dout_mest == "doutorado",
         especie != "negativa") %>%
  group_by(especie, distancia_casa_mata) %>%
  summarise(total = sum(total)) %>%
  pivot_wider(names_from = distancia_casa_mata, 
              values_from = total) %>%
  replace_na(list("0-200" = 0,
                  "200-500" = 0,
                  ">1000" = 0,
                  "500-1000" = 0))

write_csv(tab4, "./outputs/01_tab_distance_forest.csv")



# counts by ecotope -------------------------------------------------------

ecotopes <- c("intra_isca_humana" = "intra", 
             "intra_teto" = "intra", 
             "parede_interna" = "intra", 
             "bananal_isca_humana" = "isca_humana", 
             "bananal_tronco" = "bananal", 
             "floresta_tronco" = "arvore", 
             "floresta_isca_humana" = "isca_humana", 
             "peri_isca_humana" = "isca_humana", 
             "tronco" = "arvore", 
             "anexo" = "parede_externa")

tab5 <- spdata %>%
  filter(#dout_mest == "doutorado",
         especie != "negativa",
         especie != "L. sp",
         local != "A/GL",
         local != "BN/CR",
         local != "NA") %>%
  mutate(Ecotope = recode(local, !!!ecotopes),
         Species = recode(especie, !!!spnames)) %>%
  filter(Ecotope != "isca_humana") %>%
  group_by(Species, Ecotope) %>%
  summarise(total = sum(total)) %>%
  pivot_wider(names_from = Ecotope, 
              values_from = total) %>%
  replace(., is.na(.), 0) %>%
  relocate("intra", .after = "Species") %>%
  relocate("parede_externa", .after = "intra") %>%
  #select(-"isca_humana") %>%
  arrange("Species")
    
tab5_percent <- tab5 %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 0) %>%
  adorn_ns(position = "front")

write_csv(tab5, "./outputs/01_tab_ecotopes.csv")
write_csv(tab5_percent, "./outputs/01_tab_ecotopes_percent.csv")
