# Santos et al. (in preparation)
# Exploring sand fly data and summarizing results
# Bruno M. Carvalho
# brunomc.eco@gmail.com

library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(ggplot2)


# loading data ------------------------------------------------------------

spdata <- read_csv("./data/sandfly_data_paracambi_analysis.csv", 
                  col_types = cols(month = col_date(format = "%d/%m/%Y"),
                                   year = col_date(format = "%Y"),
                                   period = col_factor(levels = c("1990s", "2000s"))))

study_spcodes <- c("intermedia", "migonei", "fischeri", "whitmani")
study_spnames <- c("Nyssomyia intermedia*", "Migonemyia migonei*", 
                   "Pintomyia fischeri*", "Nyssomyia whitmani*")


# species and counts by period --------------------------------------------

tab <- spdata %>%
  filter(sp_code != "negative") %>%
  group_by(sp_name, period) %>%
  summarise("F" = sum(females),
            "M" = sum(males),
            "Total" = sum(total)) %>%
  pivot_wider(names_from = period, 
              values_from = c("F", "M", "Total")) %>%
  replace(is.na(.), 0) %>%
  mutate(Total = Total_1990s + Total_2000s) %>%
  arrange(sp_name) %>%
  adorn_totals()
  
write_csv(tab, "./outputs/01_tab_periods.csv")


# comparing abundances by period ------------------------------------------

abund <- spdata %>%
  filter(sp_code %in% study_spcodes) %>%
  rename(Species = sp_name) %>%
  mutate(Period = recode(period, "1990s" = "1992-1994", "2000s" = "2001-2003")) 

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


# testing for difference in abundande/periods -----------------------------

t_test <- matrix(nrow = length(study_spcodes), ncol=7)
colnames(t_test) <- c("species", "total_1990s", "total_2000s", 
                      "n_obs", "t", "df", "p-value")
for(i in 1:length(study_spcodes)){
  a <- filter(abund, sp_code == study_spcodes[i], period == "1990s")
  d <- filter(abund, sp_code == study_spcodes[i], period == "2000s")
  t <- t.test(a$total, d$total)
  t_test[i,1] <- study_spnames[i]
  t_test[i,2] <- sum(a$total)
  t_test[i,3] <- sum(d$total)
  t_test[i,4] <- nrow(a) + nrow(d)
  t_test[i,5] <- round(t$statistic, digits = 3)
  t_test[i,6] <- round(t$parameter, digits = 3)
  t_test[i,7] <- round(t$p.value, digits = 3)
}

write_csv(as_tibble(t_test), "./outputs/01_tab_t-test.csv")


# monthly abundance -------------------------------------------------------
month_levels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

monthly_effort <- tibble(month_number = seq(1:12),
                         month_name = factor(c("Jan", "Feb", "Mar", "Apr",
                                        "May", "Jun", "Jul", "Aug",
                                        "Sep", "Oct", "Nov", "Dec"),
                                        levels = month_levels),
                         trapnights = c(6,3,5,6,4,5,6,6,6,6,6,6))

monthly_abund_1990s <- abund %>%
  mutate(month_number = as.numeric(format(month, "%m"))) %>%
  group_by(Species, month_number) %>%
  summarise(total = sum(total)) %>%
  left_join(monthly_effort, by = "month_number") %>%
  mutate(rel_abund = total/trapnights)

write_csv(monthly_abund_1990s, "./outputs/01_monthly_abund_1990s.csv")

# all species in a single plot

levels(monthly_abund_1990s$Species) <- study_spnames

tiff(filename = "./outputs/01_plot_seasonality.tif",
     width = 1400, height = 1800, units = "px",
     res = 300,
     pointsize = 12,
     compression = "lzw")
ggplot(monthly_abund_1990s, aes(x = month_name, y = rel_abund)) +
  geom_col() + xlab("") + ylab("Relative abundance") + 
  facet_wrap(~ Species, nrow = 4, scales = 'free') +
  theme(strip.text = element_text(face = "italic"))
dev.off()

# few consistent observations to use the seasonality plot for 2000s


# counts by trap type -----------------------------------------------------

tab2 <- spdata %>%
  filter(sp_code != "negative") %>%
  group_by(sp_name, trap_type) %>%
  summarise(total = sum(total)) %>%
  drop_na() %>%
  pivot_wider(names_from = trap_type, 
              values_from = total) %>%
  replace_na(list(light = 0, manual = 0)) %>%
  mutate(total = light + manual) %>%
  arrange(desc(total)) %>%
  adorn_totals()

write_csv(tab2, "./outputs/01_tab_trap_type.csv")

plot2data <- spdata %>%
  filter(sp_code != "negative",
         sp_code != "L_sp",
         sp_code != "B_sp",
         trap_type != "NA") %>%
  group_by(sp_name, trap_type) %>%
  summarise(total = sum(total)) %>%
  rename(Species = sp_name) %>%
  mutate(percent = total/sum(total),
         Method = recode(trap_type, "light" = "Light trap", "manual" = "Manual capture")) %>%
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
  filter(sp_code != "negative") %>%
  group_by(sp_name, zone) %>%
  summarise(total = sum(total)) %>%
  drop_na() %>%
  pivot_wider(names_from = zone, 
              values_from = total) %>%
  replace_na(list(semiurban = 0, rural = 0)) %>%
  mutate(total = semiurban + rural) %>%
  arrange(desc(total))

write_csv(tab3, "./outputs/01_tab_zones.csv")


plot3data <- spdata %>%
  filter(sp_code != "negative",
         sp_code != "L_sp",
         sp_code != "B_sp",
         zone != "NA") %>%
  group_by(sp_name, zone) %>%
  summarise(total = sum(total)) %>%
  rename(Species = sp_name) %>%
  mutate(percent = total/sum(total),
         Zone = recode(zone, "rural" = "Rural", "semiurban" = "Semiurban")) %>%
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
  filter(sp_code != "negative") %>%
  group_by(sp_name, distance_to_forest) %>%
  summarise(total = sum(total)) %>%
  pivot_wider(names_from = distance_to_forest, 
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
  filter(sp_code != "negative",
         sp_code != "L_sp",
         ecotope != "mixed_peridomestic",
         ecotope != "NA") %>%
  group_by(sp_name, ecotope) %>%
  summarise(total = sum(total)) %>%
  pivot_wider(names_from = ecotope, 
              values_from = total) %>%
  replace(., is.na(.), 0) %>%
  relocate(internal_wall, .after = sp_name) %>%
  relocate(external_wall, .after = internal_wall)

tab5_percent <- tab5 %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting(digits = 0) %>%
  adorn_ns(position = "front")

write_csv(tab5, "./outputs/01_tab_ecotopes.csv")
write_csv(tab5_percent, "./outputs/01_tab_ecotopes_percent.csv")
