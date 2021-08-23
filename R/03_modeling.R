# checking effects of land cover in sand fly species via cca

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

dados <- read_csv("./data/flebs_paracambi_ginelza_completo_check.csv")
landcov <- read_csv("./data/outputs/02_land_covers_mod.csv")

# tudo daqui pra frente é somente com dados do doutorado, porque as capturas do mestrado nao tem coordenadas
# calcular esforço amostral por ponto (numero de dias/noites amostragem)

sampl <- dados %>%
  filter(dout_mest == "doutorado") %>%
  group_by(point_id) %>%
  summarize(dias = n_distinct(mes))

# calcular abundancia relativa

spp_abund <- dados %>%
  filter(dout_mest == "doutorado") %>%
  group_by(point_id, especie) %>%
  summarize(n = sum(total)) %>%
  left_join(sampl, by = "point_id") %>%
  mutate(abundr = round(n/dias, 2))


# gerar matriz de especies e abundancia

spp <- spp_abund %>%
  select(-abundr, -dias) %>%
  pivot_wider(names_from = especie, 
              values_from = n) %>%
  replace(., is.na(.), 0) %>%
  select(-negativa)

spp_rel <- spp_abund %>%
  select(-n, -dias) %>%
  pivot_wider(names_from = especie, 
              values_from = abundr) %>%
  replace(., is.na(.), 0) %>%
  select(-negativa)


# combinar com land covers e coords

landcov <- landcov %>%
  filter(year == 2002)

spp_landcov <- spp %>% #mudar aqui spp ou spp_rel para abundancias por esforço de coleta
  left_join(landcov, by = "point_id")

coords <- dados %>%
  select(point_id, lon_campo, lat_campo) %>%
  distinct()

data_cca <- spp_landcov %>%
  left_join(coords, by = "point_id")
  #rename("L.int" = "L. intermedia",
  #       "L.fis" = "L. fischeri",
  #       "L.mig" = "L. migonei",
  #       "L.pel" = "L. pelloni",
  #       "B.bru" = "B. brumpti",
  #       "L.whi" = "L. whitmani",
  #       "L.ant" = "L. antunesi",
  #       "L.cor" = "L. cortelezzii",
  #       "L.sor" = "L. sordelli",
  #       "L.sch" = "L. schreiberi",
  #       "B.nit" = "B. nitzulescui")


# inspect for visual correlations:

ggplot(data = data_cca) +
  geom_point(mapping = aes(x = forest, y = L.int))

ggplot(data = data_cca) +
  geom_point(mapping = aes(x = urban, y = L.int))

ggplot(data = data_cca) +
  geom_point(mapping = aes(x = nonveg, y = L.int))

ggplot(data = data_cca) +
  geom_point(mapping = aes(x = pasture, y = L.fis))

ggplot(data = data_cca) +
  geom_point(mapping = aes(x = agripasture, y = L.int))

# test pearson correlations:

data_cor <- data.frame(data_cca)
data_cor <- cbind(data_cor[,2:6], data_cor[,15:19])

study_spnames <- c("L.int", "L.mig", "L.fis", "L.whi")
cob <- c("forest", "pasture", "agripasture", "nonveg", "urban")

teste <- matrix(nrow=10, ncol=4)
rownames(teste) <- c("floresta_r", "floresta_p",
                     "pasto_r","pasto_p",
                     "agropecuaria_r","agropecuaria_p",
                     "sem_vegetacao_r","sem_vegetacao_p",
                     "urbano_r", "urbano_p")
colnames(teste) <- study_spnames
for(i in 1:4){
  t <- cor.test(data_cor[,i], data_cor[,6], method = "pearson")
  teste[1,i] <- round(t$estimate, digits = 3)
  teste[2,i] <- round(t$p.value, digits = 3)
  t <- cor.test(data_cor[,i], data_cor[,7], method = "pearson")
  teste[3,i] <- round(t$estimate, digits = 3)
  teste[4,i] <- round(t$p.value, digits = 3)
  t <- cor.test(data_cor[,i], data_cor[,8], method = "pearson")
  teste[5,i] <- round(t$estimate, digits = 3)
  teste[6,i] <- round(t$p.value, digits = 3)
  t <- cor.test(data_cor[,i], data_cor[,9], method = "pearson")
  teste[7,i] <- round(t$estimate, digits = 3)
  teste[8,i] <- round(t$p.value, digits = 3)
  t <- cor.test(data_cor[,i], data_cor[,10], method = "pearson")
  teste[9,i] <- round(t$estimate, digits = 3)
  teste[10,i] <- round(t$p.value, digits = 3)
}
teste <- data.frame(especie = row.names(teste), teste)

write_csv(teste, "./results/03_pearson_cor.csv")

# test some regressions

m <- lm(L.int ~ forest + pasture + agripasture + nonveg + urban, data = data_cca)
ms <- summary(m)
m2 <- step(m)
m2s <- summary(m2)


m1 <- glm(L.int ~ forest + pasture + agripasture + nonveg + urban, data = data_cca, family = "poisson")
summary(m1)

# com poisson: todos coeficientes iguais = talvez porque coberturas somam 1?


#### CCA (finalmente) ####

#dividir as tres matrizes:

data_cca2 <- data_cca %>%
  mutate(total = sum(c_across(L.int:L.cor))) %>%
  filter(total != 0)

spe <- data_cca2[,2:12]
env <- data_cca2[,15:19]
spatial <- data_cca2[,20:21]

# Apply log+1 transformation to your species occurrences data (spe matrix) in order to correct for possible statistical errors associated to rare or very common species:

spelog <- decostand(spe, "log")

## Perform CCA. We need to specify that spe (species distribution matrix) is explained by env (environmental matrix).

ccamodel <- cca(spe~., env)

## To perform a Partial CCA, we need to use a third matrix (conditioning matrix). To do that, we must first combine variables from the environmental ("env") and conditional ("spatial") matrices. We have to do that in order to later apply another function ("ordistep" function). After we have combined all variables, we will apply the Partial CCA using a formula where we specify the response ("spe" matrix), the constraint variables (each of the variables from "env" matrix), and the conditioning variables (variables from the conditioning matrix; in our case "spatial" matrix)

envspatial<-cbind(env,spatial)

nams <- names(envspatial)

partialccamodel <- formula(paste("spe ~", paste(nams[1: (length(envspatial)-(length(spatial)) )], collapse = " + "),"+ Condition(", paste(nams[(length(envspatial)-(length(spatial)-1) ):length(envspatial)], collapse ="+"),")"))

partialccamodel<-cca(partialccamodel, envspatial)

# However, we also have to automatically select variables of "env" matrix that best explain "spe" matrix. We can do that by using a stepwise model from "ordistep" function. Let us do that with our "ccamodel" (not the partial cca).

finalmodel<- ordistep(ccamodel, scope=formula(ccamodel))

# Then, we can calculate Variance Inflation Factors (VIF) for each of the constraints (variables) from the "env" matrix (environmental matrix). If we find an environmental variable with VIF>10, we'll know that this variable presents colinearity with another or other variables. In that case, we would have to delete the variable from our initial dataset and redo all the analysis. In our example, no variable is redundant with each other (all of them have VIF<10).

vif.cca(finalmodel)

# Let us now fit the stepwise model (ordistep function) using our partial cca model ("partialccamodel"). Note that we will use X (longitude) and Y (latitude) variables from the previously created "envspatial" object as our conditioning variables. After "vif.cca", note that "X" (spatial) variable has a value > 10, so one should consider to delete the variable before the analysis.

partialccamodel <- formula(paste("spe ~", paste(nams[1: (length(envspatial)-(length(spatial)) )], collapse = " + "),"+ Condition(", paste(nams[(length(envspatial)-(length(spatial)-1) ):length(envspatial)], collapse ="+"),")"))

simplemodel<-cca(partialccamodel, envspatial)

finalmodelpartial<- ordistep(simplemodel, scope=formula(partialccamodel))

vif.cca(finalmodelpartial)

### Ok, let us now call "ccamodel" object (not the "partialccamodel") to see how to interpret results.

finalmodelpartial


## Note that "Total Inertia" is the total variance in species (observations matrix) distributions. "Constrained Inertia" is the variance explained by the environmental variables (gradients matrix). The "Proportion" values represent the percentages of variance of species distributions explained by Constrained (environmental) and Unconstrained variables. Eigenvalues of constrained and unconstrained axes represent the amount of variance explained by each CCA axis (graphs usually present the first two constrained axes, so take a look at their values).


# This is a critical step when doing the CCA. When we have our final model, we must use permutation tests to observe if our whole CCA model, the CCA terms (environmental varibles), and CCA axes explain more variance of "spe" (observations) matrix than expected by chance (tests should be significant; p values lower or equal 0.05). If the tests are not significant, there is no point in using the CCA. In our example, we'll see that we can continue using the CCA results.

# Testing the significance of the CCA model:
anova.cca(finalmodel)

# Testing the significance of terms (environmental variables):
anova.cca(finalmodel, by="terms")

# Testing the significance of CCA axes (at least the first two or three should present a significant p value):
anova.cca(finalmodel, by="axis")


### Finally, we may want to generate graphs in order to better understand results. To do that, we have to take a look at "xlim" (limits of x axis), "ylim" (limits of y axis), and "display" (if we want to observe species, environmental gradients, and/or sites in the graph) arguments. For example, to show only species scores along CCA axes, use only "sp" inside display" argument.

plot(finalmodel, xlim=c(-1.5,2), ylim=c(-1,1.5), display=c("sp"))

# If you want to show species and environmental gradients in our plot, use both "sp" and "cn" inside "display" argument,

plot(finalmodel, xlim=c(-3,3), ylim=c(-3,3), display=c("sp","cn"))

# To show species, environmental gradients, and sites, use "sp", "cn", and "wa" inside "display" argument,

plot(finalmodel, xlim=c(-3,3), ylim=c(-3,3), display=c("sp","cn","wa"))



#### Por enquanto tá tudo uma merda. Essas porcentagens de coberturas são obviamente todas correlacionadas, então tem que repensar o uso delas como variáveis explicativas