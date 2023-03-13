### Titre -------------------------------------
# Nom : Calcul indices overlap
# Auteure : Perle Charlot
# Date de création : 10-03-2023
# Dates de modification : 13-03-2023

### Librairies -------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
### Fonctions -------------------------------------

### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")
input_path <- paste0(wd,"/input/")

#### Données spatiales ####

#### Autre ####
liste.mois = c("mai","juin","juillet","aout","septembre")
df.mois = data.frame(nom_mois = liste.mois, numero_mois = c("05","06","07","08","09"))
# # Liste dimensions
# liste.dim =  c("CA","B","PV","CS","D","I")
# Liste usages
liste.usages = c("Ni","Lk","Co","Pa","Rp","Vt")

### Programme -------------------------------------

#### E-space ####

##### Schoener D abond obs #####
# corrigées par disponibilité (Broennimann 2012), env grid

#TEST
usage = "Ni"
mois = "juin"
type_donnees = "ACP_avec_ponderation" # "ACP_ACP" ou "ACP_avec_ponderation" ou
# "ACP_sans_ponderation" ou "brute"
fit = "2_axes" # "2_axes" ou all_simple"
algorithme = "glm"
chemin_esp_eco = paste0(output_path,"/niches/",type_donnees,"/",usage,"/",
                        fit,"/predictions_",algorithme,"/espace_eco/")

# import data pour dessiner contour niche
load(paste0(chemin_esp_eco,"/niche_potentielle.rdata")) # load grid_usage_rdata limits

chemin_pred = paste0(output_path,"/niches/",type_donnees,"/",usage,"/",
                     fit,"/predictions_",algorithme,"/")

dt_uses_env <- fread(paste0(chemin_pred,"/dt_probUs_condiEnv_",mois,".csv"), dec=",")
dt_test = dt_uses_env[,1:5]
names(dt_test) = c("obs_usage","proba_presence","pred_presence","axe1","axe2")
# Grid 100 * 100
dt_test2 = dt_test %>% mutate(
  cut_x = cut(axe1, breaks = seq(from = limits[1], to = limits[2], length.out = 100), include.lowest = T),
  cut_y = cut(axe2, breaks = seq(from = limits[3], to = limits[4], length.out = 100), include.lowest = T)
) %>%
  group_by(cut_x, cut_y) %>% 
  mutate(n_bin = n(), 
         mean_proba = mean(proba_presence),
         sd_proba = sd(proba_presence) )


dt_grid_test = dt_test %>% mutate(
  cut_x = cut(axe1, breaks = seq(from = limits[1], to = limits[2], length.out = 100), include.lowest = T),
  cut_y = cut(axe2, breaks = seq(from = limits[3], to = limits[4], length.out = 100), include.lowest = T)
) %>%
  group_by(cut_x, cut_y) %>% 
  summarise(n_bin = n(), 
         mean_proba = mean(proba_presence),
         sd_proba = sd(proba_presence),
         sum_obs = sum(obs_usage))

# densité condition env pour les présences observées + contour niche 
grid_usage2 = grid_usage_rdata %>% filter(proba_presence >= seuil)
grid_usage2 = grid_usage_rdata

dt_test2$obs_usage[dt_test2$obs_usage == 0.5] <- 1

P = dt_test2 %>% 
  # filter(obs_usage == "1") %>%
  ggplot(aes(axe1, axe2, color=n_bin, alpha=n_bin)) +
  geom_point()+
  scale_color_distiller(palette ="Spectral")+
  xlim(limits[1], limits[2])+
  ylim(limits[3], limits[4])+
  stat_contour_filled(data=grid_usage2,
                      aes(x=axe1, y=axe2, z=proba_presence),
                      color="black",
                      size=0.55, bins=2,
                      show.legend =T,
                      alpha=0.1)+
  #scale_fill_manual(values=c("transparent")) +
  theme_minimal() +
  guides(alpha = FALSE)+
  theme(text = element_text(size=15)) +
  facet_grid(~ as.factor(obs_usage))

  #  +
  # labs(y="Environmental axe 2",x="Environmental axe 1",
  #      title =paste0(usage, " - ",mois),
  #      subtitle = paste0("density of realized environmental space, for observed presences, with niche contour"),
  #      color="Density")+
  # theme(
  #   legend.position = c(.95, .95),
  #   legend.justification = c("right", "top"),
  #   legend.box.just = "right",
  #   legend.margin = margin(6, 6, 6, 6)
  # )

P2 = dt_grid_test %>% 
  # filter(obs_usage == "1") %>%
  ggplot(aes(cut_x, cut_y, fill=n_bin)) +
  geom_raster() +
  scale_fill_distiller(palette ="Spectral")


# scale des densités
MAX_N = max(dt_grid_test$n_bin)
MAX_nobs = max(dt_grid_test$sum_obs)

dt_grid_test$e_ij <- dt_grid_test$n_bin / MAX_N
dt_grid_test$o_ij <- dt_grid_test$sum_obs / MAX_nobs
max_oe = max(dt_grid_test$o_ij/dt_grid_test$e_ij)
dt_grid_test$z_ij <- (dt_grid_test$o_ij / dt_grid_test$e_ij) / max_oe

# z = 0 lorsque obs < 3
stock = dt_grid_test

dt_grid_test$z_ij[dt_grid_test$sum_obs < 3 & dt_grid_test$z_ij == 1] <- 0

# available environment
Pe = dt_grid_test %>% 
  ggplot(aes(cut_x, cut_y, fill=e_ij, alpha = e_ij)) +
  geom_raster() +
  scale_fill_distiller(palette ="Spectral")
# density observation
Po = dt_grid_test %>% 
  ggplot(aes(cut_x, cut_y, fill=o_ij, alpha=o_ij)) +
  geom_raster() +
  scale_fill_distiller(palette ="Spectral")
# densité obs en fonction disponibilité env
Poe = dt_grid_test %>% 
  ggplot(aes(cut_x, cut_y, fill=z_ij, alpha=z_ij)) +
  geom_raster() +
  scale_fill_distiller(palette ="Spectral")

# bugs : z = 1 ? pour des pixels où peu d'obs + env rare
# quand pixel env le plus rare avec 1 seule observation de présence -> z = 1

library(patchwork)
Pe / Po  /Poe


#####  Schoener D sur probabilités issues SDM, projetées dans esp env ##### 

#####  % surface niche binarisé (en fonction cutoff) ##### 

#### G-space ####

##### % surface obs #####
# ("témoin", ce qui est basiquement fait)

##### Schoener D proba occ ##### 
# projetées dans esp géographique

#####  % surface pred occ ##### 
# (en fonction cutoff)
