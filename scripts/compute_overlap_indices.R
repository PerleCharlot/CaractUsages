### Titre -------------------------------------
# Nom : Calcul indices overlap
# Auteure : Perle Charlot
# Date de création : 10-03-2023
# Dates de modification : 15-03-2023

### Librairies -------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
### Fonctions -------------------------------------

# Compute overlap D schoener between 2 uses
DSchoener <- function(i, j){1-(0.5*(sum(abs(i - j), na.rm=T)))}

# Fonction qui retourne un dataframe avec o_ij e_ij et z_ij pour un usage, pour un mois
GridObs <- function(usage,  mois,
                    type_donnees = "ACP_avec_ponderation", 
                    algorithme = "glm", fit = "2_axes"){
  
  # #TEST
  # usage = liste.usages[1]
  # mois = liste.mois[2]
  # type_donnees = "ACP_avec_ponderation" # "ACP_ACP" ou "ACP_avec_ponderation" ou
  # # "ACP_sans_ponderation" ou "brute"
  # fit = "2_axes" # "2_axes" ou all_simple"
  # algorithme = "glm"
  
  # Conserver les chemins
  chemin_esp_eco = paste0(output_path,"/niches/",type_donnees,"/",usage,"/",
                          fit,"/predictions_",algorithme,"/espace_eco/")
  chemin_pred = paste0(output_path,"/niches/",type_donnees,"/",usage,"/",
                       fit,"/predictions_",algorithme,"/")
  
  # import data pour dessiner contour niche
  load(paste0(chemin_esp_eco,"/niche_potentielle.rdata")) # load grid_usage_rdata limits seuil
  # import data du mois en cours
  dt_uses_env <- fread(paste0(chemin_pred,"/dt_probUs_condiEnv_",mois,".csv"), dec=",")
  names(dt_uses_env)[1:3] = c("obs_usage","proba_presence","pred_presence")
  # correction noms colonnes
  if(fit == "2_axes"){
    ncols = 2
    ind_endaxe = grep("axe2", names(dt_uses_env))
    ind_x = grep("x$", names(dt_uses_env))
    ind_y = grep("y$", names(dt_uses_env))
    dt_uses_env = dt_uses_env[,c(1:ind_endaxe,ind_x, ind_y), with=FALSE]
    names(dt_uses_env)[4:grep("axe2", names(dt_uses_env))] = paste0("axe",seq(1,ncols))
  }else{
    ncols = length(grep("axe", names(dt_uses_env)))
    names(dt_uses_env)[4:grep(paste0("axe",ncols), names(dt_uses_env))] = paste0("axe",seq(1,ncols))
  }
  # Filtrer entre -1 et 1
  dt_uses_env2 = dt_uses_env %>% 
    filter( axe1 >= limits[1] ) %>%
    filter( axe2 >= limits[3] ) %>%
    filter( axe1 <= limits[2] ) %>%
    filter( axe2 <= limits[4] ) 
  
  # Grid 100 * 100
  dt_uses_env_grid2 = dt_uses_env2 %>% mutate(
    cut_x = cut(axe1, breaks = round(seq(from = limits[1], to = limits[2], length.out = 100),2),
                include.lowest = T),
    cut_y = cut(axe2, breaks = round(seq(from = limits[3], to = limits[4], length.out = 100),2),
                include.lowest = T)
  ) %>%
    group_by(cut_x, cut_y, .drop = FALSE) %>% 
    summarise(n_bin = n(), 
              mean_proba = mean(proba_presence, na.rm=T ),
              sd_proba = sd(proba_presence),
              sum_obs = sum(obs_usage))
  
  # correction des densités par la disponibilité
  MAX_N = max(dt_uses_env_grid2$n_bin,na.rm=T)
  MAX_nobs = max(dt_uses_env_grid2$sum_obs, na.rm=T)
  dt_uses_env_grid2$e_ij <- dt_uses_env_grid2$n_bin / MAX_N
  dt_uses_env_grid2$e_ij[dt_uses_env_grid2$mean_proba == "NaN"] <- NA
  dt_uses_env_grid2$o_ij <- dt_uses_env_grid2$sum_obs / MAX_nobs
  dt_uses_env_grid2$o_ij[dt_uses_env_grid2$mean_proba == "NaN"] <- NA
  max_oe = max(dt_uses_env_grid2$o_ij/dt_uses_env_grid2$e_ij, na.rm=T)
  dt_uses_env_grid2$z_ij <- (dt_uses_env_grid2$o_ij / dt_uses_env_grid2$e_ij) / max_oe
  
  # # AVEC GEOM_RASTER
  # # available environment
  # Pe = dt_grid_test %>% 
  #   ggplot(aes(cut_x, cut_y, fill=e_ij, alpha = e_ij)) +
  #   geom_raster() +
  #   scale_fill_distiller(palette ="Spectral")
  # # density observation
  # Po = dt_grid_test %>% 
  #   ggplot(aes(cut_x, cut_y, fill=o_ij, alpha=o_ij)) +
  #   geom_raster() +
  #   scale_fill_distiller(palette ="Spectral")
  # # densité obs en fonction disponibilité env
  # Poe = dt_grid_test %>% 
  #   ggplot(aes(cut_x, cut_y, fill=z_ij, alpha=z_ij)) +
  #   geom_raster() +
  #   scale_fill_distiller(palette ="Spectral")
  
  # avec des points (à la place de carrés)
  xlabels <- levels(dt_uses_env_grid2$cut_x)
  ind_remove = c(seq(2, length(xlabels)/2, 1), seq(round(length(xlabels)/2)+1, length(xlabels)-1, 1))
  xlabels[ind_remove ] <- ""
  
  P1 =  dt_uses_env_grid2 %>% 
    ggplot(aes(cut_x, cut_y, colour=e_ij)) +
    geom_point(size=2) +
    scale_colour_distiller(palette ="Spectral",na.value = "transparent") +
    theme(#axis.ticks.x = element_blank(),
      legend.position = "none") +
    labs(x="Axe 1",y="Axe 2",  title = "Environment Availibity")+
    theme(axis.text.x= element_text(angle = 45, hjust = 1))+
    scale_x_discrete(labels = xlabels)+
    scale_y_discrete(labels = xlabels)
  
  P2 =  dt_uses_env_grid2 %>% 
    ggplot(aes(cut_x, cut_y, colour = ifelse(o_ij > 0 , o_ij, NA))) +
    geom_point( size=2) +
    scale_colour_distiller(palette ="Spectral",na.value = "transparent") +
    theme(#axis.ticks.x = element_blank(),
      legend.position = "none") +
    labs(x="Axe 1",y="Axe 2",  title = "Density of Observations")+
    theme(axis.text.x= element_text(angle = 45, hjust = 1))+
    scale_x_discrete(labels = xlabels)+
    scale_y_discrete(labels = xlabels)
  
  # bugs : z = 1 ? pour des pixels où peu d'obs + env rare
  # quand pixel env le plus rare avec 1 seule observation de présence -> z = 1
  
  # z = 0 lorsque obs < 2
  stock = dt_uses_env_grid2
  nb_corr = length(dt_uses_env_grid2$z_ij[dt_uses_env_grid2$sum_obs < 2 & dt_uses_env_grid2$z_ij == 1])
  cat(paste0("\nz corrigé pour ",nb_corr," cellules."))
  dt_uses_env_grid2$z_ij[dt_uses_env_grid2$sum_obs < 2 & dt_uses_env_grid2$z_ij == 1] <- 0
  
  # avec des points (à la place de carrés)
  P3 = dt_uses_env_grid2 %>% 
    ggplot(aes(cut_x, cut_y, colour=ifelse(z_ij > 0 , z_ij, NA))) +
    geom_point(size=2) +
    scale_colour_distiller(palette ="Spectral",na.value = "transparent") +
    labs(x="Axe 1",y="Axe 2",  title="Occupancy",colour = "Corrected Density\nof Observations")+
    theme(axis.text.x= element_text(angle = 45, hjust = 1))+
    scale_x_discrete(labels = xlabels)+
    scale_y_discrete(labels = xlabels)
  
  P_all = P2 + P1 + P3+ plot_layout(nrow = 2, ncol=2) 
  
  png(file = paste0(chemin_esp_eco,"/occupancy_obs_",mois,".png"),width=1400, height=800)
  plot(P_all)
  dev.off()
  
  # densité condition env pour les présences observées + contour niche 
  dt_uses_env_grid = dt_uses_env %>% mutate(
    cut_x = cut(axe1, breaks = seq(from = limits[1], to = limits[2], length.out = 100), include.lowest = T),
    cut_y = cut(axe2, breaks = seq(from = limits[3], to = limits[4], length.out = 100), include.lowest = T)
  ) %>%
    group_by(cut_x, cut_y) %>%
    mutate(n_bin = n(),
           mean_proba = mean(proba_presence),
           sd_proba = sd(proba_presence))
  grid_usage2 = grid_usage_rdata %>% filter(proba_presence >= seuil)
  grid_usage2 = grid_usage_rdata
  dt_uses_env_grid$obs_usage[dt_uses_env_grid$obs_usage == 0.5] <- 1
  # New facet label names
  supp.labs <- c("Absence", "Presence")
  names(supp.labs) <- c(0, 1)
  P =  dt_uses_env_grid %>% 
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
    theme(panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey")) +
    guides(alpha = FALSE)+
    theme(text = element_text(size=15)) +
    facet_grid(~ obs_usage,labeller = labeller(obs_usage = supp.labs))
  
  #save plot
  png(file = paste0(chemin_esp_eco,"/niche_contour_obs_",mois,".png"),width=1400, height=800)
  plot(P)
  dev.off()
  
  # rescale occupancy
  dt_uses_env_grid2$p <- dt_uses_env_grid2$z_ij/sum(dt_uses_env_grid2$z_ij, na.rm=T)
  
  return(dt_uses_env_grid2)
}

# Schoener D sur occupancy à partir données observations
applyD_Schoener <- function(liste.usages, mois){
  
  # # TEST
  # liste.usages = c("Ni","Pa","Rp","Co","Vt")
  # mois = "juillet"
  
  # compute occupancy fpr each use
  list_dt <- lapply(liste.usages, 
                    function(x) GridObs(usage = x, mois = mois))
  # sélectionner colonne p de chaque dt et index la (same as usage order)
  df = data.frame()
  for(i in 1:length(list_dt)){
    df.u = data.frame("p" = list_dt[[i]]$p)
    names(df.u) = paste0("p_",liste.usages[i])
    df = as.data.frame(append(df, df.u))
  }
  # Schoener D pairwise computation matrix
  pairwise_D <- matrix(nrow = length(liste.usages), ncol = length(liste.usages), 
                       dimnames = list(liste.usages,liste.usages) )
  for(i in 1:length(liste.usages)){
    for(j in 1:length(liste.usages)){
      pairwise_D[j,i]  <- DSchoener(df[,i], df[,j])
    }
  }
  return(pairwise_D)
}

# Calcule moyenne et écart types de plusieurs matrices
meansd4listMatrices <- function(liste.matrices, liste.usages){
  # # TEST
  # liste.matrices = A
  # liste.usages = liste.usages
  
  pairwise_mean_D <- matrix(nrow = length(liste.usages), ncol = length(liste.usages), 
                            dimnames = list(liste.usages,liste.usages) )
  pairwise_sd_D <- matrix(nrow = length(liste.usages), ncol = length(liste.usages), 
                          dimnames = list(liste.usages,liste.usages) )
  for(i in 1:length(liste.usages)){
    for(j in 1:length(liste.usages)){
      pairwise_mean_D[j,i]  <- mean(unlist(lapply(liste.matrices, function(x) x[j,i])))
      pairwise_sd_D[j,i] <- sd(unlist(lapply(liste.matrices, function(x) x[j,i])))
    }
  }
  pairwise_smr <- list(pairwise_mean_D, pairwise_sd_D)
  names(pairwise_smr) <- c("mean","sd")
  return(pairwise_smr)
}

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

A = lapply(liste.mois[-1], function(x) 
  applyD_Schoener(liste.usages = c("Ni","Pa","Rp","Co","Vt"), mois = x) )
# Rp et VT présents tous l'été
Abis = lapply(liste.mois, function(x) 
  applyD_Schoener(liste.usages = c("Rp","Vt"), mois = x) )
# Lk que en mai
B = applyD_Schoener(liste.usages = c("Lk","Rp","Vt"), mois = "mai")
schoenerD_list = append(list(B),A)
names(schoenerD_list) <- c("mai", liste.mois[-1])
# save matrices
type_donnees = "ACP_avec_ponderation"
algorithme = "glm"
fit = "2_axes"
schoener_D_obs_path = paste0(output_path,"/niches/",type_donnees,"/niche_overlap/",
                            fit,"/",algorithme,"/Schoener_D_obs/") 
if(!dir.exists(schoener_D_obs_path)){dir.create(schoener_D_obs_path)}
# save in rdata
save(schoenerD_list,
     file = paste0(schoener_D_obs_path,"/matrices_schoener_d_obs.rdata"))
# save in csv
for(i in 1:length(schoenerD_list)){
  write.csv(schoenerD_list[i], 
            paste0(schoener_D_obs_path,"/matrice_schoener_d_obs_",names(schoenerD_list)[i],".csv"))
}
# mean + sd throught summer
M1 = meansd4listMatrices(A, c("Ni","Pa","Rp","Co","Vt"))
M2 = meansd4listMatrices(Abis, c("Rp","Vt"))

M2[1,2]

M3 = B 
M3[,2:3] <- NA
M3

# TODO : combine means & sd in the same matrix
# de la forme mean (+- sd) dans chaque cellule







#####  Schoener D sur probabilités issues SDM, projetées dans esp env ##### 

#####  % surface niche binarisé (en fonction cutoff) ##### 

#### G-space ####

##### % surface obs #####
# ("témoin", ce qui est basiquement fait)

##### Schoener D proba occ ##### 
# projetées dans esp géographique

#####  % surface pred occ ##### 
# (en fonction cutoff)
