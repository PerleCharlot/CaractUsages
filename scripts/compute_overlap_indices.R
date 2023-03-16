### Titre -------------------------------------
# Nom : Calcul indices overlap
# Auteure : Perle Charlot
# Date de création : 10-03-2023
# Dates de modification : 16-03-2023

### Librairies -------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(raster)
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
applyD_Schoener_obs <- function(liste_usages, mois){
  
  # # TEST
  # liste.usages = c("Ni","Pa","Rp","Co","Vt")
  # mois = "juillet"
  
  # compute occupancy fpr each use
  list_dt <- lapply(liste_usages, 
                    function(x) GridObs(usage = x, mois = mois))
  # sélectionner colonne p de chaque dt et index la (same as usage order)
  df = data.frame()
  for(i in 1:length(list_dt)){
    df.u = data.frame("p" = list_dt[[i]]$p)
    names(df.u) = paste0("p_",liste_usages[i])
    df = as.data.frame(append(df, df.u))
  }
  # Schoener D pairwise computation matrix
  pairwise_D <- matrix(nrow = length(liste_usages), ncol = length(liste_usages), 
                       dimnames = list(liste_usages,liste_usages) )
  for(i in 1:length(liste_usages)){
    for(j in 1:length(liste_usages)){
      pairwise_D[j,i]  <- DSchoener(df[,i], df[,j])
    }
  }
  return(pairwise_D)
}

# Schoener D sur occupancy à partir données observations
applyD_Schoener_proba_pred <- function(liste_usages, 
                                       mois,
                                       space){

  # # TEST
  # liste_usages = liste.usages
  # mois = NULL # NULL "juin"
  # space = "E" # "G"
  
  # sortir usage nom
  noms_us_ord <- liste_usages[order(liste_usages)]
  
  if(space == "G"){
    # load raster of proba for month studied
    list_rast <- base::list.files(path = paste0(output_path,"/niches/", type_donnees,"/"),
                                  pattern = ".tif$", 
                                  recursive = T, full.names = T)
    # trier le mois
    list_rast <- list_rast[grep(mois,list_rast)]
    #trier les usages
    list_rast <- list_rast[grep(paste(liste_usages, collapse="|"), list_rast)]
    # Load rasters
    stack_us = stack(list_rast[grep(mois,list_rast)])
    # rename with use names
    names(stack_us) <- unlist(lapply(noms_us_ord, function(x) paste0(x,c("_obs","_proba","_pred"))))
    # conserve proba layer = 2nd
    stack_proba <- stack_us %>% subset(grep("proba", names(stack_us)))
    # transform to df
    df_proba_us <- data.frame(data.table(stack_proba[]))
    
    

  }
  if(space == "E"){
    
    # load raster of proba for month studied
    list_rast <- base::list.files(path = paste0(output_path,"/niches/", type_donnees,"/"),
                                  pattern = "niche_potentielle", 
                                  recursive = T, full.names = T)
    list_rast <- list_rast[grep(".rdata",list_rast)]
    # fonction load niche data for a use
    loadProbaNiche <- function(i){
      load(list_rast[i])
      names(grid_usage_rdata)[3] <- paste0(names(grid_usage_rdata)[3],"_",noms_us_ord[i])
      return(grid_usage_rdata)
    }
    df_proba_us <- lapply(1:length(liste_usages), loadProbaNiche)
    # Merge
    df_proba_us <- Reduce(function(x, y) merge(x, y, all=FALSE), df_proba_us)
    df_proba_us <- df_proba_us %>% subset(select=-grep("axe", names(df_proba_us)))
  }

  # rescale chaque proba (pour que ça fasse une distrib de proba => sum = 1)
  df_proba_us_scale <- as.data.frame(apply(df_proba_us,
                                           2,
                                           function(x) x/sum(x, na.rm=T)))
  #apply(df_proba_us_scale, 2, function(x) sum(x, na.rm=T))
  
  # Schoener D pairwise computation matrix
  pairwise_D <- matrix(nrow = length(liste_usages), ncol = length(liste_usages), 
                       dimnames = list(liste_usages,liste_usages) )
  for(u1 in 1:length(liste_usages)){
    for(u2 in 1:length(liste_usages)){
      pairwise_D[u1,u2]  <- DSchoener(df_proba_us_scale[,u1], df_proba_us_scale[,u2])
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


schoenerD.stats <- function(fonction_applyschoener, chemin_save){
  # # TEST
  # fonction_applyschoener = "E_proba" # "E_obs" "E_proba"
  # chemin_save = path_save
  
  if(fonction_applyschoener == "G_proba"){
    A = lapply(liste.mois[-1], function(x) 
      applyD_Schoener_proba_pred(liste_usages = c("Ni","Pa","Rp","Co","Vt"), mois = x, space="G") )
    # Rp et VT présents tout l'été
    Abis = lapply(liste.mois, function(x) 
      applyD_Schoener_proba_pred(liste_usages = c("Rp","Vt"), mois = x) )
    # Lk que en mai
    B = applyD_Schoener_proba_pred(liste_usages = c("Lk","Rp","Vt"), mois = "mai")
  }
  if(fonction_applyschoener == "E_obs"){
    A = lapply(liste.mois[-1], function(x) 
      applyD_Schoener_obs(liste_usages = c("Ni","Pa","Rp","Co","Vt"), mois = x) )
    # Rp et VT présents tous l'été
    Abis = lapply(liste.mois, function(x) 
      applyD_Schoener_obs(liste_usages = c("Rp","Vt"), mois = x) )
    # Lk que en mai
    B = applyD_Schoener_obs(liste_usages = c("Lk","Rp","Vt"), mois = "mai")
  }
  # pas de mean ni sd car niches construites sur aggrégation mois
  if(fonction_applyschoener == "E_proba"){
      D = applyD_Schoener_proba_pred(liste_usages = liste.usages, mois=NULL, space= "E")
      D[upper.tri(D)] <- NA
      diag(D) <- 1
      df_D = as.data.frame(D)
      write.csv(df_D, 
                paste0(chemin_save,"/mat_schoener_d_summer.csv"))
  }
  
  # quand D schoener mensuel, calcul mean + sd
  if(any(fonction_applyschoener == "E_obs" | fonction_applyschoener == "G_proba")){
    # Combine all 3
    schoenerD_list = append(list(B),append(A, Abis))
    names(schoenerD_list) <- c("mai", paste0(liste.mois[-1],"_5uses"),paste0(liste.mois,"_2uses"))
    # save matrices
    save(schoenerD_list,
         file = paste0(chemin_save,"/matrices_schoener_d.rdata"))
    # save in csv
    for(i in 1:length(schoenerD_list)){
      write.csv(schoenerD_list[i], 
                paste0(chemin_save,"/matrice_schoener_d_",names(schoenerD_list)[i],".csv"))
    }
    # mean + sd throught summer
    M1 = meansd4listMatrices(A, c("Ni","Pa","Rp","Co","Vt"))
    names(M1) = paste0(names(M1),"_5uses")
    M2 = meansd4listMatrices(Abis, c("Rp","Vt"))
    names(M2) = paste0(names(M2),"_2uses")
    M12 = append(M1, M2)
    for(i in 1:length(M12)){
      write.csv(M12[i], 
                paste0(chemin_save,"/schoener_d_",names(M12)[i],".csv"))
    }
    M3 = B 
    # Combiner à la main means & sd in the same matrix
    M_mean = cbind(M1$mean,matrix(rep(NA,5), dimnames =list(c("Ni","Pa","Rp","Co","Vt"),"Lk")) )
    N = matrix(c(NA,NA,M3[2,1],NA,M3[3,1],1), dimnames = list(c("Ni","Pa","Rp","Co","Vt","Lk"),"Lk"))
    M_mean = rbind(M_mean, t(N))
    # remplacer la valeur pour paire Rp/Vt
    M_mean[5,3] <- M2$mean[2,1]
    # NA dans triangle du haut
    M_mean[upper.tri(M_mean)] <- NA
    
    M_sd = cbind(M1$sd,matrix(rep(NA,5), dimnames =list(c("Ni","Pa","Rp","Co","Vt"),"Lk")) )
    Nsd = matrix(rep(NA,6), dimnames = list(c("Ni","Pa","Rp","Co","Vt","Lk"),"Lk"))
    M_sd = rbind(M_sd, t(Nsd))
    # remplacer la valeur pour paire Rp/Vt
    M_sd[5,3] <- M2$sd[2,1]
    # NA dans triangle du haut
    M_sd[upper.tri(M_sd)] <- NA
    # de la forme mean (+- sd) dans chaque cellule
    M = matrix(paste0(round(M_mean, 2), " (", 
                      round(M_sd, 2), ")"),
               6,6,
               dimnames = list(c("Ni","Pa","Rp","Co","Vt","Lk"),c("Ni","Pa","Rp","Co","Vt","Lk")))
    M[upper.tri(M)] <- NA
    diag(M) <- 1
    M[6,c(1,2,4)] <- NA
    df_mean_sd_D_obs = as.data.frame(M)
    df_mean_sd_D_obs
    write.csv(df_mean_sd_D_obs, 
              paste0(chemin_save,"/mat_schoener_d_mean_sd_summer.csv"))
  }
  
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

type_donnees = "ACP_avec_ponderation"
algorithme = "glm"
fit = "2_axes"

### Programme -------------------------------------

#### E-space ####
schoener_D_Espace_path = paste0(output_path,"/niches/",type_donnees,"/niche_overlap/",
                                fit,"/",algorithme,"/Schoener_D/E_space") 
if(!dir.exists(schoener_D_Espace_path)){dir.create(schoener_D_Espace_path, recursive = T)}

##### Schoener D abond obs #####
# corrigées par disponibilité (Broennimann 2012), env grid
path_save <- paste0(schoener_D_Espace_path,"/obs/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}
schoenerD.stats("E_obs", path_save)

#####  Schoener D sur probabilités issues SDM, projetées dans esp env ##### 
path_save <- paste0(schoener_D_Espace_path,"/proba/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}
schoenerD.stats("E_proba", path_save)

#####  % surface niche binarisé (en fonction cutoff) ##### 

path_save <- paste0(schoener_D_Espace_path,"/surface_niche/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}

#### G-space ####
schoener_D_Gspace_path = paste0(output_path,"/niches/",type_donnees,"/niche_overlap/",
                                fit,"/",algorithme,"/Schoener_D/G_space") 
if(!dir.exists(schoener_D_Gspace_path)){dir.create(schoener_D_Gspace_path)}
##### % surface obs #####
# ("témoin", ce qui est basiquement fait)
path_save <- paste0(schoener_D_Gspace_path,"/surface_obs/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}

# monthly + stats

##### Schoener D proba occ ##### 
# projetées dans esp géographique
path_save <- paste0(schoener_D_Gspace_path,"/proba/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}
schoenerD.stats("G_proba", path_save)

#####  % surface pred occ ##### 
# (en fonction cutoff)

path_save <- paste0(schoener_D_Gspace_path,"/surface_pred/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}
