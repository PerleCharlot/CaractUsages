### Titre -------------------------------------
# Nom : Calcul indices overlap
# Auteure : Perle Charlot
# Date de création : 10-03-2023
# Dates de modification : 17-03-2023

### Librairies -------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(raster)
library(sf)
library(stringr)
library(tidyverse)
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
  
  cat("\nUsage : ", usage)
  
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
  
  cat("\nMois : ", mois)
  
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
      pairwise_mean_D[j,i]  <- mean(unlist(lapply(liste.matrices, function(x) x[j,i])), na.rm=T)
      pairwise_sd_D[j,i] <- sd(unlist(lapply(liste.matrices, function(x) x[j,i])), na.rm=T)
    }
  }
  pairwise_smr <- list(pairwise_mean_D, pairwise_sd_D)
  names(pairwise_smr) <- c("mean","sd")
  return(pairwise_smr)
}

schoenerD.stats <- function(fonction_applyschoener, chemin_save){
  # # TEST
  # fonction_applyschoener = "G_proba" # "E_obs" "E_proba"
  # chemin_save = path_save
  
  if(fonction_applyschoener == "G_proba"){
    # en mai : que 3 usages
    D_mai <- applyD_Schoener_proba_pred(liste_usages = sort(c("Lk","Rp","Vt")), mois = c("mai"),space="G")
    # en juin : tous les usages
    D_juin <- applyD_Schoener_proba_pred(liste_usages = liste.usages, mois = c("juin"),space="G")
    # de juillet à septembre : les 5 mêmes usages
    D_jui_ao_sep <- lapply(c("juillet","aout","septembre"), function(x) 
      applyD_Schoener_proba_pred(liste_usages = sort(c("Ni","Pa","Rp","Co","Vt")), mois = x,space="G"))
  }
  if(fonction_applyschoener == "E_obs"){

    # Appliquer fonction en fonction des usages en présence
    # ! nécessite de le savoir en amont ...
    # en mai : que 3 usages
    D_mai <- applyD_Schoener_obs(liste_usages = sort(c("Lk","Rp","Vt")), mois = c("mai"))
    # en juin : tous les usages
    D_juin <- applyD_Schoener_obs(liste_usages = sort(liste.usages), mois = c("juin"))
    # de juillet à septembre : les 5 mêmes usages
    D_jui_ao_sep <- lapply(c("juillet","aout","septembre"), function(x) 
      applyD_Schoener_obs(liste_usages = sort(c("Ni","Pa","Rp","Co","Vt")), mois = x) )
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
    
    # Rassembler les matrices au cours de l'été
    D_summer <- append(list(D_mai, D_juin), D_jui_ao_sep)
    # Combler les vides quand absence usage pendant un mois
    a = unlist(lapply(1:length(D_summer ), function(x) dim(D_summer[[x]])[1]))
    higher_M <- D_summer[[which.max(a)]]
    D_summer2 <- lapply(1:length(D_summer), 
                        function(x) MatchMatrixDims(biggest_mat = higher_M, 
                                                    mat_to_expand =  D_summer[[x]]))
    # remove upper triangle
    D_summer2 <- lapply(1:length(D_summer2), function(x) {
      M <- D_summer2[[x]]
      M[upper.tri(M)] <- NA
      return(M)})
    names(D_summer2) <- liste.mois
    
    # save matrices
    save(D_summer2,
         file = paste0(chemin_save,"/matrices_schoener_d.rdata"))
    # save in csv
    for(i in 1:length(D_summer2)){
      write.csv(D_summer2[i], 
                paste0(chemin_save,"/matrice_schoener_d_",names(D_summer2)[i],".csv"))
    }
    
    # mean + sd throught summer
    M_mean_sd <- meansd4listMatrices(D_summer2, sort(liste.usages))
    for(i in 1:length(M_mean_sd)){
      write.csv(M_mean_sd[i], 
                paste0(chemin_save,"/schoener_d_",names(M_mean_sd)[i],".csv"))
    }
    
    M = matrix(paste0(round(M_mean_sd$mean, 2), " (", 
                      round(M_mean_sd$sd, 2), ")"),
               6,6,
               dimnames = list(sort(liste.usages),sort(liste.usages)))
    M[upper.tri(M)] <- NA
    diag(M) <- 1
    df_mean_sd_D_obs = as.data.frame(M)
    df_mean_sd_D_obs
    write.csv(df_mean_sd_D_obs, 
              paste0(chemin_save,"/mat_schoener_d_mean_sd_summer.csv"))
  }
}

# Extand a matrix to match a biffer one
MatchMatrixDims <- function(biggest_mat, mat_to_expand){
  # # TEST
  # biggest_mat = higher_M
  # mat_to_expand = M_inter_obs[[1]]
  
  df_sup <- as.data.frame(biggest_mat)
  df <- as.data.frame(mat_to_expand)
  
  # quand c'est la matrice la plus grande vs elle même
  if(sum(dim(df_sup) == dim(df)) == 2){
    return(as.matrix(df))
  }else{
    # extract absent row.s of df from df_sup
    abs.rows <- df_sup[!rownames(df_sup) %in% rownames(df),]
    # remove abs col.s
    abs.rows <- abs.rows[,colnames(abs.rows) %in% colnames(df)]
    abs.rows[] <- NA
    df.2 <- rbind(df, abs.rows)
    # extract absent col.s of df from df_sup
    abs.cols <- as.data.frame(df_sup[,!colnames(df_sup) %in% colnames(df)])
    names(abs.cols) <- names(df_sup)[!colnames(df_sup) %in% colnames(df)]
    rownames(abs.cols) <- rownames(df_sup)
    abs.cols[] <- NA
    df <- cbind(df.2, abs.cols)
    # reorder
    df.sorted <- df[order(rownames(df)),]
    df <- df[match(rownames(df.sorted ), rownames(df)),]
    df <- df[,order(colnames(df))]
    
    return(as.matrix(df))
  }
}

# Fonction qui calcule l'intersection entre obs usages et sort matrice
inter_mois <- function(mois){
  # # TEST
  # mois = "juin"
  
  list_r_us = list.files(paste0(output_path,"/par_periode/",mois),
                         recursive=TRUE, ".tif$", full.names=TRUE)
  stack_us = stack(list_r_us)
  
  # masquer par contour zone d'étude car usages VTT et rando ont une emprise plus large
  stack_us <- raster::mask(stack_us, limiteN2000)
  
  surf_1_pix <- res(stack_us)[1]*res(stack_us)[2] # en m²
  
  M <- matrix(nrow =  dim(stack_us)[3],ncol = dim(stack_us)[3],
              dimnames = list(names(stack_us),names(stack_us)))
  
  for(i in 1:dim(stack_us)[3]){
    for(j in 1:dim(stack_us)[3]){
      nb_pix_u <- sum(values(stack_us[[i]]), na.rm=T) # nb pixels d'un usage
      A_u <- nb_pix_u * surf_1_pix # en m²
      pair <- stack_us[[i]]+stack_us[[j]]
      A_inter <- length(pair[pair == 2]) * surf_1_pix # en m²
      C_u <- A_inter/A_u 
      name_i <- names(stack_us[[i]])
      name_j <- names(stack_us[[j]])
      cat(paste0("\nC(",name_i,"-",name_j,") : ", C_u))
      M[i,j] <- C_u
    }
  }
  return(M)
}

# Calcul mean et sd sur une liste de matrice mensuelle contenant ratio surface
MeanSd_Surface_Month <- function(liste_de_matrices_mensuelle,
                                 chemin_save){
  # # TEST
  # liste_de_matrices_mensuelle = M_inter_obs
  # chemin_save = path_save
  
  a = unlist(lapply(1:length(liste_de_matrices_mensuelle ), 
                    function(x) dim(liste_de_matrices_mensuelle[[x]])[1]))
  higher_M <- liste_de_matrices_mensuelle[[which.max(a)]]
  liste_de_matrices_mensuelle <- lapply(1:length(liste_de_matrices_mensuelle), 
                                        function(x) MatchMatrixDims(biggest_mat = higher_M, 
                                                                    mat_to_expand =  liste_de_matrices_mensuelle[[x]]))
  names(liste_de_matrices_mensuelle) <- liste.mois
  for(i in 1:length(liste_de_matrices_mensuelle)){
    write.csv(liste_de_matrices_mensuelle[i], 
              paste0(path_save,"/mat_surf_",names(liste_de_matrices_mensuelle)[i],".csv"))
  }
  # compute mean + sd
  
  M_mean_sd <-meansd4listMatrices(liste_de_matrices_mensuelle, 
                                  sort(colnames(liste_de_matrices_mensuelle[[1]])))
  for(i in 1:length(M_mean_sd)){
    write.csv(M_mean_sd[i], 
              paste0(chemin_save,"/ratio_surf_",names(M_mean_sd)[i],".csv"))
  }
  
  M = matrix(paste0(round(M_mean_sd$mean, 2), " (", 
                    round(M_mean_sd$sd, 2), ")"),
             6,6,
             dimnames = list(sort(colnames(liste_de_matrices_mensuelle[[1]])),
                             sort(colnames(liste_de_matrices_mensuelle[[1]]))))
  diag(M) <- 1
  df_mean_sd_D_obs = as.data.frame(M)
  df_mean_sd_D_obs
  write.csv(df_mean_sd_D_obs, 
            paste0(chemin_save,"/mat_ratio_surf_mean_sd_summer.csv"))
  
}

### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")
input_path <- paste0(wd,"/input/")

#### Données spatiales ####
dos_var_sp <- "C:/Users/perle.charlot/Documents/PhD/DATA/Variables_spatiales_Belledonne/"
limiteN2000 <- paste0(dos_var_sp, "/limites_etude/cembraie_N2000_limites.gpkg")

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
limiteN2000 <- st_read(limiteN2000)

M_inter_obs <- lapply(liste.mois,inter_mois)
MeanSd_Surface_Month(M_inter_obs, path_save)

##### Schoener D proba occ ##### 
# projetées dans esp géographique
path_save <- paste0(schoener_D_Gspace_path,"/proba/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}
schoenerD.stats("G_proba", path_save)

#####  % surface pred occ ##### 
# (en fonction cutoff)

path_save <- paste0(schoener_D_Gspace_path,"/surface_pred/")
if(!dir.exists(path_save)){dir.create(path_save, recursive = T)}

#% superficie intersection en fonction cutoff
pInter_fCutOff <- function(mois, space){
  
  # # TEST
  # mois = "juin" # NULL "juin"
  # space = "G" # "G"
  
  # sortir usage nom
  noms_us_ord <- liste.usages[order(liste.usages)]
  
  if(space == "G"){
    # load raster of proba for month studied
    list_rast <- list.files(path = paste0(output_path,"/niches/", type_donnees,"/"),
                                  pattern = ".tif$", 
                                  recursive = T, full.names = T)
    # trier le mois
    list_rast <- list_rast[grep(mois,list_rast)]
    # obtenir noms usages
    A <- list.files(path = paste0(output_path,"/niches/", type_donnees,"/"),
               pattern = ".tif$", 
               recursive = T)
    A  <- A [grep(mois,A)]
    noms_us_ord <- str_sub(A,start = 1,end = 2)

    # Load rasters
    stack_us <- stack(list_rast[grep(mois,list_rast)])
    # rename with use names
    names(stack_us) <- unlist(lapply(noms_us_ord, function(x) paste0(x,c("_obs","_proba","_pred"))))
    # conserve proba layer = 2nd
    stack_proba <- stack_us %>% subset(grep("proba", names(stack_us)))
    # transform to df
    df_proba_us <- data.frame(data.table(stack_proba[]))
    
    # creer un df par cutoff puis liste de df
    convertCutoff <- function(cutoff, index_col){
      b <- data.frame(b = ifelse(df_proba_us[,index_col] > cutoff, 1 ,0))
      names(b) <- paste0(names(df_proba_us[index_col]),"_",cutoff)
      return(b)
    }
    liste_preds <- lapply(seq(0.05,0.95,0.05), function(x) 
      do.call(cbind,
              lapply(1:dim(df_proba_us)[2], 
                     function(i) convertCutoff(cutoff=x, index_col = i))))
    
    # une matrice de calcul (de % aire intersection / aire usage ) / cutoff
    surf_1_pix <- res(stack_us)[1]*res(stack_us)[2] # en m²
    
    inter_cutoff <- function(df){
      # # TEST
      # df = liste_preds[[11]]
      # cat(names(df))
      
      M <- matrix(nrow =  length(noms_us_ord),ncol = length(noms_us_ord),
                  dimnames = list(noms_us_ord,noms_us_ord))
      
      for(i in 1:dim(df)[2]){
        for(j in 1:dim(df)[2]){
          
          nb_pix_u <- sum(df[,i], na.rm=T) # nb pixels d'un usage
          A_u <- nb_pix_u * surf_1_pix # en m²
          #cat(paste0("\nA_u : ", A_u/1000000))
          tb <- table(data.frame(inter = df[,i] + df[,j]))
          A_inter <- tb["2"] * surf_1_pix # en m²
          # si 0 pixels d'intersection => chevauchement nul
          if(is.na(A_inter)){
            A_inter <- 0
            C_u <- 0}
          # si usage jamais prédit => NA
          if(A_u == 0){
            C_u <- NA
          }
          if(A_inter != 0 & A_u != 0){
            C_u <- A_inter/A_u
          }
          #cat(paste0("\n C_u : ",C_u))
          name_i <- names(df)[i]
          name_j <- names(df)[j]
          #cat(paste0("\nC(",name_i,"-",name_j,") : ", C_u))
          M[i,j] <- C_u
        }
      }
      return(M)
    }
    
    liste_M <- lapply(liste_preds,inter_cutoff)
    names(liste_M) <- paste0("Ainter_",seq(0.05,0.95,0.05))
    
    out_mois <- paste0(path_save,"/",mois,"/") 
    if(!dir.exists(out_mois)){dir.create(out_mois, recursive = T)}
    
    for(i in 1:length(liste_M)){
      write.csv(liste_M[i], 
                paste0(out_mois,"/mat_surf_",names(liste_M)[i],".csv"))
    }
    
    # TODO : pas tous les usages présents au mois donné !!
    # f2 doit prendre ça en considération
    
    # Graphs variation overla^p en fonction cutoff
    f2 <- function(USE){
      cat(paste0("\nUse : ", USE, " pour mois ",mois ))
      f <- function(i,us){
        M = liste_M[[i]]
        return(M[us,])
      }
      # USE <- "Ni"
      df <- as.data.frame(do.call(rbind,lapply(1:length(liste_M), function(x)f(i=x, us=USE))))
      df$cutoff <- seq(0.05,0.95,0.05)
      
      df2 <- df %>% pivot_longer(liste.usages,
                                 names_to = "use2",
                                 values_to = "C")
      P = df2 %>% ggplot(aes(x=cutoff, y=C, col=use2))+
        geom_point()+geom_line()+
        labs(title=USE,y="Overlap")+
        theme(text = element_text(size=15))
      
      png(file = paste0(out_mois,"/overlap_surface_cutoff_",USE,"_",mois,".png"),width=1400, height=800)
      plot(P)
      dev.off()
    }
    lapply(liste.usages, f2)
    
    return(liste_M) # pour un mois donné, avec tous les cutoffs possibles
  }
  
  if(space == "E"){
    
  }
}

# liste de 6 éléments (1/mois), overlap=f(cutoff)
overlap_pred_G <- lapply(liste.mois,function(x) pInter_fCutOff(mois = x, space="G"))

# TODO : mean + sd ?