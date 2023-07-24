### Titre -------------------------------------
# Nom : Calcul indices overlap
# Auteure : Perle Charlot
# Date de création : 26-06-2023
# Dates de modification : 24-07-2023

### Librairies -------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(raster)
library(sf)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(hexbin)
library(scico)
### Fonctions -------------------------------------

# Fonction qui sort une table avec coordonnées xy et proba de tous les usages, pour un mois donné
GetGspaceTable <- function(mois,
                           usages_etudies = NULL, #si on veut une sélection d'usages particulière, 
                           # de base, la fonction fait pour tous les usages rencontrés le mois étudié
                           data_type,  # "proba" or "obs"
                           type_donnees_glm = "brute",
                           fit = "all_simple",
                           algorithme = "glm"){
  # # TEST
  # mois = "juin"
  # usages_etudies = NULL
  # data_type = "proba"
  
  # Create directories
  table_fold_path_month <- paste0(table_fold_path,mois,"/") 
  if(!dir.exists(table_fold_path_month)){dir.create(table_fold_path_month, recursive = T)}
  
  # Translate french <-> english (for figures)
  english_month <- df_time$english_month[df_time$mois == mois]
  
  # load raster of proba for month studied
  list_rast <- base::list.files(path = paste0(output_path,"/niches/", type_donnees_glm,"/"),
                                pattern = ".tif$", 
                                recursive = T, full.names = T)
  # Trier le threshold
  list_rast <- list_rast[grep("TSS",list_rast)]
  
  # #######
  # name_months <- df_time$english_month[order(df_time$english_month, df_time$mois)]
  # # ANALYSE
  # stack_Vt <- raster::stack(list_rast[grep("Vt",list_rast)])
  # names(stack_Vt) <-  paste0("Vt_",
  #        unlist(lapply(name_months, function(x) paste0(x,c("_obs","_proba","_pred")))))
  # stack_Vt_obs <- stack_Vt %>% subset(grep("obs", names(stack_Vt)))
  # stack_Vt_obs <- dropLayer(stack_Vt_obs,4)
  # df_Vt <- data.frame(data.table(stack_Vt_obs[]))
  # df_Vt <- cbind(coordinates(stack_Vt_obs),df_Vt)
  # names(df_Vt)[3:6] <- name_months[-4]
  # df_Vt <- df_Vt %>% pivot_longer(cols=3:6, names_to = "Month",
  #                                                  values_to = "Vt")
  # name_months2 <- name_months[-4]
  # stack_Co <- raster::stack(list_rast[grep("Co",list_rast)])
  # names(stack_Co) <-  paste0("Co_",
  #                            unlist(lapply(name_months2, function(x) paste0(x,c("_obs","_proba","_pred")))))
  # stack_Co_obs <- stack_Co %>% subset(grep("obs", names(stack_Co)))
  # df_Co <- data.frame(data.table(stack_Co_obs[]))
  # df_Co  <- cbind(coordinates(stack_Co_obs),df_Co)
  # names(df_Co)[3:6] <- name_months[-4]
  # df_Co <- df_Co %>% pivot_longer(cols=3:6,
  #                                 names_to = "Month",
  #                                 values_to = "Co")
  # df_Vt_Co <- merge(df_Vt, df_Co, by=c("x","y","Month"))
  # df_Vt_Co <-  df_Vt_Co %>% pivot_longer(cols=4:5,
  #                                 names_to = "Use",
  #                                 values_to = "Obs")
  # df_Vt_Co$Month <- factor(df_Vt_Co$Month,
  #                             levels = c("May","June","July","August","September"))
  # 
  # 
  # map_expts <- na.omit(df_Vt_Co) %>%
  #   ggplot() +
  #   geom_raster(aes(x = x, y = y, fill = as.factor(Obs))) +
  #   scale_fill_manual(values = c("1" = "#FDE725FF",
  #                                "0" = "#440154FF"),
  #                     labels=c('Absence', 'Presence'))+
  #   theme(panel.background = element_rect(fill="white"),
  #         plot.background = element_rect(fill="white"),
  #         panel.grid.major = element_line(colour="grey"),
  #         legend.position = "right",
  #         text = element_text(size=15),
  #         axis.text.x = element_text(angle=45)) +
  #   labs(fill="Probability of\noccurrence",
  #        x="Longitude",y="Latitude")+
  #   facet_grid(Month~Use,labeller = labeller(Use = uses.labs.obs),drop=F)+
  #   coord_equal()
  # png(file = paste0(table_fold_path,"/map_obs_Vt_Co.png"),width=1400, height=800)
  # plot(map_expts)
  # dev.off()
  # write.csv(df_Vt_Co, paste0(table_fold_path,"xy_obs_Vt_Co.csv"))
  # #######
  
  # trier le mois
  list_rast <- list_rast[grep(mois,list_rast)]
  
  #trier les usages
  if(!is.null(usages_etudies)){
    list_rast <- list_rast[grep(paste(liste_usages, collapse="|"), list_rast)]}
  # Load rasters
  stack_uses <- raster::stack(list_rast[grep(mois,list_rast)])
  # Récupérer les noms des usages présents au mois étudié + renommer stack
  n <- nchar(mois) + nchar(fit)+  nchar(algorithme) + 47
  uses_names <- str_sub(string = list_rast, start=-n-1, end=-n)
  names(stack_uses) <- unlist(lapply(uses_names, function(x) paste0(x,c("_obs","_proba","_pred"))))
  
  # conserve proba
  stack_uses2 <- stack_uses %>% subset(grep(data_type, names(stack_uses)))
  plot(stack_uses2)

  # raster -> df
  df_xy_uses <- data.frame(data.table(stack_uses2[]))
  df_xy_uses <- cbind(coordinates(stack_uses2),df_xy_uses)
  
  # pivot
  df_xy_uses <- df_xy_uses %>% pivot_longer(cols=ends_with(data_type),
                                                  names_to = "Use", 
                                                  values_to = data_type)
  
  # Remove "obs" or "proba" radical
  df_xy_uses$Use <- str_split(df_xy_uses$Use,'_',simplify = T)[,1]

  # Reordering group factor levels
  df_xy_uses$Use <- factor(df_xy_uses$Use,      
                              levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt")))
  df_xy_uses <- as.data.frame(df_xy_uses) 

  #names(uses.labs) <- paste0(names(uses.labs))
  
  # Plot map (G space) PROBA PRESENCE
  if(data_type == "proba"){
    map <- na.omit(df_xy_uses) %>%
      ggplot() +
      geom_raster(aes(x = x, y = y, fill = proba)) +
      scale_fill_viridis(limits=c(0,1)) +
      theme(panel.background = element_rect(fill="white"),
            plot.background = element_rect(fill="white"),
            panel.grid.major = element_line(colour="grey"),
            legend.position = "right",
            text = element_text(size=15),
            axis.text.x = element_text(angle=45)) +
      labs(fill="Probability of\noccurrence", title = english_month,
           x="Longitude",y="Latitude")+
      facet_wrap(.~Use,labeller = labeller(Use = uses.labs),drop=F)+
      coord_equal()
  }
  if(data_type == "obs"){
    map <- na.omit(df_xy_uses) %>%
      ggplot() +
      geom_raster(aes(x = x, y = y, fill = as.factor(obs))) +
      scale_fill_manual(values = c("1" = "#FDE725FF",
                                   "0" = "#440154FF"),
                        labels=c('Absence', 'Presence'))+
      theme(panel.background = element_rect(fill="white"),
            plot.background = element_rect(fill="white"),
            panel.grid.major = element_line(colour="grey"),
            legend.position = "right",
            text = element_text(size=15),
            axis.text.x = element_text(angle=45)) +
      labs(fill="Observed\noccurrence",title = english_month,
           x="Longitude",y="Latitude")+
      facet_wrap(.~Use,labeller = labeller(Use = uses.labs),drop=F)+
      coord_equal()
    
    df_xy_uses2 <- df_xy_uses
  }
  map
  #save plot
  png(file = paste0(table_fold_path_month,"/map_",data_type,"_uses.png"),width=1400, height=800)
  plot(map)
  dev.off()

  # if work with probabilities, need to make binary pres/abs with TSS threshold
  if(data_type == "proba"){
    #select threshold for each use
    df_xy_uses2 <- summary_models %>%
      dplyr::select(use,threshold_kept) %>%
      rename(Use = use) %>%
      mutate(Use = paste0(Use)) %>% 
      merge(df_xy_uses, by="Use",all=T)
    df_xy_uses2$pres_pred <- ifelse(df_xy_uses2$proba > df_xy_uses2$threshold,1,0)
    
    # Reordering group factor levels
    df_xy_uses2$Use <- factor(df_xy_uses2$Use,      
                              levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt")))
    
    # Plot map (G space) PRESENCE BINERISE (pred)
    map2 <- na.omit(df_xy_uses2) %>%
      ggplot() +
      geom_raster(aes(x = x, y = y, fill = as.factor(pres_pred))) +
      scale_fill_manual(values = c("1" = "#FDE725FF",
                                   "0" = "#440154FF"),
                        labels=c('Absence', 'Presence')) +
      theme(panel.background = element_rect(fill="white"),
            plot.background = element_rect(fill="white"),
            panel.grid.major = element_line(colour="grey"),
            legend.position = "right",
            text = element_text(size=15),
            axis.text.x = element_text(angle=45)) +
      labs(fill="Predicted of\noutcome", title = english_month,
           x="Longitude",y="Latitude")+
      facet_wrap(.~Use,labeller = labeller(Use = uses.labs),drop=F)+
      coord_equal()
    map2
    png(file = paste0(table_fold_path_month,"/map_pred_uses.png"),width=1400, height=800)
    plot(map2)
    dev.off()
  }

  # save df
  df_xy_uses2$Month <- english_month
  # make different name if not on all uses
  df_name <- paste0("/df_coords_",data_type,"_uses",paste0(usages_etudies,collapse="_"),".csv")
  write.csv(df_xy_uses2, paste0(table_fold_path_month,df_name))
  
  return(df_xy_uses2)
}



# Fonction qui sort une table avec coordonnées xy et valeurs PCA1 PCA2 pour ACP globale, pour un mois donné
GetEspaceTableGlob <- function(mois,
                               usages_etudies = NULL, #si on veut une sélection d'usages particulière, 
                               # de base, la fonction fait pour tous les usages rencontrés le mois étudié
                               type_donnees_glm = "brute",
                               fit = "all_simple",
                               algorithme = "glm"){ # or "summer"
  # # TEST
  # mois = "juin"
  # mois = NULL
  # usages_etudies = NULL
  
  # Create directories
  table_fold_path_month <- paste0(table_fold_path,mois,"/") 
  if(!dir.exists(table_fold_path_month)){dir.create(table_fold_path_month, recursive = T)}
  
  # Translate french <-> english (for figures)
  english_month <- df_time$english_month[df_time$mois == mois]
  # Get PCA axes from global analysis

  if(is.null(mois)){
    PCA1 <- stack(list.files(path=paste0(gitCaractMilieu,"/output/ACP/ACP_avec_ponderation/summer/"),
                             "axe1", full.names = T))
    
    PCA2 <- stack(list.files(path=paste0(gitCaractMilieu,"/output/ACP/ACP_avec_ponderation/summer/"),
                             "axe2", full.names = T))
  }else{
    PCA1 <- stack(list.files(path=paste0(gitCaractMilieu,"/output/ACP/ACP_avec_ponderation/summer/pred_month/",mois),
                             "axe1", full.names = T))
    
    PCA2 <- stack(list.files(path=paste0(gitCaractMilieu,"/output/ACP/ACP_avec_ponderation/summer/pred_month/",mois),
                             "axe2", full.names = T))
  }

  PCA_stack <- stack(PCA1,PCA2)
  names(PCA_stack) <- c("PCA1","PCA2")
  plot(PCA_stack)
  df_env <- data.frame(data.table(PCA_stack[]))
  df_env <- cbind(coordinates(PCA_stack),df_env)
  # Plot map (G space)
  df_env2 <- df_env %>% 
    pivot_longer(cols=starts_with("PCA"),values_to = "PCA_value" , names_to = "PCA_axe")
  map <- na.omit(df_env2) %>%
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = PCA_value)) +
    # scale_fill_viridis(limits=c(0,1)) +
    scale_fill_distiller(palette ="Spectral",direction=1)+
    theme(panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"),
          legend.position = "right",
          text = element_text(size=15),
          axis.text.x = element_text(angle=45)) +
    labs(fill="PCA value", title = english_month,
         x="Longitude",y="Latitude")+
    facet_wrap(.~PCA_axe)+
    coord_equal()
  map
  #save plot
  png(file = paste0(table_fold_path_month,"/map_PCA_axes.png"),width=1400, height=800)
  plot(map)
  dev.off()

  # save df
  if(is.null(mois)){}else{df_env$Month <- english_month}
  # make different name if not on all uses
  df_name <- "/df_coords_PCA_axes.csv"
  write.csv(df_env, paste0(table_fold_path_month,df_name))
  
  # df_env2$Month <- english_month
  # df_env2$dimension <- "global"
  
  return(df_env)
}
#  function to make a data.frame for geom_hex that can be used with stat_identity
# from
# https://stackoverflow.com/questions/39296198/operation-between-stat-summary-hex-plots-made-in-ggplot2
makeHexData <- function(df) {
  # # TEST
  # df <- A
  
  h <- hexbin(df$x, df$y, nbins, xbnds = xbnds, ybnds = ybnds, IDs = TRUE)
  data.frame(hcell2xy(h),
             z = tapply(df$z, h@cID, FUN = function(z) median(z)), # proba médiane
             n = tapply(df$z, h@cID, FUN = function(z) length(z)), # nb pixels géo
             # sd = tapply(df$z, h@cID, FUN = function(z) sd(z)),
             # mean = tapply(df$z, h@cID, FUN = function(z) mean(z)),
             # max = tapply(df$z, h@cID, FUN = function(z) max(z)),
             cid = h@cell)
}
# Compute hexagonal summary of proba for each use
f_hex <- function(df2, u){
  # # TEST
  # df2 <- df_for_hex
  # u <- 3
  
  # Subset pour l'usage u
  A <- df2[,c(1:2,u)]
  name_u <- names(A)[3]
  names(A)[3] <- "z"
  # compute grille hexagonale
  Ahex <- makeHexData(A)
  
  # rename and add vars
  if(exists("mois_run")){
    Ahex <- Ahex %>%
      rename(PCA1 =x, PCA2 = y) %>%
      mutate(Month = mois_run)
  }else{
    Ahex <- Ahex %>%
      rename(PCA1 =x, PCA2 = y)
  }
  
  # Filtrer par threshold
  th_u <- summary_models$threshold_kept[which(substr(name_u,1,2) %in% summary_models$use)] 
  Ahex$median_proba_pres <- ifelse(Ahex$z > th_u,
                                   1,0)
  names(Ahex)[3] <- name_u
  names(Ahex)[dim(Ahex)[2]] <- paste0(name_u,"_pred")
  return(Ahex)
}
# Fonction qui retourne une table, pour un mois donné, des valeurs env aggrégés
f_hex_month <- function(mois_run, df){
  # # TEST
  # mois_run <- "July"
  # df <- df_proba_PCA_xy
  # 
  # head(df)
  
  # Pour un mois donné, tableau avec tous les usages
  test_all <- df %>%
    filter(Month == mois_run) %>%
    dplyr::select(Use,PCA1,PCA2, Proba) %>%
    pivot_wider(names_from = Use, values_from = Proba) %>%
    rename(x=PCA1, y =PCA2)
  # # Compute hexagonal summary of proba for each use
  # f_hex <- function(df2, u){
  #   # # TEST
  #   # df2 <- df_for_hex
  #   # u <- 3
  #   
  #   # Subset pour l'usage u
  #   A <- df2[,c(1:2,u)]
  #   name_u <- names(A)[3]
  #   names(A)[3] <- "z"
  #   # compute grille hexagonale
  #   Ahex <- makeHexData(A)
  #   
  #   # rename and add vars
  #   if(exists("mois_run")){
  #     Ahex <- Ahex %>%
  #       rename(PCA1 =x, PCA2 = y) %>%
  #       mutate(Month = mois_run)
  #   }else{
  #     Ahex <- Ahex %>%
  #       rename(PCA1 =x, PCA2 = y)
  #   }
  #   
  #   # Filtrer par threshold
  #   th_u <- summary_models$threshold_kept[which(substr(name_u,1,2) %in% summary_models$use)] 
  #   Ahex$median_proba_pres <- ifelse(Ahex$z > th_u,
  #                                    1,0)
  #   names(Ahex)[3] <- name_u
  #   names(Ahex)[dim(Ahex)[2]] <- paste0(name_u,"_pred")
  #   return(Ahex)
  # }
  hex_uses <- lapply(3:dim(test_all)[2], function(x) f_hex(df2=test_all, u=x))
  
  hex_uses <- Reduce(function(x,y) merge(x = x, y = y, by = c('PCA1','PCA2','n','cid')), 
                     hex_uses)
  
  hex_uses <- hex_uses %>%
    pivot_longer(cols=ends_with("_proba"),
                 names_to = "Use" ,values_to = "median_proba") %>%
    pivot_longer(cols=ends_with("_proba_pred"),
                 names_to = "Use2" ,values_to = "pred")
  
  hex_uses$Month <- mois_run
  
  return(hex_uses)
}

# crée les figures :
# - G space : probabilité et prédiction binaires, 6 uses X 5 mois
# - E space : probabilité et prédiction binaires, 6 uses X 5 mois
Create.Plots.Monthly <- function(){
  # For all uses, across months, probabilities, in E space
  plot_allMonths_proba_E <- na.omit(df_proba_PCA_xy) %>%
    ggplot() +
    stat_summary_hex(aes(x=PCA1, y=PCA2, z= Proba),
                     fun = function(x) median(x), bins=40,colour='grey')+
    scale_fill_viridis(na.value = "transparent",
                       limits=c(0,1)) +
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed") +
    facet_grid(Month ~ Use, scales = "free",labeller = labeller(Use = uses.labs))+
    labs(title = "Probability of Occurrence projected in Ecological Space",
         fill = "Median probability\nof occurrence")+
    theme(text = element_text(size=18),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"))
  png(file = paste0(table_fold_path,"/proba_E_space_across_months_uses.png"),
      width=2100, height=1200)
  plot(plot_allMonths_proba_E)
  dev.off()
  # For all uses, across months, probabilities, in G space
  plot_allMonths_proba_G <- na.omit(df_proba_PCA_xy) %>%
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Proba)) +
    scale_fill_viridis(limits=c(0,1)) +
    theme(panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"),
          legend.position = "right",
          text = element_text(size=15),
          axis.text.x = element_text(angle=45)) +
    labs(fill="Probability of\noccurrence",
         x="Longitude",y="Latitude")+
    facet_grid(Month ~ Use,labeller = labeller(Use = uses.labs),drop=F)+
    coord_equal()
  png(file = paste0(table_fold_path,"/proba_G_space_across_months_uses.png"),
      width=2100, height=1200)
  plot(plot_allMonths_proba_G)
  dev.off()
  # For all uses, across months, prediction of presence/absence, in G space
  plot_allMonths_pred_G <- na.omit(df_proba_PCA_xy) %>%
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = as.factor(pres_pred))) +
    scale_fill_manual(values = c("1" = "#FDE725FF",
                                 "0" = "#440154FF"),
                      labels=c('Absence', 'Presence')) +
    theme(panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"),
          legend.position = "right",
          text = element_text(size=15),
          axis.text.x = element_text(angle=45)) +
    labs(fill="Predicted of\noutcome",
         x="Longitude",y="Latitude")+
    facet_grid(Month ~ Use,labeller = labeller(Use = uses.labs),drop=F)+
    coord_equal()
  png(file = paste0(table_fold_path,"/pred_G_space_across_months_uses.png"),
      width=2100, height=1200)
  plot(plot_allMonths_pred_G)
  dev.off()
  
  str(as.data.frame(as.data.table(df_proba_PCA_xy)))
  
  # Faire tourner les fonctions qui font un summary des proba utilisable pour hex
  df_hex_uses <- lapply(as.character(unique(df_proba_PCA_xy$Month)), 
                        function(x) f_hex_month(df=df_proba_PCA_xy, mois_run=x))
  df_hex_uses <- do.call(rbind, df_hex_uses)
  # Ordonner les facteurs dans le bon ordre pour la figure
  df_hex_uses$Use2 <- factor(df_hex_uses$Use2,      
                             levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba_pred"))
  df_hex_uses$Month <- factor(df_hex_uses$Month,      
                              levels = c("May","June","July","August","September"))
  
  ##  plot the results (all uses across months)
  plot_allMonths_pred_E <- ggplot(df_hex_uses) +
    geom_hex(aes(x = PCA1, y = PCA2, fill = pred),
             stat = "identity", alpha = 0.8,colour='grey') +
    scale_fill_viridis(limits=c(0,1)) +
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed") +
    labs(title = "Probability of Occurrence projected in Ecological Space",
         fill = "Predicted\noutcome")+
    facet_grid(Month ~ Use2, scales = "free",
               labeller = labeller(Use2 = uses.labs2)
    ) +
    theme(text = element_text(size=18),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"))
  # Save figure
  png(file = paste0(table_fold_path,"/pred_E_space_across_months_uses.png"),
      width=2100, height=1200)
  plot(plot_allMonths_pred_E)
  dev.off()
  
  # Plot density of geographics summarized in E space
  P <- ggplot(df_hex_uses) +
    geom_hex(aes(x = PCA1, y = PCA2, fill = n),
             stat = "identity", 
             colour='grey',
             alpha = 0.8) +
    scico::scale_fill_scico(palette = "vik")+
    #scale_fill_viridis(limits=c(0,1)) +
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed") +
    labs(fill = "Density of\ngeographic pixels")+
    facet_grid(Month ~ . ) +
    theme(text = element_text(size=18),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"))+
    coord_equal()
  # Save
  png(file = paste0(table_fold_path,"/density_geopix_E_space_months.png"),
      width=2100, height=1200)
  plot(P)
  dev.off()
  
  # Save csv to compute overlap
  write.csv(df_hex_uses,paste0(table_fold_path,"/df_pred_ACP_months.csv"))
}

# crée la figure :
# - G space : prédiction binaires, 6 uses X 1 summer
# - E space : prédiction binaires, 6 uses X 1 summer
Create.Plots.Summer <- function(){
  # Mettre en forme et aggréger sur les 5 mois
  df_pred_xy_summer <- df_proba_PCA_xy %>%
    filter(!is.na(PCA1)) %>%
    dplyr::select(-c(PCA1, PCA2, Proba, threshold_kept)) %>% 
    group_by(x,y,Use) %>%
    pivot_wider(names_from = Month, values_from = pres_pred) %>%
    mutate(sum_pred_accross_time = sum(August, July, June, May, September,na.rm=T))
  # Si au cours ne serait-ce que d'1 seul mois un usage est présent,
  # alors il est considéré présent pour l'été
  df_pred_xy_summer$pred_accross_time <- ifelse(df_pred_xy_summer$sum_pred_accross_time >0,1,0)
  # à conserver pour pouvoir ensuite projeter dans ACP summer
  write.csv(df_pred_xy_summer,paste0(table_fold_path,"/df_pred_xy_sumtimes.csv"))
  # Plot predicted outcome, for the summer, in G space
  plot_accrosstime_pred_G <- df_pred_xy_summer %>%
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = as.factor(pred_accross_time))) +
    scale_fill_manual(values = c("1" = "#FDE725FF",
                                 "0" = "#440154FF"),
                      labels=c('Absence', 'Presence')) +
    theme(panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"),
          legend.position = "right",
          text = element_text(size=15),
          axis.text.x = element_text(angle=45)) +
    labs(fill="Predicted\noutcome",
         #title = "Probability of Occurrence projected in Geographic Space",
         x="Longitude",y="Latitude")+
    facet_grid( ~ Use,labeller = labeller(Use = uses.labs),drop=F)+
    coord_equal()
  # Save
  png(file = paste0(table_fold_path,"/pred_G_space_sumtime_uses.png"),
      width=2100, height=1200)
  plot(plot_accrosstime_pred_G)
  dev.off()
  
  
  df_PCA_xy_summer <- na.omit(GetEspaceTableGlob(mois = NULL))
  # reduce table with predicted presence in the summer
  df_pred_xy_summer2 <- df_pred_xy_summer %>%
    dplyr::select(-c(August,June,July,May,September, sum_pred_accross_time)) %>%
    group_by(x,y) %>%
    pivot_wider(names_from = Use, values_from = pred_accross_time) %>%
    as.data.frame()
  
  # croiser df_pred_xy_summer et df_PCA_xy_summer : passer par le stack de raster
  # car merge() marche pas
  PCA1_r <- rasterFromXYZ(data.frame(x=df_PCA_xy_summer$x,y=df_PCA_xy_summer$y,z=df_PCA_xy_summer$PCA1), crs=EPSG_2154)
  PCA2_r <- rasterFromXYZ(data.frame(x=df_PCA_xy_summer$x,y=df_PCA_xy_summer$y,z=df_PCA_xy_summer$PCA2), crs=EPSG_2154)
  Pa_r <- rasterFromXYZ(data.frame(x=df_pred_xy_summer2$x,y=df_pred_xy_summer2$y,z=df_pred_xy_summer2$Pa_proba), crs=EPSG_2154)
  Rp_r <- rasterFromXYZ(data.frame(x=df_pred_xy_summer2$x,y=df_pred_xy_summer2$y,z=df_pred_xy_summer2$Rp_proba), crs=EPSG_2154)
  Vt_r <- rasterFromXYZ(data.frame(x=df_pred_xy_summer2$x,y=df_pred_xy_summer2$y,z=df_pred_xy_summer2$Vt_proba), crs=EPSG_2154)
  Ni_r <- rasterFromXYZ(data.frame(x=df_pred_xy_summer2$x,y=df_pred_xy_summer2$y,z=df_pred_xy_summer2$Ni_proba), crs=EPSG_2154)
  Lk_r <- rasterFromXYZ(data.frame(x=df_pred_xy_summer2$x,y=df_pred_xy_summer2$y,z=df_pred_xy_summer2$Lk_proba), crs=EPSG_2154)
  Co_r <- rasterFromXYZ(data.frame(x=df_pred_xy_summer2$x,y=df_pred_xy_summer2$y,z=df_pred_xy_summer2$Co_proba), crs=EPSG_2154)
  # Stack 6 uses + 2 PCAs
  PCA_pred_r <- stack(PCA1_r, PCA2_r,
                      Ni_r,Pa_r,Rp_r,Lk_r,Co_r,Vt_r)
  names(PCA_pred_r) <- c("PCA1","PCA2",paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba_pred"))
  plot(PCA_pred_r, colNA='black')
  # From rasters to table
  df_PCA_pred <- data.frame(data.table(PCA_pred_r[]))
  # df_PCA_pred_xy <- cbind(coordinates(PCA_pred_r),df_PCA_pred)
  # Aggregate geographic pixels values in E space
  xbnds <- range(c(df_PCA_pred$PCA1), na.rm = T)
  ybnds <- range(c(df_PCA_pred$PCA2), na.rm = T)
  #nbins <- 60
  nbins <- 100
  # rename PCA axis
  df_for_hex <- df_PCA_pred %>%
    #dplyr::select(-c(x,y)) %>%
    filter(!is.na(PCA1)) %>%
    rename(x=PCA1, y =PCA2)
  #df_for_hex <- filter(df_for_hex,rowSums(is.na(df_for_hex)) != ncol(df_for_hex))
  # Faire tourner les fonctions qui font un summary des proba utilisable pour hex
  hex_uses <- lapply(3:dim(df_for_hex)[2], function(x) f_hex(df2=df_for_hex, u=x))
  # Merge for all uses
  hex_uses <- Reduce(function(x,y) merge(x = x, y = y, by = c('PCA1','PCA2','n','cid')), 
                     hex_uses)
  # mettre en forme pour ggplot
  df_hex_uses <- hex_uses %>%
    pivot_longer(cols=ends_with("_proba_pred"),
                 names_to = "Use" ,values_to = "median_proba") %>%
    pivot_longer(cols=ends_with("_proba_pred_pred"),
                 names_to = "Use2" ,values_to = "pred")
  df_hex_uses$Use2 <- factor(df_hex_uses$Use2,      
                             levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba_pred_pred"))
  # nouveaux labels
  uses.labs3 <- uses.labs
  names(uses.labs3) <- paste0(names(uses.labs3),"_pred_pred")
  # Figure
  plot_accrosstime_pred_E <- ggplot(df_hex_uses) +
    geom_hex(aes(x = PCA1, y = PCA2, fill = pred),
             stat = "identity", 
             colour='grey',
             alpha = 0.8) +
    scale_fill_viridis(limits=c(0,1)) +
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed") +
    labs(title = "Probability of Occurrence projected in Ecological Space",
         fill = "Predicted\noutcome")+
    facet_grid( ~ Use2,# scales = "free",
                labeller = labeller(Use2 = uses.labs3)
    ) +
    theme(text = element_text(size=18),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"))+
    coord_equal()
  # Save
  png(file = paste0(table_fold_path,"/pred_E_space_sumtime_uses.png"),
      width=2100, height=1200)
  plot(plot_accrosstime_pred_E)
  dev.off()
  
  # Plot density of geographic pixels summarized in E space
  P <- ggplot(df_hex_uses) +
    geom_hex(aes(x = PCA1, y = PCA2, fill = n),
             stat = "identity", 
             colour='grey',
             alpha = 0.8) +
    scico::scale_fill_scico(palette = "vik")+
    #scale_fill_viridis(limits=c(0,1)) +
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_vline(xintercept = 0,linetype="dashed") +
    labs(fill = "Density of\ngeographic pixels")+
    theme(text = element_text(size=18),
          panel.background = element_rect(fill="white"),
          plot.background = element_rect(fill="white"),
          panel.grid.major = element_line(colour="grey"))+
    coord_equal()
  # Save
  png(file = paste0(table_fold_path,"/density_geopix_E_space_sumtime.png"),
      width=2100, height=1200)
  plot(P)
  dev.off()
  # Save csv to compute overlap
  write.csv(df_hex_uses,paste0(table_fold_path,"/df_pred_ACP_sumtimes.csv"))
}

# Fonction qui retourne le % d'aire partagée, pour u1, avec u2
ComputeOverlap <- function(table_pres_abs, liste_u1_u2){
  # # TEST
  # table_pres_abs <- df_4O_run
  # liste_u1_u2 <- c("Lk","Rp")
  
  nom_u1 <- liste_u1_u2[1]
  nom_u2 <- liste_u1_u2[2]
  
  A_u1 <- sum(table_pres_abs[,grep(nom_u1, names(table_pres_abs))], na.omit=T)
  A_u2 <- sum(table_pres_abs[,grep(nom_u2, names(table_pres_abs))], na.omit=T)
  
  table_pres_abs$inter <- table_pres_abs[, grep(nom_u1, names(table_pres_abs))] + 
    table_pres_abs[, grep(nom_u2, names(table_pres_abs))]
  
  A_u1_inter_u2 <- as.numeric(table(table_pres_abs$inter)["2"])
  
  overlap_u1 <- A_u1_inter_u2/A_u1
  overlap_u2 <- A_u1_inter_u2/A_u2
  
  df_overlap <- data.frame(u1=nom_u1,
                           u2=nom_u2,
                           A_u1 = A_u1,
                           A_u2 = A_u2,
                           A_inter = A_u1_inter_u2,
                           overlap_u1=overlap_u1,
                           overlap_u2=overlap_u2)
  
  return(df_overlap)
}

# Fonction qui calcule overlap pour un mois donné
Run_ComputeOverlap_monthly <- function(df , mois_run){
  # # TEST
  # mois_run <- 'May'
  # df <- df_pred_PCA_4O
  
  # Subset selon le mois
  df_4O_run <- df %>%
    filter(Month == mois_run) 
  #dplyr::select(-Month)
  
  df_overlap_run <- do.call(rbind, lapply(liste_paires_usages, function(x) 
    ComputeOverlap(table_pres_abs = df_4O_run, liste_u1_u2 = x)))
  df_overlap_run$Month <- mois_run
  str(df_overlap_run)
  return(df_overlap_run)
}

# met en forme une table d'overlap
makePaire <- function(df, mois_run,space){
  # #TEST
  # df <- df_overlap_G_monthly
  # mois_run <- 'June'
  # space <- "G"
  
  df_month <- df %>%
    filter(Month == mois_run) %>%
    dplyr::select(-c(A_u1,A_u2,A_inter)) %>%
    mutate(paire1 = paste(u1, u2,sep="_"),
           paire2 = paste(u2, u1,sep="_"))
  
  df_overlap <- data.frame(paire = c(df_month$paire1, df_month$paire2),
                           overlap = c(df_month$overlap_u1, df_month$overlap_u2),
                           Month= mois_run)
  names(df_overlap)[2] <- paste0("overlap_",space)
  return(df_overlap)
}
### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Arguments pour le choix du modèle à utiliser
type_donnees_glm <- "brute"
algorithme <- "glm"
fit <- "all_simple"
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")
gitCaractMilieu <- "C:/Users/perle.charlot/Documents/PhD/DATA/R_git/CaractMilieu"
input_path <- paste0(wd,"/input/")
table_fold_path <- paste0(output_path,"/niches/",type_donnees_glm,"/niche_overlap/",
                          fit,"/",algorithme,"/tables_pres/")
path_summary_models <- paste0(output_path,"/niches/",type_donnees_glm,
                              "/summary_model_kept.csv")
#### Données spatiales ####
dos_var_sp <- "C:/Users/perle.charlot/Documents/PhD/DATA/Variables_spatiales_Belledonne/"
limiteN2000 <- paste0(dos_var_sp, "/limites_etude/cembraie_N2000_limites.gpkg")
# Projection Lambert 93 (EPSG : 2154)
EPSG_2154 =  "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs "
#### Autre ####
liste.mois <- c("mai","juin","juillet","aout","septembre")
df.mois <- data.frame(nom_mois = liste.mois, numero_mois = c("05","06","07","08","09"))
# Translate french <-> english (for figures)
df_time <- data.frame(mois = c("mai","juin","juillet","aout","septembre"),
                      english_month = c("May","June","July","August","September"))

# Liste dimensions
liste.dim =  c("CA","B","PV","CS","D","I")
# Liste usages
liste.usages = c("Ni","Lk","Co","Pa","Rp","Vt")

# Tables
path_table_variables <- paste0(gitCaractMilieu,"/input/liste_variables.csv")
path_table_variables_dummies <- paste0(gitCaractMilieu,"/input/liste_variables_dummies.csv")
# Table correspondance entre dimension et couleur à utiliser dans les graphs
corresp_col <- data.frame(dim_name = c(liste.dim,"toutes"),
                         colour_dim = c("dodgerblue","darkgoldenrod1","darkgreen",
                                        "brown","blueviolet","darkgray","antiquewhite"))
# Labels names for plots
uses.labs <- c("Nesting","Sheep Grazing","Hiking",
               "Lek","Sheep Night Camping" ,"Mountain Bike")
names(uses.labs) <- paste0(c("Ni","Pa","Rp","Lk","Co","Vt"))

labs.env <- c("Biomass","Abiotic Conditions","Spatial Context",
              "Dynamic","Infrastructure" ,"Vegetation Physionomy")
names(labs.env) <- c("B","CA","CS",
                     "D","I", "PV")
### Programme -------------------------------------

# PARTIE A FAIRE TOUJOURS TOURNER ####
table_variables <- as.data.frame(fread(path_table_variables))
table_variable_dummies <- as.data.frame(fread(path_table_variables_dummies, dec=","))
summary_models <- as.data.frame(fread(path_summary_models,drop="V1"))
table_variables$Nom_var_initial <- table_variables$Nom
col_dim <- merge(rbind(table_variables, table_variable_dummies), 
                 corresp_col,by.x="Dimension", by.y="dim_name")

# Construction ou lecture des tables de présence #

# 1 - E space (indpdt proba ou obs)
t = try(as.data.frame(fread(paste0(table_fold_path,"/df_PCA_glob_xy_all_months.csv"),drop="V1")),
        silent = TRUE)
if(inherits(t, "try-error")){
  df_PCA_glob <- as.data.frame(do.call(rbind,lapply(liste.mois, function(x) GetEspaceTableGlob(mois=x))))
  write.csv(df_PCA_glob,paste0(table_fold_path,"/df_PCA_glob_xy_all_months.csv"))
} else{
  df_PCA_glob <- as.data.frame(fread(paste0(table_fold_path,"/df_PCA_glob_xy_all_months.csv"),drop="V1"))
  }

# 2 - G space
#   a - proba
t = try(as.data.frame(fread(paste0(table_fold_path,"/df_proba_uses_xy_all_months.csv"),drop="V1")),
        silent = TRUE)
if(inherits(t, "try-error")){
  df_uses_proba <- as.data.frame(do.call(rbind,lapply(liste.mois, 
                                                      function(x) GetGspaceTable(mois=x,
                                                                                 data_type="proba"))))
  write.csv(df_uses_proba ,paste0(table_fold_path,"/df_proba_uses_xy_all_months.csv"))
}else{
  df_uses_proba <- as.data.frame(fread(paste0(table_fold_path,"/df_proba_uses_xy_all_months.csv"),drop="V1"))
}
#   b - obs
t = try(as.data.frame(fread(paste0(table_fold_path,"/df_obs_uses_xy_all_months.csv"),drop="V1")),
        silent = TRUE)
if(inherits(t, "try-error")){
  df_uses_obs <- as.data.frame(do.call(rbind,lapply(liste.mois, 
                                                      function(x) GetGspaceTable(mois=x,
                                                                                 data_type="obs"))))
  write.csv(df_uses_obs,paste0(table_fold_path,"/df_obs_uses_xy_all_months.csv"))
}else{
  df_uses_obs <- as.data.frame(fread(paste0(table_fold_path,"/df_obs_uses_xy_all_months.csv"),drop="V1"))
}
  
# 3 - merge E space - G space proba - G space obs ?
t = try(as.data.frame(fread(paste0(table_fold_path,"/df_PCA_glob_proba_obs_xy_all_months.csv"),drop="V1")),
        silent = TRUE)
if(inherits(t, "try-error")){
  df_proba_PCA_xy <- merge(df_uses_proba,df_PCA_glob, by=c("x","y","Month"),all=T)
  df_proba_PCA_xy <- na.omit(df_proba_PCA_xy)
  df_obs_PCA_xy <- merge(df_uses_obs,df_PCA_glob, by=c("x","y","Month"),all=T)
  df_obs_PCA_xy <- na.omit(df_obs_PCA_xy)  
  df_proba_obs_PCA_xy <- merge(df_proba_PCA_xy, df_obs_PCA_xy, by=c("x","y","Month","Use","PCA1","PCA2"))
  write.csv(df_proba_obs_PCA_xy,paste0(table_fold_path,"/df_PCA_glob_proba_obs_xy_all_months.csv"))
}else{
  df_proba_obs_PCA_xy <- as.data.frame(fread(paste0(table_fold_path,"/df_PCA_glob_proba_obs_xy_all_months.csv"),drop="V1"))
}

df_proba_obs_PCA_xy$Month <- factor(df_proba_obs_PCA_xy$Month,      
                                    levels = c("May","June","July",
                                               "August","September"))
# TODO
# 4 - display in E and G spaces

# TODO
# 5 - compute overlaps
#   a - with proba
#   b - with obs

# TODO : 
# - dans Create.Plots.Monthly, ajouter un argument en proba ou obs pour calculer chvch

### Plots des proba dans espaces G et E ####

## find the bounds for the complete data 
xbnds <- range(c(df_proba_PCA_xy$PCA1), na.omit=T)
ybnds <- range(c(df_proba_PCA_xy$PCA2), na.omit=T)
#nbins <- 60
nbins <- 100

Create.Plots.Monthly()

Create.Plots.Summer()

### Calcul des chevauchements dans espaces G et E ####

# liste des paires d'usages
df_combi <- combn(liste.usages, 2)
liste_paires_usages <- lapply(1:dim(df_combi)[2], function(x) paste(df_combi[,x]))

## Par mois

# G space (4O = for overlap analysis)
df_pred_xy_4O <- df_proba_PCA_xy %>%
  dplyr::select(-c('threshold_kept','Proba','PCA1','PCA2')) %>%
  pivot_wider(names_from = 'Use', values_from = 'pres_pred')
# répéter pour chaque mois
df_overlap_G_monthly <- do.call(rbind, 
                                lapply(df_time$english_month, 
                                       function(x) Run_ComputeOverlap_monthly(df = df_pred_xy_4O, 
                                                                              mois_run = x)))
# E space
df_hex_uses <- fread(paste0(table_fold_path,"/df_pred_ACP_months.csv"),drop="V1")
df_pred_PCA_4O <- df_hex_uses %>%
  group_by(Month, PCA1, PCA2) %>%
  dplyr::select(-c('cid','Use','median_proba')) %>%
  distinct() %>%
  pivot_wider(names_from = 'Use2', values_from = 'pred') %>%
  as.data.frame
# répéter pour chaque mois
df_overlap_E_monthly <- do.call(rbind, lapply(df_time$english_month, 
                      function(y) Run_ComputeOverlap_monthly(df = df_pred_PCA_4O, 
                                                             mois_run=y)))

# Plot similarité E vs intensité chevauchement G
df_overlap_G_monthly_2 <- do.call(rbind,lapply(df_time$english_month, function(x) makePaire(df=df_overlap_G_monthly,
                                                                  mois_run = x,
                                                                  space="G")))
df_overlap_E_monthly_2 <- do.call(rbind,lapply(df_time$english_month, function(x) makePaire(df=df_overlap_E_monthly,
                                                                  mois_run = x,
                                                                  space="E")))
df_overlap_monthly <- merge(df_overlap_G_monthly_2, df_overlap_E_monthly_2,
      by=c('paire', 'Month'))
# Ordonner les mois
df_overlap_monthly$Month <- factor(df_overlap_monthly$Month,      
                            levels = c("May","June","July","August","September"))

# Classifier les paires selon appartenance activités humaines ou faune sauvage
test <- as.data.frame(t(combn(liste.usages, 2)))
test2 <- data.frame(use = liste.usages,
                    category = c("FS","FS","PA","PA","PR","PR"), #FS = faune sauvage, PA = pastoralisme, PR = pratique récréative
                    category_big = c("FS","FS","AH","AH","AH","AH")) # AH = activité humaine

a <- merge(test,test2, by.x="V1",by.y="use",sort=F)
a <- a %>%
  dplyr::select(-V2)
b <- merge(test,test2, by.x="V2",by.y="use",sort=F)
b <- b %>%
  dplyr::select(-V1)

d1 <- data.frame(paire = paste(a$V1,b$V2,sep="_"),
           category = paste(a$category,b$category,sep="_"),
           category_big = paste(a$category_big,b$category_big,sep="_"))

a <- merge(test,test2, by.x="V2",by.y="use",sort=F)
a <- a %>%
  dplyr::select(-V1)
b <- merge(test,test2, by.x="V1",by.y="use",sort=F)
b <- b %>%
  dplyr::select(-V2)

d2 <- data.frame(paire = paste(a$V2,b$V1,sep="_"),
                 category = paste(a$category,b$category,sep="_"),
                 category_big = paste(a$category_big,b$category_big,sep="_"))

corresp_classif_paire <- rbind(d1,d2)





# au global (tous les mois mélanges, toutes les paires)
plot_chvch_glob <- df_overlap_monthly %>%
  ggplot(aes(x=overlap_G, y = overlap_E)) +
  geom_point() +
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Intensity of Geographic Overlap",
       y = "Ecological Similarity") +
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))
plot_chvch_glob

#save plot
png(file = paste0(table_fold_path,"/global_overlap.png"),width=1400, height=800)
plot(plot_chvch_glob)
dev.off()


Pmonth <- df_overlap_monthly %>%
  ggplot(aes(x=overlap_G, y = overlap_E, col=Month)) +
  geom_point() +
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Intensity of Geographic Overlap",
       y = "Ecological Similarity") +
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))
png(file = paste0(table_fold_path,"/global_overlap_colorMonth.png"),width=1400, height=800)
plot(Pmonth)
dev.off()

# TODO : comparer si on trouve une diff avec le graph summer
# normalement il devrait y avoir moins de points sur l'autre,
# mais est-ce que même tendance ?

df_overlap_monthly2 <- merge(df_overlap_monthly, corresp_classif_paire, by='paire')

df_overlap_monthly2 <- df_overlap_monthly2 %>%
  distinct()

df_overlap_monthly2 %>%
  ggplot(aes(x=overlap_G, y = overlap_E, col=as.factor(category_big))) +
  geom_point() +
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Intensity of Geographic Overlap",
       y = "Ecological Similarity") +
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))

df_overlap_monthly2 %>%
  ggplot(aes(x=overlap_G, y = overlap_E, col=as.factor(category))) +
  geom_point() +
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Intensity of Geographic Overlap",
       y = "Ecological Similarity") +
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))

df_overlap_monthly2 %>%
  ggplot(aes(x=overlap_G, y = overlap_E, col=as.factor(category_big))) +
  geom_point() +
  ylim(0,1)+
  facet_grid( ~ Month)+
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Intensity of Geographic Overlap",
       y = "Ecological Similarity") +
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))

# par mois
df_overlap_monthly %>%
  ggplot(aes(x=overlap_G, y = overlap_E)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid( ~ Month)+
  labs(x = "Intensity of Geographic Overlap",
       y = "Ecological Similarity") +
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))

# par paire
plot_chvch_bypaire <- df_overlap_monthly %>%
  ggplot(aes(x=overlap_G, y = overlap_E, col=Month)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Intensity of Geographic Overlap",
       y = "Ecological Similarity") +
  facet_wrap(~paire) +
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))
  
png(file = paste0(table_fold_path,"/bypair_overlap.png"),width=1400, height=800)
plot(plot_chvch_bypaire)
dev.off()



df_overlap_monthly
df_overlap_monthly2

write.csv(df_overlap_monthly,paste0(table_fold_path,"/df_overlap_monthly.csv"))

  write.csv(df_overlap_monthly2,paste0(table_fold_path,"/df_overlap_monthly2.csv"))


# TODO

## Pour tout l'été

# G space
head(df_proba_PCA_xy)
df_proba_PCA_xy %>%
  filter
# E space
df_hex_uses <- fread(paste0(table_fold_path,"/df_pred_ACP_sumtimes.csv"),drop="V1")









# Analyse points extrêmes ######
head(df_Vt_Co)

df <- df_proba_PCA_xy %>% 
  dplyr::select('x','y','Month','PCA1','PCA2') %>%
  distinct()
df_Vt_Co2 <- merge(df_Vt_Co, df, by=c("x","y","Month"))

xbnds <- range(c(df_Vt_Co2$PCA1), na.omit=T)
ybnds <- range(c(df_Vt_Co2$PCA2), na.omit=T)
#nbins <- 60
nbins <- 100


df_Vt_Co2$Month <- factor(df_Vt_Co2$Month,
                          levels = c("May","June","July","August","September"))




test <- na.omit(df_Vt_Co2) %>%
  ggplot() +
  stat_summary_hex(aes(x=PCA1, y=PCA2, z= Obs),
                   #fun = function(x) sum(x), 
                   bins=60,
                   colour='grey'
  )+
  scale_fill_viridis(na.value = "transparent",
                     limits=c(0,1)) +
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed") +
  facet_grid(Use~ Month,labeller = labeller(Use = uses.labs.obs),drop=F)+
  labs(title = "Mean Observed Occurrence projected in Ecological Space"
       #fill = "Median probability\nof occurrence"
  )+
  coord_equal()+
  theme(text = element_text(size=18),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour="grey"))

png(file = paste0(table_fold_path,"/Vt_Co_Espace2.png"),
    width=2100, height=1200)
plot(test)
dev.off()

######