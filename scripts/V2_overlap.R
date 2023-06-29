### Titre -------------------------------------
# Nom : Calcul indices overlap
# Auteure : Perle Charlot
# Date de création : 26-06-2023
# Dates de modification : -2023

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
### Fonctions -------------------------------------

# Fonction qui sort une table avec coordonnées xy et proba de tous les usages, pour un mois donné
GetGspaceTable <- function(mois,
                           usages_etudies = NULL, #si on veut une sélection d'usages particulière, 
                           # de base, la fonction fait pour tous les usages rencontrés le mois étudié
                           type_donnees = "brute",
                           fit = "all_simple",
                           algorithme = "glm"){
  # # TEST
  # mois = "juin"
  # usages_etudies = NULL
  
  # Create directories
  table_fold_path_month <- paste0(table_fold_path,mois,"/") 
  if(!dir.exists(table_fold_path_month)){dir.create(table_fold_path_month, recursive = T)}
  
  # Translate french <-> english (for figures)
  english_month <- df_time$english_month[df_time$mois == mois]
  
  # load raster of proba for month studied
  list_rast <- base::list.files(path = paste0(output_path,"/niches/", type_donnees,"/"),
                                pattern = ".tif$", 
                                recursive = T, full.names = T)
  # trier le mois
  list_rast <- list_rast[grep(mois,list_rast)]
  
  # Trier le threshold
  list_rast <- list_rast[grep("TSS",list_rast)]
  
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
  stack_uses_proba <- stack_uses %>% subset(grep("proba", names(stack_uses)))
  plot(stack_uses_proba)
  
  # raster -> df
  df_uses_proba <- data.frame(data.table(stack_uses_proba[]))
  df_uses_proba <- cbind(coordinates(stack_uses_proba),df_uses_proba)
  
  # pivot
  df_uses_proba <- df_uses_proba %>% pivot_longer(cols=ends_with("proba"),
                                                  names_to = "Use", 
                                                  values_to = "Proba")
  # Reordering group factor levels
  df_uses_proba$Use <- factor(df_uses_proba$Use,      
                              levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba"))
  
  
  df_uses_proba <- as.data.frame(df_uses_proba) 
  
  # Plot map (G space) PROBA PRESENCE
  map <- na.omit(df_uses_proba) %>%
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = Proba)) +
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
  map
  # si nécessaire, ne garder que certains usages sur les cartes
  if(!is.null(usages_etudies)){
    map <- na.omit(df_uses_proba) %>%
      ggplot() +
      geom_raster(aes(x = x, y = y, fill = Proba)) +
      scale_fill_viridis(limits=c(0,1)) +
      theme(panel.background = element_rect(fill="white"),
            plot.background = element_rect(fill="white"),
            panel.grid.major = element_line(colour="grey"),
            legend.position = "right",
            text = element_text(size=15),
            axis.text.x = element_text(angle=45)) +
      labs(fill="Probability of\noccurrence", title = english_month,
           x="Longitude",y="Latitude")+
      facet_wrap(.~Use,labeller = labeller(Use = uses.labs),drop=T)+
      coord_equal()
  }
  #save plot
  png(file = paste0(table_fold_path_month,"/map_proba_uses.png"),width=1400, height=800)
  plot(map)
  dev.off()

  #select threshold for each use
  df_uses_proba2 <- summary_models %>%
    select(use,threshold_kept) %>%
    rename(Use = use) %>%
    mutate(Use = paste0(Use,"_proba")) %>% 
    merge(df_uses_proba, by="Use",all=T)
  df_uses_proba2$pres_pred <- ifelse(df_uses_proba2$Proba > df_uses_proba2$threshold,1,0)
  
  # Reordering group factor levels
  df_uses_proba2$Use <- factor(df_uses_proba2$Use,      
                              levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba"))
  
  # Plot map (G space) PRESENCE BINERISE (pred)
  map2 <- na.omit(df_uses_proba2) %>%
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
  # save df
  df_uses_proba2$Month <- english_month
  # make different name if not on all uses
  df_name <- paste0("/df_coords_proba_uses",paste0(usages_etudies,collapse="_"),".csv")
  write.csv(df_uses_proba2, paste0(table_fold_path_month,df_name))
  
  return(df_uses_proba2)
}



# Fonction qui sort une table avec coordonnées xy et valeurs PCA1 PCA2 pour ACP globale, pour un mois donné
GetEspaceTableGlob <- function(mois,
                               usages_etudies = NULL, #si on veut une sélection d'usages particulière, 
                               # de base, la fonction fait pour tous les usages rencontrés le mois étudié
                               type_donnees = "brute",
                               fit = "all_simple",
                               algorithme = "glm"){
  # # TEST
  # mois = "juin"
  # usages_etudies = NULL
  
  # Create directories
  table_fold_path_month <- paste0(table_fold_path,mois,"/") 
  if(!dir.exists(table_fold_path_month)){dir.create(table_fold_path_month, recursive = T)}
  
  # Translate french <-> english (for figures)
  english_month <- df_time$english_month[df_time$mois == mois]
  # Get PCA axes from global analysis
  PCA1 <- stack(list.files(path=paste0(gitCaractMilieu,"/output/ACP/ACP_avec_ponderation/summer/pred_month/",mois),
                           "axe1", full.names = T))
  
  PCA2 <- stack(list.files(path=paste0(gitCaractMilieu,"/output/ACP/ACP_avec_ponderation/summer/pred_month/",mois),
                           "axe2", full.names = T))
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
  df_env$Month <- english_month
  # make different name if not on all uses
  df_name <- "/df_coords_PCA_axes.csv"
  write.csv(df_env, paste0(table_fold_path_month,df_name))
  
  # df_env2$Month <- english_month
  # df_env2$dimension <- "global"
  
  return(df_env)
}

### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")
gitCaractMilieu <- "C:/Users/perle.charlot/Documents/PhD/DATA/R_git/CaractMilieu"
input_path <- paste0(wd,"/input/")

type_donnees <- "brute"
algorithme <- "glm"
fit <- "all_simple"

table_fold_path <- paste0(output_path,"/niches/",type_donnees,"/niche_overlap/",
                          fit,"/",algorithme,"/tables_pres/")

# # directories paths
# Gspace_path <- paste0(output_path,"/niches/",type_donnees,"/niche_overlap/",
#                       fit,"/",algorithme,"/G_space/") 
# Espace_path <- paste0(output_path,"/niches/",type_donnees,"/niche_overlap/",
#                       fit,"/",algorithme,"/E_space/") 

#### Données spatiales ####
dos_var_sp <- "C:/Users/perle.charlot/Documents/PhD/DATA/Variables_spatiales_Belledonne/"
limiteN2000 <- paste0(dos_var_sp, "/limites_etude/cembraie_N2000_limites.gpkg")

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



# English names
labs.vars <- c("GDD","NDVI","Water Balance",
               "No Foliage","Rare Foliage","Medium Foliage","Abundant Foliage",
               "Easting","Snow Height","Erosion","Number of Landforms",
               "Freezing Days","Low Nebulosity","High Nebulosity","No Thawing",
               "Northing","Slope","No Rain","Landforms Diversity",
               "Solar Radiation","Low Temperature","High Temperature","TWI",
               "Low Wind","High Wind","Land","Open Water",
               "Canyons","Shallow Valleys","Headwaters","U-Shape Valleys",
               "Plains","Open Slopes","Upper Slopes","Local Ridges","Midslope Ridges",
               "Mountain Tops",
               "Distance from Open Waters","Distance from Forest",
               "Distance from Infrastructure","Similar Habitat (1km)",
               "Similar Habitat (100m)","Similar Habitat (250m)","Similar Habitat (500m)",
               "Size of Habitat Patch","Access Time","Visibility",
               "Temperature Shift","No Avalanche","Avalanche Area",
               "Infrastructure Cover","Natural Environment","Modified Environment",
               "Road","Building",
               "Outside Protected Area","Protected Area","Restricted Area",
               "Canopy Maximum Height","Number of Vegetation Stratum",
               "High Vegetation Cover Penetrability","Medium Vegetation Cover Penetrability",
               "Low Vegetation Cover Penetrability")
names(labs.vars) <- c("GDD","NDVI","P_ETP",
                      "abondance_feuillage0","abondance_feuillage1",
                      "abondance_feuillage2","abondance_feuillage3",
                      "easting_25m","htNeigmean","LS_factor","nb_distinct_landform",
                      "nbJgel","nbJneb10","nbJneb90","nbJssdegel",
                      "northing_25m","pente_25m","rain0","shannon_landform",
                      "SWDmean","t10","t90","TWI_25m","wind10","wind90",
                      "presence_eau0","presence_eau1",
                      "landform_25m1","landform_25m2","landform_25m3",
                      "landform_25m4","landform_25m5","landform_25m6",
                      "landform_25m7","landform_25m8","landform_25m9",
                      "landform_25m10",
                      "distance_eau","distance_foret_IGN_cout_pente",
                      "distance_infrastructure","habitat_similaire_1000m",
                      "habitat_similaire_100m","habitat_similaire_250m","habitat_similaire_500m",
                      "taille_patch_habitat_m2","temps_acces","visibilite_mediane",
                      "diffT","presence_avalanche0","presence_avalanche1",
                      "pourcentage_infrastructures","degre_artif0","degre_artif1",
                      "degre_artif2","degre_artif3",
                      "degre_interdiction0","degre_interdiction1","degre_interdiction2",
                      "ht_physio_max","nb_strates","penetrabilite1",
                      "penetrabilite2","penetrabilite3")

uses.labs <- c("Nesting","Sheep Grazing","Hiking",
               "Lek","Sheep Night Camping" ,"Mountain Bike")
names(uses.labs) <- paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba")

labs.env <- c("Biomass","Abiotic Conditions","Spatial Context",
              "Dynamic","Infrastructure" ,"Vegetation Physionomy")
names(labs.env) <- c("B","CA","CS",
                     "D","I", "PV")

path_summary_models <- paste0(output_path,"/niches/",type_donnees,
                              "/summary_model_kept.csv")

### Programme -------------------------------------

table_variables <- as.data.frame(fread(path_table_variables))
table_variable_dummies <- as.data.frame(fread(path_table_variables_dummies, dec=","))
summary_models <- as.data.frame(fread(path_summary_models,drop="V1"))
table_variables$Nom_var_initial <- table_variables$Nom
col_dim <- merge(rbind(table_variables, table_variable_dummies), 
                 corresp_col,by.x="Dimension", by.y="dim_name")




# 1 : récupérer prédctions spatialisées, par mois

df_proba_uses <- as.data.frame(do.call(rbind,lapply(liste.mois, function(x) GetGspaceTable(mois=x))))
str(df_proba_uses)
write.csv(df_proba_uses,paste0(table_fold_path,"/df_proba_uses_xy_all_months.csv"))

# 2 : par paire, calculer indice chvch = INTENSITE CHVCH ESP GEO (peut se faire en restant en raster)
# 3 : transformer raster en table et projeter dans esp de niche

df_PCA_glob <- as.data.frame(do.call(rbind,lapply(liste.mois, function(x) GetEspaceTableGlob(mois=x))))
write.csv(df_PCA_glob,paste0(table_fold_path,"/df_PCA_glob_xy_all_months.csv"))

df_proba_PCA_xy <- merge(df_proba_uses,df_PCA_glob, by=c("x","y","Month"),all=T)
df_proba_PCA_xy$Month <- factor(df_proba_PCA_xy$Month,      
                                  levels = c("May","June","July",
                                             "August","September"))
write.csv(df_proba_PCA_xy,paste0(table_fold_path,"/df_PCA_glob_proba_xy_all_months.csv"))
# Reordering group factor levels
df_proba_PCA_xy$Use <- factor(df_proba_PCA_xy$Use,      
                             levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba"))


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
  theme(text = element_text(size=18))
png(file = paste0(table_fold_path,"/proba_E_space_across_months_uses.png"),
    width=2100, height=1200)
plot(plot_allMonths_proba_E)
dev.off()
# For all uses, across months, probabilities, in G space
plot_allMonths_proba_G <-na.omit(df_proba_PCA_xy) %>%
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

# For all uses, across months, prediction of presence/absence, in E space
# Comment faire ?
# Car un hexagone = plusieurs pixels donc plusieurs probabilités ou plusieurs pres/abs
# 2 options :
# - faire la médiane des probabilités puis binariser --> semble faire des trucs chelous
# - partir des pres/abs déjà seuillé, faire la médiane, binariser


# https://stackoverflow.com/questions/39296198/operation-between-stat-summary-hex-plots-made-in-ggplot2

test <- na.omit(df_proba_PCA_xy)

## find the bounds for the complete data 
xbnds <- range(c(test$PCA1))
ybnds <- range(c(test$PCA2))
nbins <- 60
#  function to make a data.frame for geom_hex that can be used with stat_identity
makeHexData <- function(df) {
  h <- hexbin(df$x, df$y, nbins, xbnds = xbnds, ybnds = ybnds, IDs = TRUE)
  data.frame(hcell2xy(h),
             z = tapply(df$z, h@cID, FUN = function(z) median(z)), # proba médiane
             n = tapply(df$z, h@cID, FUN = function(z) length(z)), # nb pixels géo
             # sd = tapply(df$z, h@cID, FUN = function(z) sd(z)),
             # mean = tapply(df$z, h@cID, FUN = function(z) mean(z)),
             # max = tapply(df$z, h@cID, FUN = function(z) max(z)),
             cid = h@cell)
}



f_hex_month <- function(mois_run, df){
  # TEST
  mois_run <- "July"
  df <- test
  
  # Pour un mois donné, tableau avec tous les usages
  test_all <- df %>%
    filter(Month == mois_run) %>%
    select(Use,PCA1,PCA2, Proba) %>%
    pivot_wider(names_from = Use, values_from = Proba) %>%
    rename(x=PCA1, y =PCA2)
  # Compute hexagonal summary of proba for each use
  f_hex <- function(df2, u){
    # Subset pour l'usage u
    A <- df2[,c(1:2,u)]
    name_u <- names(A)[3]
    names(A)[3] <- "z"
    # compute grille hexagonale
    Ahex <- makeHexData(A)
    
    # rename and add vars
    Ahex <- Ahex %>%
      rename(PCA1 =x, PCA2 = y) %>%
      mutate(Month = mois_run)
    # Filtrer par threshold
    th_u <- summary_models$threshold_kept[which(substr(name_u,1,2) %in% summary_models$use)] 
    Ahex$median_proba_pres <- ifelse(Ahex$z > th_u,
                                     1,0)
    names(Ahex)[3] <- name_u
    names(Ahex)[dim(Ahex)[2]] <- paste0(name_u,"_pred")
    return(Ahex)
  }
  hex_uses <- lapply(3:8, function(x) f_hex(df2=df, u=x))
  
  # TODO
  
  hex_uses <- Reduce(function(x,y) merge(x = x, y = y, by = c('PCA1','PCA2','Month','n','cid')), 
                     hex_uses)
  
  hex_uses <- hex_uses %>%
    pivot_longer(cols=ends_with("_proba"),
                 names_to = "Use" ,values_to = "median_proba") %>%
    pivot_longer(cols=ends_with("_proba_pred"),
                 names_to = "Use2" ,values_to = "pred")
  return(hex_uses)
}

f_hex_month("July", test)
  
df_hex_uses <- lapply(as.character(unique(test$Month)), 
                   function(x) f_hex_month(df=test, mois_run=x))

hex_uses <- Reduce(function(x,y) merge(x = x, y = y, by = c('PCA1','PCA2','Month','n','cid')), 
                   hex_uses)

##  plot the results
ggplot(testAB) +
  geom_hex(aes(x = PCA1, y = PCA2, fill = pred),
           stat = "identity", alpha = 0.8,colour='grey') +
  # scale_fill_manual(values = c("1" = "#FDE725FF",
  #                              "0" = "#440154FF"),
  #                   labels=c('Absence', 'Presence')) +
  scale_fill_viridis(limits=c(0,1)) +
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed") +
  labs(title = "Probability of Occurrence projected in Ecological Space",
       fill = "Predicted\noutcome")+
  
  facet_grid(Month ~ Use2, scales = "free"
             #labeller = labeller(Use2 = uses.labs)
             ) +
  
  theme(text = element_text(size=18))




testB$Co_proba_pred <- as.integer(testB$Co_proba_pred)
ggplot(testB) +
  geom_hex(aes(x = PCA1, y = PCA2, fill = Co_proba_pred ),
           stat = "identity", alpha = 0.8,colour='grey')
# THIS WORKS

testB %>%
  mutate(pred = factor(Co_proba_pred, levels =c(0,1))) %>%
  ggplot() +
  geom_hex(aes(x = PCA1, y = PCA2, fill = as.character(Co_proba_pred)),
           stat = "identity", alpha = 0.8,colour='grey')
# POURQUOI CA MARCHE PAS PUTAAAAIN







#
# Reordering group factor levels
df_uses_proba2$Use <- factor(df_uses_proba2$Use,      
                             levels = paste0(c("Ni","Pa","Rp","Lk","Co","Vt") ,"_proba"))


# 4 : par paire, calculer indice chvch = SIMILARITE DE NICHE


ComputeOverlap <- function(table_pres_abs, nom_u1, nom_u2){
  
  A_u1 <- sum(table_pres_abs["nom_u1"])
  A_u2 <- sum(table_pres_abs["nom_u2"])
  
  table_pres_abs$inter <- table_pres_abs["nom_u1"] + table_pres_abs["nom_u2"]
  
  table(table_pres_abs$inter)
  
}
