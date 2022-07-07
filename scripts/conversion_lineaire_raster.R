### Titre -------------------------------------
# Nom : Conversion usages ponctuels en raster binaire présence/absence
# Auteure : Perle Charlot
# Date de création : 24-06-2022
# Dates de modification : 03-07-2022

### Librairies -------------------------------------
library(raster)
library(fasterize)
library(sf)
library(rgdal)

# library(purrr)
# library(dplyr)
# library(data.table)
# library(stringi)
### Fonctions -------------------------------------

# Fonction qui rasterize un usage ponctuel
RasterizeUsage <- function(chemin_ponctuel_usage, 
                           raster_de_ref, 
                           nom_sortie,
                           liste_mois){
  
  # # TEST
  # chemin_ponctuel_usage = path_gpkg_randonnee
  # raster_de_ref = raster_ref
  # nom_sortie = "randonnee_pedestre" 
  # liste_mois = c("mai","juin","juillet","aout","septembre")
  
  
  # Import donné linéaire
  usage_ponctuel <- st_read(chemin_ponctuel_usage)
  #usage_ponctuel <- st_transform(usage_ponctuel,EPSG_2154)
  
  # # Application d'un buffer pour transfo linéaire -> polygone
  # # st_buffer n'a pas l'air de fonctionner avec des lignes
  # usage_ponctuel_buf <- raster::buffer(as_Spatial(usage_ponctuel), 10) 
  # usage_ponctuel <- st_cast(st_as_sf( usage_ponctuel_buf), "POLYGON")
  # FAIT SOUS QGIS car ça marche pas sur R ... buffer de 10 m
  
  usage_raster <- fasterize(usage_ponctuel, raster_de_ref, background=0)
  plot(usage_raster, colNA="black")
  
  # La carte des usages est à écrire dans les mois correspondant
  liste_chemins_outputs <- paste0(output_path,"/usages_recreatifs/par_periode/",liste_mois,"/", nom_sortie,".tif")
  
  for(i in liste_chemins_outputs){ # ok c'est moche de faire une boucle for, mais il n'y a pas bcp de rasters à écrire ...
    writeRaster(usage_raster, i, overwrite=TRUE)
  }
  cat(paste0("\nUsage ", nom_sortie, " a été rasterisé."))
  
}
### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")
# Projection Lambert 93 (EPSG : 2154)
EPSG_2154 =  "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs "

#### Données spatiales ####

# Dossier des variables spatiales & chemins des fichiers
dos_var_sp <- "C:/Users/perle.charlot/Documents/PhD/DATA/Variables_spatiales_Belledonne/"
# MNT 25m CALé sur le grille de REFERENCE
path_MNT <- paste0(dos_var_sp ,"/Milieux/IGN/mnt_25m_belledonne_cale.tif")

# Fichier spatial des tracés de randonnée pédestre
# durant tout l'été
path_gpkg_randonnee <- paste0(output_path, "/usages_recreatifs/randonnee/randonnee_buffer.gpkg")
# seulement durant l'ouverture des remontées mécaniques
path_gpkg_randonneeRMK <- paste0(output_path, "/usages_recreatifs/randonnee_RMK/randonnee_RMK_buffer.gpkg")

# Fichier spatial des tracés de VTT
path_gpkg_VTT <- paste0(output_path, "/usages_recreatifs/VTT/VTT_buffer.gpkg")
path_gpkg_VTTRMK <- paste0(output_path, "/usages_recreatifs/VTT_RMK/VTT_RMK_buffer.gpkg")

#### Tables ####

### Programme -------------------------------------

# Import raster de référence
raster_ref <- raster(path_MNT)

#### Randonnée pédestre ####

RasterizeUsage(path_gpkg_randonnee, raster_ref, "randonnee_pedestre", 
               c("mai","juin","juillet","aout","septembre"))
RasterizeUsage(path_gpkg_randonneeRMK, raster_ref, "randonnee_pedestreRMK", c("juillet","aout"))

# TODO : réfléchir à comment/quand bind les usages rando et randoRMK en juillet/août

path_rast_juillet <- list.files(paste0(output_path,"/usages_recreatifs/par_periode/juillet/"),
           ".tif", full.names = TRUE)
rast_juillet <- stack(path_rast_juillet)
plot(rast_juillet)

usage_rando_juillaout <- sum(rast_juillet)
usage_rando_juillaout[usage_rando_juillaout ==2 ] <- 1
#plot(usage_rando_juillaout)# ça rajoute vraiment pas grand chose, mais bon

# Supprimer les rasters calculés en disjoints
file.remove(path_rast_juillet)
file.remove(list.files(paste0(output_path,"/usages_recreatifs/par_periode/aout/"),
                                ".tif", full.names = TRUE))
# Enregistrer le nouveau raster combiné
writeRaster(usage_rando_juillaout, 
            paste0(output_path,"/usages_recreatifs/par_periode/juillet/randonnee_pedestre_juillet.tif"))
writeRaster(usage_rando_juillaout, 
            paste0(output_path,"/usages_recreatifs/par_periode/aout/randonnee_pedestre_aout.tif"))

#### VTT ####

RasterizeUsage(path_gpkg_VTT, raster_ref, "VTT", c("mai","juin","juillet","aout","septembre"))
RasterizeUsage(path_gpkg_VTTRMK, raster_ref, "VTT_RMK", c("juillet","aout"))


path_rast_juillet <- list.files(paste0(output_path,"/usages_recreatifs/par_periode/juillet/"),
                                ".tif", full.names = TRUE)
rast_juillet <- stack(path_rast_juillet)
plot(rast_juillet, colNA='black')

usage_rando_juillaout <- sum(rast_juillet)
usage_rando_juillaout[usage_rando_juillaout ==2 ] <- 1
plot(usage_rando_juillaout, colNA='black')# ça rajoute vraiment pas grand chose, mais bon

# Supprimer les rasters calculés en disjoints
file.remove(path_rast_juillet)
file.remove(list.files(paste0(output_path,"/usages_recreatifs/par_periode/aout/"),
                       ".tif", full.names = TRUE))
# Enregistrer le nouveau raster combiné
writeRaster(usage_rando_juillaout, 
            paste0(output_path,"/usages_recreatifs/par_periode/juillet/VTT_juillet.tif"))
writeRaster(usage_rando_juillaout, 
            paste0(output_path,"/usages_recreatifs/par_periode/aout/VTT_aout.tif"))


