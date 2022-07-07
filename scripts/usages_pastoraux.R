### Titre -------------------------------------
# Nom : Cartographie usages pastoraux
# Auteure : Perle Charlot
# Date de création : 25-06-2022
# Dates de modification : 03-07-2022

### Librairies -------------------------------------
library(sf)
library(fasterize)
library(raster)
### Fonctions -------------------------------------

# Fonction qui rastérise des polygones d'usages, par période donnée
RasterizePolyUsage <- function(mois, nom_usage, chemin_spat_usage, raster_de_ref){
  
  # # TEST
  # mois = "juin"
  # nom_usage = "couchade"
  # chemin_spat_usage = list(path_couchade_Pra_67, path_couchade_Gab_679)
  # raster_de_ref = rast_ref

  list_sp_usage <- lapply(chemin_spat_usage, st_read)
  sp_usage <- st_as_sf(data.table::rbindlist(list_sp_usage, fill=TRUE))

  # 0 pour absence (pseudo absence), 1 pour présence
  rast_sp_usage <- fasterize(sp_usage, raster_de_ref, background = 0)
  #plot(rast_sp_usage)
  
  writeRaster(rast_sp_usage, 
              paste0(output_path,"/usages_pastoraux/par_periode/",mois,"/",nom_usage,".tif"),
              overwrite=TRUE)
  
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

# Couches usages pastoraux (paturage et couchade)
path_paturage_juin_Pra <- paste0(dos_var_sp, "Usages/pastoralisme/La Pra/paturage_juin_la_pra.gpkg")
path_paturage_juin_Gab <- paste0(dos_var_sp, "Usages/pastoralisme/Gaboureaux/paturage_juin.gpkg")
path_paturage_juill_Pra <- paste0(dos_var_sp, "Usages/pastoralisme/La Pra/paturage_juillet_la_pra_sans_couchade.gpkg")
path_paturage_juill_Gab <- paste0(dos_var_sp, "Usages/pastoralisme/Gaboureaux/paturage_juillet.gpkg")
path_paturage_aout_Pra <- paste0(dos_var_sp, "Usages/pastoralisme/La Pra/paturage_aout_la_pra_sans_couchade.gpkg")
path_paturage_aout_Gab <- paste0(dos_var_sp, "Usages/pastoralisme/Gaboureaux/paturage_aout.gpkg")
path_paturage_sept_Pra <- paste0(dos_var_sp, "Usages/pastoralisme/La Pra/paturage_septembre_la_pra_sans_couchade.gpkg")
path_paturage_sept_Gab <- paste0(dos_var_sp, "Usages/pastoralisme/Gaboureaux/paturage_septembre.gpkg")

path_couchade_Pra_789 <- paste0(dos_var_sp, "Usages/pastoralisme/La Pra/couchade_juillet_aout_septembre.gpkg")
path_couchade_Pra_67 <- paste0(dos_var_sp, "Usages/pastoralisme/La Pra/couchade_juin_juillet_v2.gpkg")
path_couchade_Gab_679 <- paste0(dos_var_sp, "Usages/pastoralisme/Gaboureaux/couchade_juin_juillet_septembre.gpkg")
path_couchade_Gab_8<- paste0(dos_var_sp, "Usages/pastoralisme/Gaboureaux/couchade_aout.gpkg")

#### Tables ####

### Programme -------------------------------------

rast_ref <- raster(path_MNT)

#### Pâturage ####

RasterizePolyUsage("juin",
                   "paturage",
                   list(path_paturage_juin_Pra,path_paturage_juin_Gab),
                   rast_ref)

RasterizePolyUsage("juillet",
                   "paturage",
                   list(path_paturage_juill_Pra,path_paturage_juill_Gab),
                   rast_ref)

RasterizePolyUsage("aout",
                   "paturage",
                   list(path_paturage_aout_Pra,path_paturage_aout_Gab),
                   rast_ref)

RasterizePolyUsage("septembre",
                   "paturage",
                   list(path_paturage_sept_Pra,path_paturage_sept_Gab),
                   rast_ref)


#### Couchade ####

RasterizePolyUsage("juin",
                   "couchade",
                   list(path_couchade_Pra_67, path_couchade_Gab_679),
                   rast_ref)
RasterizePolyUsage("juillet",
                   "couchade",
                   list(path_couchade_Pra_67,path_couchade_Pra_789, path_couchade_Gab_679),
                   rast_ref)
RasterizePolyUsage("aout",
                   "couchade",
                   list(path_couchade_Pra_789 ,path_couchade_Gab_8),
                   rast_ref)
RasterizePolyUsage("septembre",
                   "couchade",
                   list(path_couchade_Pra_789 ,path_couchade_Gab_679),
                   rast_ref)

