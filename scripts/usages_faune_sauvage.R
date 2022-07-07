### Titre -------------------------------------
# Nom : Cartographie usages Faune Sauvage
# Auteure : Perle Charlot
# Date de création : 15-06-2022
# Dates de modification : 25-06-2022

### Librairies -------------------------------------
library(sf)
library(fasterize)
library(raster)
### Fonctions -------------------------------------

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

# Couche nidification TLY
path_nidif_TLY <- paste0(dos_var_sp, "Usages/faune_sauvage/donnes_galliformes_FDCI/Zone de nichées TLY.shp")
# Couche zone de parade TLY
path_chant_TLY <- paste0(dos_var_sp, "Usages/faune_sauvage/donnes_galliformes_FDCI/Secteur avec place de chant répertorié.shp")

# MNT 25m CALé sur le grille de REFERENCE
path_MNT <- paste0(dos_var_sp ,"/Milieux/IGN/mnt_25m_belledonne_cale.tif")

#### Tables ####

### Programme -------------------------------------

#### Nidification TLY ####
nidif_TLY <- st_read(path_nidif_TLY)
raster_nidif_TLY <- fasterize(nidif_TLY, raster(path_MNT), background = 0)
# 0 pour absence (pseudo absence), 1 pour présence
plot(raster_nidif_TLY)
writeRaster(raster_nidif_TLY, paste0(output_path,"/usages_faune_sauvage/nidification_TLY.tif"),overwrite=TRUE)

#### Parade TLY ####
chant_TLY <- st_read(path_chant_TLY)
raster_chant_TLY <- fasterize(chant_TLY, raster(path_MNT), background = 0)
# 0 pour absence (pseudo absence), 1 pour présence
plot(raster_chant_TLY)
writeRaster(raster_chant_TLY, paste0(output_path,"/usages_faune_sauvage/parade_TLY.tif"), overwrite=TRUE)
