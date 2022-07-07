### Titre -------------------------------------
# Nom : Extraction traces gpx depuis web
# Auteure : Perle Charlot
# Date de création : 24-06-2022
# Dates de modification : 03-07-2022

### Librairies -------------------------------------
library(sf)
library(purrr)
library(dplyr)
library(data.table)
library(stringi)
### Fonctions -------------------------------------

# fonction qui permet de lire un gpx (de Singletrack) et le mettre propre
ReadBunchGPX <- function(x) {
  
  # # # # TEST
  # x = liste_gpx[[1]]
  
  # # # # # Infos sur les couches dans le gpx
  # st_layers(x)

  format = stri_sub(x,-1)
  
  if(format == "g"){
    cat("\nAssemblage des traces gpkg.")
    gpx_file <- st_read(x)
  } 
  if(format == "l"){
    cat("\nAssemblage des traces kml.")
    gpx_file <- st_read(x)
    # Ne garder que les linestring s'il y a aussi des points (d'entrée, sortie, nuitée)
    index = which(st_geometry_type(gpx_file) == "LINESTRING" | 
            st_geometry_type(gpx_file) == "MULTILINESTRING")
    gpx_file = gpx_file[index,]
    gpx_file <- gpx_file %>% select(1)
    gpx_file <- st_cast(gpx_file,"MULTILINESTRING")
  }  
  
  if(format == "x"){
    cat("\nAssemblage des traces gpx.")
    gpx_file <- st_read(x,layer='tracks')
    
    # Si 'tracks' vide, aller piocher dans 'routes'
    if(dim(gpx_file)[1] == 0){
      gpx_file <- st_read(x,layer='routes')
    }
    gpx_file <- gpx_file %>% select(1)
    gpx_file <- st_cast(gpx_file,"MULTILINESTRING")
    
  }
  
  return(gpx_file)
}

# Fonction qui combine toutes traces gpx en un unique objet spatial
CombineGPX <- function(path_dos_gpx, pratique, source=NULL){
  
  # # # # TEST
  # path_dos_gpx = path_dos_gpkg_randonnee
  # pratique = "randonnee"
  # source = NULL
  
  if(is.null(source)){
    cat("Assemblage des traces des différentes sources.")
    liste_gpx <- as.list(list.files(path_dos_gpx, ".gpkg", full.names = TRUE))
    liste_sp_gpx <- map(liste_gpx,ReadBunchGPX)
    # fusionner toutes les traces
    liste_sp_gpx <- st_as_sf(data.table::rbindlist(liste_sp_gpx, fill=TRUE))
    # Projeter bien
    liste_sp <- st_transform(liste_sp_gpx,EPSG_2154)
  } else{
    cat("Assemblage des traces du site ",source,".")
    liste_gpx <- as.list(list.files(path_dos_gpx, ".gpx", full.names = TRUE))
    liste_kml <- as.list(list.files(path_dos_gpx, ".kml", full.names = TRUE))
    
    if(length(liste_kml)>0){
      
      liste_sp_kml <- map(liste_kml,ReadBunchGPX)
      liste_sp_kml <- st_as_sf(data.table::rbindlist(liste_sp_kml))
      liste_sp_kml <- st_zm(liste_sp_kml)
      names(liste_sp_kml)[1] = "name"
      liste_sp_kml <- st_transform(liste_sp_kml,EPSG_2154)
    }
    
    liste_sp_gpx <- map(liste_gpx,ReadBunchGPX)
    # fusionner toutes les traces
    liste_sp_gpx <- st_as_sf(data.table::rbindlist(liste_sp_gpx, fill=TRUE))
    # Projeter bien
    liste_sp_gpx <- st_transform(liste_sp_gpx,EPSG_2154)
    
    if(length(liste_kml)>0){
      # Rassembler kml et gpx si existent
      liste_sp <- rbind(liste_sp_gpx, liste_sp_kml)
    } else {liste_sp <- liste_sp_gpx}
  }

  # Sauvegarder
  st_write(liste_sp, paste0(output_path, "/usages_recreatifs/",pratique,"/",
                                pratique,"_",source,".gpkg"))
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
# vecteur emprise carrée
path_emprise <- paste0(dos_var_sp ,"/limites_etude/emprise.gpkg")
# MNT 25m CALé sur le grille de REFERENCE
path_MNT <- paste0(dos_var_sp ,"/Milieux/IGN/mnt_25m_belledonne_cale.tif")

# Dossier des traces gpx de VTT
path_dos_gpx_singletrack = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/Singletrack/")
path_topo_VTTour <- paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/VTTour/liens.txt")
path_dos_gpx_vttour = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/VTTour/")
path_topo_VisuGPX <- paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/VisuGPX/liens_VTT.txt")
path_dos_gpx_VisuGPX = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/VisuGPX/")
path_dos_gpx_station = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/Chamrousse_station/")
# Dossier général
path_dos_gpkg_VTT <- paste0(output_path, "/usages_recreatifs/VTT/")
path_dos_gpkg_VTTRMK <- paste0(output_path, "/usages_recreatifs/VTT_RMK/")

# Dossier des traces pour randonnée pédestre
path_dos_gpx_camptocamp = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/randonnee/camptocamp/")
path_dos_gpx_visorando = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/randonnee/visorando/") 
path_dos_gpx_altituderando = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/randonnee/altituderando/")
path_dos_gpx_station_rando = paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/randonnee/chamrousse_station/")
# Dossier général
path_dos_gpkg_randonnee <- paste0(output_path, "/usages_recreatifs/randonnee/")
path_dos_gpkg_randonneeRMK <- paste0(output_path, "/usages_recreatifs/randonnee_RMK/")

#### Tables ####

### Programme -------------------------------------

#### VTT ####

### Singletrack ###
CombineGPX(path_dos_gpx_singletrack ,"VTT","singletrack")
# LIE REMONTE MECANIQUE : concerne très peu d'itinéraire (1)
CombineGPX(paste0(path_dos_gpx_singletrack,"/dpdt_remontee_mecanique/"),
           "VTT_RMK","singletrack")

### VTTour ###
liste_topo <-fread(path_topo_VTTour, header=FALSE)
liste_url_gpx <- paste0("http://www.vttour.fr/topos/gpx/download.php?id=",sub(".html","",liste_topo$V2))
liste_destfile <- paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/VTTour/vttour_",
                   sub(".html","",liste_topo$V2),".gpx")
# Téléchargement de toutes les traces gpx présentes sur la liste de liens
map2(liste_url_gpx, liste_destfile,download.file)
# Combinaison des traces en un seul fichier
CombineGPX(path_dos_gpx_vttour ,"VTT","vttour")

### VisuGPX ###
liste_topo <-fread(path_topo_VisuGPX, header=FALSE, sep="/")
liste_url_gpx <- paste0("http://www.visugpx.com/download.php?id=",liste_topo$V4)
liste_destfile <- paste0(dos_var_sp ,"/Usages/Pratiques_récréatives/VTT/VisuGPX/visugpx_",
                         liste_topo$V4,".gpx")
# Téléchargement de toutes les traces gpx présentes sur la liste de liens
map2(liste_url_gpx, liste_destfile,download.file)
# Combinaison des traces en un seul fichier
CombineGPX(path_dos_gpx_VisuGPX ,"VTT","visugpx")
# LIE REMONTE MECANIQUE : concerne très peu d'itinéraire (1)
CombineGPX(paste0(path_dos_gpx_VisuGPX,"/dpdt_remontee_mecanique/"),
           "VTT_RMK","visugpx")

### Site Chamrousse ###
CombineGPX(path_dos_gpx_station ,"VTT","station")
# LIE REMONTE MECANIQUE : concerne très peu d'itinéraire (1)
CombineGPX(paste0(path_dos_gpx_station,"/dpdt_remontee_mecanique/"),
           "VTT_RMK","station")


### Combiner tous les tracés ###
CombineGPX(path_dos_gpx = path_dos_gpkg_VTT , pratique = "VTT", source=NULL)
CombineGPX(path_dos_gpx = path_dos_gpkg_VTTRMK , pratique = "VTT_RMK", source=NULL)


# TODO : regarder selon altitude si c'est descente ou montée, et enregistrer séparément
# Dans singletrack, on voit le D+/D- et on en déduit si la trace était montée+descente ou juste descente (si pas une boucle)
# idem dans VisuGPx
# Mais le hic c'est que je vois pas comment le faire autrement qu'à la main..

# TODO : réfléchir si appliquer méthode de Byczek pour aire d'influence

## Strava : VTT ###

# library(curl)
# # page API de skitour
# # https://skitour.fr/api/
# 
# # Notice package curl
# # https://cran.r-project.org/web/packages/curl/vignettes/intro.html#Setting_handle_options
# 
# # Notice package httr
# # https://cran.r-project.org/web/packages/httr/vignettes/quickstart.html 
# 
# 
# # API pour les sorties
# # cle: yYOlbMcps8rOuUCdVDTwbawLXO26IHLM
# 
# req <- curl_fetch_memory("https://skitour.fr/api/sorties?a=2022&s=84")
# str(req)
# parse_headers(req$headers)
# jsonlite::prettify(rawToChar(req$content))
# 
# library(httr)
# 
# url      <- "https://skitour.fr/api/sorties?"
# query    <- paste0("a=2022&s=84") # a = année, s = sommet
# full_url <- paste0(url, query)
# key      <- "yYOlbMcps8rOuUCdVDTwbawLXO26IHLM"
# key_test      <- "cle: yYOlbMcps8rOuUCdVDTwbawLXO26IHLM"
# res      <- GET(full_url, add_headers(Authorization = paste("Bearer", key)))
# res      <- GET(full_url, add_headers(Authorization = paste("Bearer", key_test)))
# res      <- GET(full_url, add_headers(cle = paste(key)))
# str(res)
# res_text <- content(res, "text")
# jsonlite::prettify(rawToChar(res$content))
# 
# # ça ça marche étrangement
# res_sortie <- GET("https://skitour.fr/sorties/571", add_headers(cle = paste(key)), content_type(".gpx"))
# str(res_sortie)
# 
# content(res_sortie, "text")
# a = res_sortie$all_headers
# 
# jsonlite::prettify(rawToChar(res_sortie$content))
# 
# fromJSON(file = "input.json")

#### Randonnée pédestre ####

### CampToCamp ###
# Combinaison des traces en un seul fichier
CombineGPX(path_dos_gpx_camptocamp ,"randonnee","camptocamp")

### Visorando ###
# Combinaison des traces en un seul fichier
CombineGPX(path_dos_gpx_visorando ,"randonnee","visorando")
# LIE REMONTE MECANIQUE : concerne très peu d'itinéraire (1)
CombineGPX(paste0(path_dos_gpx_visorando,"/dpdt_remontee_mecanique/"),
           "randonnee_RMK","visorando")

### Altituderando ###
# Combinaison des traces en un seul fichier
CombineGPX(path_dos_gpx_altituderando ,"randonnee","altituderando")

### Station Chamrousse ###
# Combinaison des traces en un seul fichier
CombineGPX(path_dos_gpx_station_rando ,"randonnee","station_Chamrousse")
# Différencier randonnée selon la saison (ouvertiure remontée mécanique)
# quand précisé dans le topo de la randonnée
# concerne très peu d'itinéraire (4)
CombineGPX(paste0(path_dos_gpx_station_rando,"/dpdt_remontee_mecanique/"),
           "randonnee_RMK","station_Chamrousse")

### Combiner tous les tracés ###
CombineGPX(path_dos_gpx = path_dos_gpkg_randonnee , pratique = "randonnee", source=NULL)
CombineGPX(path_dos_gpx = path_dos_gpkg_randonneeRMK , pratique = "randonnee_RMK", source=NULL)


