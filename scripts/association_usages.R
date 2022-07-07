### Titre -------------------------------------
# Nom : Interactions spatiales entre usages
# Auteure : Perle Charlot
# Date de création : 05-07-2022
# Dates de modification : 06-07-2022

### Librairies -------------------------------------

library(raster)
library(data.table)
library(sf)
library(magick)
library(corrplot)

# library(fasterize)

### Fonctions -------------------------------------

# Transforme stack d'usages en dataframe, pour un mois donné
DfUsages <- function(mois){
  # # TEST
  # mois = "juin"
  
  list_usages = list.files(paste0(output_path,"/par_periode/",mois),recursive=TRUE, ".tif", full.names=TRUE)
  usages = stack(list_usages)
  usages_masked <- raster::mask(usages, limiteN2000)
  df_usages_masked = as.data.frame(as.data.table(usages_masked[]))
  df_usages_masked = na.omit(df_usages_masked)
  
  # Statistiques spatiales, par usage
  # pixels de 25*25m  = 625m²
  f = function(colonne){
    (sum(colonne) * 625) / 10000
  }

  superficie_usage = apply(df_usages_masked, 2, f)
  aire_N2000 = (dim(df_usages_masked)[1]*625)/10000
  pourcent_superficieN2000 = (superficie_usage/aire_N2000)*100
  
  df_stat_usage = data.frame(superficie_usage, pourcent_superficieN2000)
  write.csv(df_stat_usage, paste0(output_path,"/par_periode/",mois,"/statistiques_spatiales_usages_",mois,".csv"))

  
  cor <- cor(df_usages_masked, method = "spearman")
 
  #"square", "ellipse", "number", "shade", "color", "pie""circle"
  
  png(file=paste0(output_path,"/correlation/plot_spearman_",mois,".png"),
      width=800, height=800)
  corrplot(cor, method = "circle", type="lower", title = mois,tl.srt=60,mar=c(0,1,1,0))
  dev.off()
  
  return(df_usages_masked)
}

# Transforme stack d'usages en dataframe, pour un mois donné
CarteMultiUsages <- function(mois){
  # # TEST
  # mois = "juin"
  
  list_usages = list.files(paste0(output_path,"/par_periode/",mois),recursive=TRUE, ".tif", full.names=TRUE)
  usages = stack(list_usages)
  usages_masked <- mask(usages, limiteN2000)

  return(usages_masked)
}

### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")
input_path <- paste0(wd,"/input/")
# Projection Lambert 93 (EPSG : 2154)
EPSG_2154 =  "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs "

#### Données spatiales ####

# Dossier des variables spatiales & chemins des fichiers
dos_var_sp <- "C:/Users/perle.charlot/Documents/PhD/DATA/Variables_spatiales_Belledonne/"

limiteN2000 <- paste0(dos_var_sp, "/limites_etude/cembraie_N2000_limites.gpkg")

#### Tables ####

### Programme -------------------------------------

# Charger les cartes des usages

# Transformer en tableau de 0(absence) et 1 (présence)

# Répéter l'opération pour chaque mois

limiteN2000 <- st_read(limiteN2000)

dfUsages_04 <- DfUsages("avril")
dfUsages_05 <- DfUsages("mai")
dfUsages_06 <- DfUsages("juin")
dfUsages_07 <- DfUsages("juillet")
dfUsages_08 <- DfUsages("aout")
dfUsages_09 <- DfUsages("septembre")


#### Coefficient d'association de Yule #### 
library(psych)

# # TEster le comportement du coefficient
# Yule(matrix(nrow=2, ncol=2 ,c(0,100,100,0))) # -1
# Yule(matrix(nrow=2, ncol=2 ,c(100,100,100,100))) # 0
# Yule(matrix(nrow=2, ncol=2 ,c(100,0,0,100))) # 1

Yule(table(dfUsages_06$paturage, dfUsages_06$nidification_TLY_juin)) # 0.81 -> synergie entre usages
Yule(matrix(nrow=2, ncol=2 ,c(0,1939,3078,1578)))

Yule(table(dfUsages_06$paturage, dfUsages_06$randonnee_pedestre_juin))

#### Indice et distance de Jaccard #### 

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard(dfUsages_06$paturage, dfUsages_06$nidification_TLY_juin) # ~0 : pas de similarité entre les deux vecteurs

#### Odds ratio #### 
library(questionr)
table(dfUsages_06$paturage, dfUsages_06$nidification_TLY_juin)
odds.ratio(table(dfUsages_06$paturage, dfUsages_06$nidification_TLY_juin))
# OR = 9.5 >1 -> synergies entre usages (la présence de l'un favorise celle de l'autre)

#### Local correlation (raster::) #### 
d = corLocal(b_masked$paturage, b_masked$randonnee_pedestre_juin, method="spearman", test=TRUE)
plot(d, colNa="black")
# only cells where the p-value < 0.1
xm <- mask(d[[1]], d[[2]] < 0.1, maskvalue=FALSE)
plot(xm,colNa="black")
mean(values(xm), na.omit=TRUE)

tt = values(xm)
tt = na.omit(tt)
mean(tt)

#### Association between paired samples ####

cor.test(dfUsages_06$paturage, dfUsages_06$nidification_TLY_juin, method = "kendall")
# tau = 0.326 (pvalue < 0.05) -> 
# plutôt pas d'asociation, mais similiaire (tau proche de 0 mais positif)
cor.test(dfUsages_06$paturage, dfUsages_06$nidification_TLY_juin, method = "spearman")
# rho = 0.326 IDEM
cor.test(dfUsages_06$paturage, dfUsages_06$nidification_TLY_juin, method = "pearson")
# cor = 0.326 IDEM


#### Superposition des usages ####
liste.mois = c("avril","mai","juin","juillet","aout","septembre")
liste.rast.usage = lapply(liste.mois, CarteMultiUsages)
liste.rast.multiusage =lapply(liste.rast.usage, function(r)sum(r))
rast.multiusage = stack(liste.rast.multiusage)
names(rast.multiusage) = paste0("sumUsage_0",seq(4,9,1),liste.mois)
writeRaster(rast.multiusage, bylayer=TRUE, 
            paste0(output_path,"/multiusage/.tif"),suffix=names(rast.multiusage))

# TODO : Statistiques multiusages au cours du temps

f2 = function(r){  return(table(values(r)))}
lapply(liste.rast.multiusage, f2)

# Création d'un gif pour visualiser les zones multiusage au cours des mois
imgs <- list.files(input_path, full.names = TRUE, '.png')
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, delay = 100)
image_write(image = img_animated,
            path = paste0(output_path,"/multiusage/multiusages.gif"))

