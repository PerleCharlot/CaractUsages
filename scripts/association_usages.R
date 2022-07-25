### Titre -------------------------------------
# Nom : Interactions spatiales entre usages
# Auteure : Perle Charlot
# Date de création : 05-07-2022
# Dates de modification : 25-07-2022

### Librairies -------------------------------------

library(raster)
library(data.table)
library(sf)
library(magick)
library(corrplot)
library(ggplot2)
library(dplyr)
library(tidyverse)

# library(fasterize)

### Fonctions -------------------------------------

# Transforme stack d'usages en dataframe, pour un mois donné
DfUsages <- function(mois){
  # # TEST
  # mois = "juin"
  
  list_usages = list.files(paste0(output_path,"/par_periode/",mois),recursive=TRUE, ".tif$", full.names=TRUE)
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
  
  list_usages = list.files(paste0(output_path,"/par_periode/",mois),recursive=TRUE, ".tif$", full.names=TRUE)
  usages <- stack(list_usages)
  usages_masked <- raster::mask(usages, limiteN2000)

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

limiteN2000 <- st_read(limiteN2000)

# dfUsages_04 <- DfUsages("avril")
# dfUsages_05 <- DfUsages("mai")
# dfUsages_06 <- DfUsages("juin")
# dfUsages_07 <- DfUsages("juillet")
# dfUsages_08 <- DfUsages("aout")
# dfUsages_09 <- DfUsages("septembre")


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
liste.rast.multiusage = lapply(liste.rast.usage, function(r)sum(r))
rast.multiusage = stack(liste.rast.multiusage)
names(rast.multiusage) = paste0("sumUsage_0",seq(4,9,1),liste.mois)
writeRaster(rast.multiusage, bylayer=TRUE, 
            paste0(output_path,"/multiusage/.tif"),suffix=names(rast.multiusage))

plot(rast.multiusage)

dt_stack <- as.data.frame(data.table(as.data.frame(rast.multiusage)))
dt_stack <- cbind(dt_stack,coordinates(rast.multiusage))
dt_stack <- dt_stack[complete.cases(dt_stack),]

write.csv(dt_stack, paste0(output_path,"/multiusage/table_points_multiusage.csv"))

# TODO : Statistiques multiusages au cours du temps

stats_multiusage = lapply(liste.rast.multiusage, function(r){  return(table(values(r)))})
max_len <- max(sapply(stats_multiusage,length))
stats_multiusage = as.data.frame(do.call("rbind",lapply(stats_multiusage,function(x) x[seq(max_len)])))
names(stats_multiusage) = c(0:5)
stats_multiusage[is.na(stats_multiusage)] <- 0
# Aire relative à l'aire totale
aire_N2000 = (sum(stats_multiusage[1,1:6])*25*25)/10000
stats_multiusage_rel = (((stats_multiusage*25*25)/10000) /aire_N2000 )*100
stats_multiusage_rel$mois = paste0("0",4:9,"_",liste.mois)
stats_multiusage_rel <- stats_multiusage_rel %>% 
  group_by(mois) %>% mutate(sum_multiusage = sum(c_across('2':'5')))

write.csv(stats_multiusage_rel, paste0(output_path,"/multiusage/pourcent_multiusage.csv"))

# stats_multiusage$mois = paste0("0",4:9,"_",liste.mois)
# df_stats_multiusage <- stats_multiusage %>% 
#   pivot_longer(cols=1:6, names_to = "somme_usages", values_to = "nb_pix")
# df_stats_multiusage %>%
#   subset(somme_usages >0) %>%
#   ggplot(aes(x=somme_usages,y=(nb_pix*25^2)/10000, fill=mois)) +
#   geom_col() +
#   facet_wrap(~mois)+
#   labs(y="Superficie (ha)")

df_stats_multiusage_rel <- stats_multiusage_rel %>% 
  pivot_longer(cols=1:6, names_to = "somme_usages", values_to = "pourcent_aire_N2000")

df_stats_multiusage_rel %>%
  subset(somme_usages >0) %>%
  ggplot(aes(x=somme_usages,y=pourcent_aire_N2000, fill=somme_usages)) +
  geom_col() +
  facet_wrap(~mois)+
  labs(y="Superficie de la zone Natura2000 (en %)", fill="Nombre d'usages")+
  theme(axis.title.x = element_blank())+
  geom_text(aes(label = round(pourcent_aire_N2000,2)), vjust = -0.5)


df_stats_multiusage_rel %>%
  subset(somme_usages >1) %>%
  ggplot(aes(x=somme_usages,y=pourcent_aire_N2000, fill=somme_usages)) +
  geom_col() +
  facet_wrap(~mois)+
  labs(y="Superficie de la zone Natura2000 (en %)", fill="Nombre d'usages\nsimultanés")+
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size = 34),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 28),
        legend.text = element_text(size=22),
        legend.title = element_text(size=24)
        
       # ,legend.position = "none"
        )+
  geom_text(aes(label = round(pourcent_aire_N2000,2)), vjust = 0.5, size=12)+
  scale_fill_manual(values=c("springgreen3",
                             "lightgoldenrodyellow",
                             "sandybrown","red3"))


p1 = df_stats_multiusage_rel %>%
  subset(somme_usages >1) %>%
  ggplot(aes(x=somme_usages,y=pourcent_aire_N2000, fill=somme_usages)) +
  geom_col() +
  facet_wrap(~mois)+
  labs(y="Superficie de la zone Natura2000 (en %)", fill="Nombre d'usages")+
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size = 34),
        axis.text = element_text(size = 28),
        axis.title = element_text(size = 28),
        legend.position = "none")+
  geom_text(aes(label = round(pourcent_aire_N2000,2)), vjust = 0.5, size=12)+
  scale_fill_manual(values=c("springgreen3",
                             "lightgoldenrodyellow",
                             "sandybrown","red3"))

p1
png(file=paste0(output_path,"/multiusage/pourcent_multiusage.png"), 
    width=1400, height=800)
p1
dev.off()

# Création d'un gif pour visualiser les zones multiusage au cours des mois
imgs <- list.files(input_path, full.names = TRUE, '.png')
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, delay = 100)
image_write(image = img_animated,
            path = paste0(output_path,"/multiusage/multiusages.gif"))

# TODO: Quels sont les usages qui se superposent le plus ?
liste.rast.usage

liste.rast.usage = map2(liste.rast.usage, liste.mois, function(r,lm) {
  names(r) = paste0(names(r),"_",lm)
  return(r)
  })

liste.dt.rast = lapply(liste.rast.usage, function(r){
  a = as.data.frame(data.table(as.data.frame(r)))
  a = a[complete.cases(a),]
  return(a)
  })

# TEst sur JUIN
bb = as.data.frame(liste.dt.rast[2])
# REcodage
# C = couchade, N = nidification, PR = parade, PT = paturage, R = randonnée, V=VTT
bbb = data.frame(
couchade_juin = ifelse(bb$couchade_juin == 1, "C",0),
nidification_TLY_juin = ifelse(bb$nidification_TLY_juin == 1, "N",0),
parade_TLY_juin = ifelse(bb$parade_TLY_juin == 1, "P",0),
paturage_juin = ifelse(bb$paturage_juin == 1, "A",0),
randonnee_pedestre_juin = ifelse(bb$randonnee_pedestre_juin == 1, "R",0),
VTT_juin = ifelse(bb$VTT_juin == 1, "V",0),
total = apply(bb,1,sum)
)
# garder que pixels avec multiusage
b4 = bbb %>% subset(total > 1) 
b4$quels_us = apply(b4, 1,function(df) paste0(df[], collapse=""))

qui_multi = table(gsub('[0-9]',"",b4$quels_us))

ok = as.data.frame(sort(qui_multi))
# Pour le mois de juin, le multiusage a lieu le plus entre
# 2 usages =
# nidification et parade (1524 pixels)
# randonnée et VTT (1520 pixels)
# nidification et paturage (1173 pixels)
# 3 usages =
# NPR (250 pixels) : Nidification Parade Randonnée
# NPA (171 pixels) : Nidification Parade Paturage
# NRV (82 pixels) : Nidification Randonnée VTT
# 4 usages =
# NPRV (164 pixels) : Nidification Parade Randonnée VTT
# NARV (35 pixels) : Nidification Paturage Randonnée VTT
# NPAR (25 pixels) : Nidification Parade Paturage Randonnée

# surface de chaque usage

(apply(bb,2,sum)*625)/10000
