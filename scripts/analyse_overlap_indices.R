### Titre -------------------------------------
# Nom : analyser indices overlap
# Auteure : Perle Charlot
# Date de création : 20-03-2023
# Dates de modification : 17-04-2023

### Librairies -------------------------------------
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(reshape2)
### Fonctions -------------------------------------

### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")

#### Autre ####
# liste.mois = c("mai","juin","juillet","aout","septembre")
# df.mois = data.frame(nom_mois = liste.mois, numero_mois = c("05","06","07","08","09"))
# Liste dimensions
liste.dim <-  c("CA","B","PV","CS","D","I")
# # Liste usages
# liste.usages = c("Ni","Lk","Co","Pa","Rp","Vt")

type_donnees <- "brute"
algorithme <- "glm"
fit <- "all_simple"

niche_overlap_path <- paste0(output_path,"/niches/",type_donnees, "/niche_overlap/",fit,"/",algorithme,"/")

Espace_path <- paste0(niche_overlap_path,"E_space/")
Gspace_path <- paste0(niche_overlap_path,"G_space/")

# path_analyse_probaD <- paste0(output_path,"/link_G_E_spaces.csv")
# path_analyse_overlap <- paste0(output_path,"/link_G_E_spaces_overlap_area.csv")

### Programme -------------------------------------

### Comparaison entre schoener D dans l'espace Ecologique, sur ACP global et par dimensions ####

# Data G
load(paste0(Gspace_path,"/matrices_schoener_D/matrices_schoener_d.rdata"))
SchoenerD_G <- setNames(melt(matrix_M_std), c('Use1', 'Use2', 'D_G','Month'))
SchoenerD_G <- na.omit(SchoenerD_G)
# Data E global
load(paste0(Espace_path,"global/matrices_schoener_D/matrices_schoener_d.rdata"))
SchoenerD_Eglob <- setNames(melt(matrix_M_std), c('Use1', 'Use2', 'D_Eglob','Month'))
SchoenerD_Eglob <- na.omit(SchoenerD_Eglob)
# Data E dimensions
extract_D_Edim <- function(dim_i){
  #dim_i = "CA"
  load(paste0(Espace_path,"dimensions/matrices_schoener_D/",dim_i,"/matrices_schoener_d.rdata"))
  df <- as.data.frame(setNames(melt(matrix_M_std), c('Use1', 'Use2', 'D_Edim','Month')))
  df <- na.omit(df)
  df$Dimension <- dim_i
  return(df)
}
SchoenerD_Edim <- do.call(rbind,lapply(liste.dim, extract_D_Edim)) 
# Merge E glob and E dim
SchoenerD_E <- merge(SchoenerD_Eglob, SchoenerD_Edim, by=c("Use1","Use2","Month"))


# Comparaison éco(global)/géo
SchoenerD_E_G <- merge(SchoenerD_Eglob, SchoenerD_G, by=c("Use1","Use2","Month"))
SchoenerD_E_G <- subset(SchoenerD_E_G, SchoenerD_E_G$Use1 != SchoenerD_E_G$Use2)

SchoenerD_E_G %>%
  # filter(Use1 == "Nesting" | Use1 == "Hiking") %>%
  # filter(Use2 == "Hiking" | Use2 ==  "Nesting") %>%
  ggplot(aes(x=D_Eglob, y =D_G, colour=Month)) +
  geom_point() +xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1)

SchoenerD_E_G %>%
  filter(Use1 == "Nesting" | Use1 == "Hiking") %>%
    filter(Use2 == "Hiking" | Use2 ==  "Nesting") %>%
  ggplot(aes(x=D_Eglob, y =D_G, colour=Month)) +
  geom_point() +xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1)




# df <- fread(path_analyse_probaD, dec=",")
# df2 <- fread(path_analyse_overlap, dec=",")
# 
# df$Pair <- paste0(df$Pair_Us1, " - ",df$Pair_Us2)
# 
# df %>%
#   ggplot(aes(x= mean_schoenerD_Gspace_proba , 
#              y= mean_schoenerD_Espace_proba , 
#              color= Pair_Us1, shape =Pair_Us2)) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1)
# 
# df %>%
#   ggplot(aes(x=mean_schoenerD_Gspace_proba, y= mean_schoenerD_Espace_proba , color= pair)) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1)
# 
# # 1 plot = 5 points (by use)
# library(stringr)
# P1 <- df[grep("Nesting", df$Pair),] %>%
#   ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
#   labs(x="E space", y = "G space", title = "Nesting", color='Pair')+
#   geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
#                     ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
#   theme(legend.key.height=unit(3, "line"),
#         legend.key.width=unit(1, "line"))
# 
# P2 <- df[grep("Night", df$Pair),] %>%
#   ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
#   labs(x="E space", y = "G space", title = "Sheep Night Camping", color='Pair')+
#   geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
#                     ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
#   theme(legend.key.height=unit(3, "line"),
#         legend.key.width=unit(1, "line"))
# 
# P3 <- df[grep("Grazing", df$Pair),] %>%
#   ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
#   labs(x="E space", y = "G space", title = "Sheep Grazing", color='Pair')+
#   geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
#                     ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
#   theme(legend.key.height=unit(3, "line"),
#         legend.key.width=unit(1, "line"))
# 
# P4 <- df[grep("Hiking", df$Pair),] %>%
#   ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
#   labs(x="E space", y = "G space", title = "Hiking", color='Pair')+
#   geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
#                     ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
#   theme(legend.key.height=unit(3, "line"),
#         legend.key.width=unit(1, "line"))
# 
# P5 <- df[grep("Bike", df$Pair),] %>%
#   ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
#   labs(x="E space", y = "G space", title = "Mountain Bike", color='Pair')+
#   geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
#                     ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
#   theme(legend.key.height=unit(3, "line"),
#         legend.key.width=unit(1, "line"))
# 
# P6 <- df[grep("Lek", df$Pair),] %>%
#   ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
#   geom_point() +
#   xlim(0,1) + ylim(0,1) +
#   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
#   labs(x="E space", y = "G space", title = "Lek", color='Pair')+
#   geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
#                     ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
#   theme(legend.key.height=unit(3, "line"),
#         legend.key.width=unit(1, "line"))
# 
# 
# P = P1 + P2 + P3 +P4 +P5 +P6 + plot_layout(nrow=2, ncol=3)
# 
# P + 
#   plot_annotation(title = "Schoener D Similarity on probability of occurrence") 
# 
# ###############################
# 
# chemin <- paste0(wd,"/output/niches/ACP_avec_ponderation/niche_overlap/2_axes/glm/overlap_metrics/test_comparaison_D.csv")
# path_save <- paste0(wd,"/output/niches/exploration_comparaison_EG/")
# table_D <- fread(chemin)
# str(table_D)
# 
# 
# table_D$Mois <- factor(table_D$Mois,      # Reordering group factor levels
#                                    levels = c("mai","juin","juillet","aout","septembre"))
# 
# # comparaison chevauchement G/E, 
# # en groupant par type de paire d'usage
# # (Activité humaine - activité humaine et faune - activité humaine)
# P <- table_D %>%
#   ggplot(aes(x=D_geo,y=D_eco,colour=faune_vs_hum)) +
#   geom_point(size=3)+
#   xlim(0,1) + ylim(0,1)+
#   geom_abline(intercept = 0, slope = 1, linetype="dashed")+
#   labs(x="Geographic Overlap", y = "Ecological Overlap", colour="Type of uses pair")
# # le détail mois par mois
# P1 <- table_D %>%
#   ggplot(aes(x=D_geo,y=D_eco,colour=faune_vs_hum)) +
#   geom_point(size=3)+
#   xlim(0,1) + ylim(0,1)+
#   geom_abline(intercept = 0, slope = 1, linetype="dashed")+
#   facet_wrap(~Mois)+
#   theme(legend.position = "none")+
#   labs(x="Geographic Overlap", y = "Ecological Overlap")
# 
# png(file = paste0(path_save,"/comparison_overlap_EGspaces_typeofPairUses.png"),width=1500, height=800)
# plot(P + P1 + plot_layout(ncol=2,nrow=1))
# dev.off()
# 
# # comparaison chevauchement G/E, coloré par mois
# Pglob <- table_D %>%
#   ggplot(aes(x=D_geo,y=D_eco,colour=Mois)) +
#   geom_point(size=3)+
#   xlim(0,1) + ylim(0,1)+
#   geom_abline(intercept = 0, slope = 1, linetype="dashed")+
#   labs(x="Geographic Overlap", y = "Ecological Overlap")
# 
# png(file = paste0(path_save,"/comparison_overlap_EGspaces_colorMonth.png"),width=1500, height=800)
# plot(Pglob)
# dev.off()
# # comparaison chevauchement G/E, coloré par mois, séparé par paire
# P2 <- table_D %>%
#   ggplot(aes(x=D_geo,y=D_eco,colour=Mois))+
#   geom_point(size=3) +
#   xlim(0,1) + ylim(0,1)+
#   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
#   facet_wrap(~Paire)+
#   labs(x="Geographic Overlap", y = "Ecological Overlap")
# png(file = paste0(path_save,"/comparison_overlap_EGspaces_colorMonth_byPair.png"),width=1500, height=800)
# plot(P2)
# dev.off()
