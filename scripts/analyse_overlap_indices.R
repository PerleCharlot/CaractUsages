### Titre -------------------------------------
# Nom : analyser indices overlap
# Auteure : Perle Charlot
# Date de création : 20-03-2023
# Dates de modification : -03-2023

### Librairies -------------------------------------
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
### Fonctions -------------------------------------

### Constantes -------------------------------------

# Espace de travail
wd <- getwd()
# Dossier des outputs (dans le git)
output_path <- paste0(wd,"/output/")

#### Autre ####
# liste.mois = c("mai","juin","juillet","aout","septembre")
# df.mois = data.frame(nom_mois = liste.mois, numero_mois = c("05","06","07","08","09"))
# # # Liste dimensions
# # liste.dim =  c("CA","B","PV","CS","D","I")
# # Liste usages
# liste.usages = c("Ni","Lk","Co","Pa","Rp","Vt")

path_analyse_probaD <- paste0(output_path,"/link_G_E_spaces.csv")
path_analyse_overlap <- paste0(output_path,"/link_G_E_spaces_overlap_area.csv")

# type_donnees = "ACP_avec_ponderation"
# algorithme = "glm"
# fit = "2_axes"

### Programme -------------------------------------

df <- fread(path_analyse_probaD, dec=",")
df2 <- fread(path_analyse_overlap, dec=",")

df$Pair <- paste0(df$Pair_Us1, " - ",df$Pair_Us2)

df %>%
  ggplot(aes(x= mean_schoenerD_Gspace_proba , 
             y= mean_schoenerD_Espace_proba , 
             color= Pair_Us1, shape =Pair_Us2)) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1)

df %>%
  ggplot(aes(x=mean_schoenerD_Gspace_proba, y= mean_schoenerD_Espace_proba , color= pair)) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1)

# 1 plot = 5 points (by use)
library(stringr)
P1 <- df[grep("Nesting", df$Pair),] %>%
  ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  labs(x="E space", y = "G space", title = "Nesting", color='Pair')+
  geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
                    ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
  theme(legend.key.height=unit(3, "line"),
        legend.key.width=unit(1, "line"))

P2 <- df[grep("Night", df$Pair),] %>%
  ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  labs(x="E space", y = "G space", title = "Sheep Night Camping", color='Pair')+
  geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
                    ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
  theme(legend.key.height=unit(3, "line"),
        legend.key.width=unit(1, "line"))

P3 <- df[grep("Grazing", df$Pair),] %>%
  ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  labs(x="E space", y = "G space", title = "Sheep Grazing", color='Pair')+
  geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
                    ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
  theme(legend.key.height=unit(3, "line"),
        legend.key.width=unit(1, "line"))

P4 <- df[grep("Hiking", df$Pair),] %>%
  ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  labs(x="E space", y = "G space", title = "Hiking", color='Pair')+
  geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
                    ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
  theme(legend.key.height=unit(3, "line"),
        legend.key.width=unit(1, "line"))

P5 <- df[grep("Bike", df$Pair),] %>%
  ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  labs(x="E space", y = "G space", title = "Mountain Bike", color='Pair')+
  geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
                    ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
  theme(legend.key.height=unit(3, "line"),
        legend.key.width=unit(1, "line"))

P6 <- df[grep("Lek", df$Pair),] %>%
  ggplot(aes(y=mean_schoenerD_Gspace_proba, x= mean_schoenerD_Espace_proba, color= str_wrap(Pair,15))) +
  geom_point() +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  labs(x="E space", y = "G space", title = "Lek", color='Pair')+
  geom_errorbar(aes(ymin= mean_schoenerD_Gspace_proba - sd_schoenerD_Gspace_proba , 
                    ymax=mean_schoenerD_Gspace_proba + sd_schoenerD_Gspace_proba)) +
  theme(legend.key.height=unit(3, "line"),
        legend.key.width=unit(1, "line"))


P = P1 + P2 + P3 +P4 +P5 +P6 + plot_layout(nrow=2, ncol=3)

P + 
  plot_annotation(title = "Schoener D Similarity on probability of occurrence") 


