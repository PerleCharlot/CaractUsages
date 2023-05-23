### Titre -------------------------------------
# Nom : exploration co-occurrence
# Auteure : Perle Charlot
# Date de création : 23-05-2023
# Dates de modification : -05-2023

### Librairies -------------------------------------
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(raster)
library(viridis)

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


uses.proba <- c("Ni_proba","Pa_proba","Rp_proba",
                "Lk_proba","Co_proba", "Vt_proba")
names(uses.proba) <- c("Nesting","Sheep Grazing","Hiking",
                       "Lek","Sheep Night Camping" ,"Mountain Bike")

label_test <- c("Nesting","Sheep Grazing","Hiking",
                "Lek","Sheep Night Camping" ,"Mountain Bike")
names(label_test) <- c("Ni_proba","Pa_proba","Rp_proba",
                       "Lk_proba","Co_proba", "Vt_proba")

### Programme -------------------------------------


# Carte de co-occurrence = multiplication probabilité prédiction

example_month <- "June"


load(paste0(niche_overlap_path,"/tables_data.rdata"))
#load df_PCA_dim + df_PCA_glob + df_proba_uses

# Mise en forme
df_proba_PCA_glob <- merge(df_proba_uses,df_PCA_glob, by=c("x","y","Month"),all=T)
df_proba_PCA_glob$Month <- factor(df_proba_PCA_glob$Month,      
                                  levels = c("May","June","July",
                                             "August","September"))
# filter to get for one month
df_proba_PCA_glob_sub <- na.omit(df_proba_PCA_glob) %>%
  filter(Month == example_month)



t <- df_proba_PCA_glob_sub %>%
  dplyr::select(-c(PCA1,PCA2,Month)) %>%
  group_by(x,y) %>%
  pivot_wider(names_from = Use,values_from = Proba) %>%
  group_by(x,y) %>%
  summarise(Rp_Pa = Rp_proba*Pa_proba,
            Rp_Co = Rp_proba*Co_proba,
            Rp_Vt = Rp_proba*Vt_proba,
            Rp_Ni = Rp_proba*Ni_proba,
            Rp_Lk = Rp_proba*Lk_proba)

t %>%
  pivot_longer(cols=starts_with("Rp_"),names_to="Pair",values_to="Proba") %>%
  #filter(co_occ > 0.2) %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = Proba)) +
  scale_fill_viridis() +
    facet_grid(~Pair)+
  coord_equal()+
  theme_bw()

# Calcule produit proba par paire d'usages
m <- data.frame(combn(names(label_test),2), stringsAsFactors = F)
df0 <- df_proba_PCA_glob_sub %>%
  dplyr::select(-c(PCA1,PCA2,Month)) %>%
  group_by(x,y) %>%
  pivot_wider(names_from = Use,values_from = Proba) %>%
  group_by(x,y) %>%
  dplyr::select(-c(x,y))
df <- df0[,-c(1:2)]
df1 <- data.frame(lapply(m, function(x) df[[x[2]]] * df[[x[1]]]))
colnames(df1) <- lapply(m, function(x) paste(x[2], x[1], sep="_"))
df3 <- cbind(df0[,1:2],df1)

# Mise en forme pour fact_grid
df4 <- df3 %>% 
  pivot_longer(cols=ends_with("proba"),names_to="pair",values_to="Proba") %>%
  mutate(u1 = substr(pair,1,2),
            u2 = substr(pair,10,11)) %>%
  dplyr::select(-pair) %>%
  as.data.frame()

df4$u1 <- factor(df4$u1, levels = c(substr(names(label_test),1,2)))
df4$u2 <- factor(df4$u2, levels = c(substr(names(label_test),1,2)))

df4 %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = Proba)) +
  scale_fill_viridis() +
  facet_grid(u2~u1, drop=F)+
  coord_equal()+
  theme_bw()

png(file = paste0(Gspace_path,"/probability_by_site_all_uses_",example_month,".png"),width=1400, height=800)
plot(PP)
dev.off()

# Distribution des probabilités, par usage
# bar
na.omit(df_proba_PCA_glob_sub) %>%
  ggplot(aes(x = Proba)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 0.01) +
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(Use~.,labeller = labeller(Use = label_test), scales="free")

na.omit(df_proba_PCA_glob_sub) %>%
  ggplot(aes(x = Proba)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 0.05) +
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(Use~.,labeller = labeller(Use = label_test), nrow=1)+
  labs(y="Frequency")

# lignes
df_proba_PCA_glob_sub %>%
  ggplot(aes(x=Proba, after_stat(density))) +
  geom_freqpoly()+
  facet_wrap(Use~.,labeller = labeller(Use = label_test))+
  xlim(0,1)

