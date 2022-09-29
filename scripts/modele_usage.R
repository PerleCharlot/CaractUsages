### Titre -------------------------------------
# Nom : Modélisation des usages
# Auteure : Perle Charlot
# Date de création : 09-09-2022
# Dates de modification : 12-09-2022

### Librairies -------------------------------------

library(glmmfields)
library(MASS)

# library(raster)
library(data.table)
# library(sf)
# library(magick)
# library(corrplot)
library(ggplot2)
# library(dplyr)
# library(tidyverse)
# library(fasterize)

### Fonctions -------------------------------------

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

#### Tables ####
gitCaractMilieu = "C:/Users/perle.charlot/Documents/PhD/DATA/R_git/CaractMilieu"


#### Autre ####
liste.mois = c("mai","juin","juillet","aout","septembre")
# Liste dimensions
liste.dim=  c("CA","B","PV","CS","D","I")

### Programme -------------------------------------

# Présence/absence usage ~ variables du milieu
# Quelles variables du milieu prendre ? 
# A - toutes (~ 40)
# B - les axes des AFDM par dimension (~20)

# GLM : binomiale
# Option B = axes des AFDM par dimension

MergeTbl_for1Month <- function(NOM_MOIS1){
  # #TEST
  # NOM_MOIS1 = liste.mois[1]
  
  MergeTbl_for1Month_accrossDim <- function(lettre_dimension,
                                          nom_mois=NOM_MOIS1,
                                          ACPclimat = "oui"){
    # #TEST
    # nom_mois = "juillet"
    # lettre_dimension = "CA"
    # ACPclimat = "oui"
    
    corresp = data.frame(numero_saison=c("05","06","07","08","09"),
                         periode = c("mai",'juin','juillet','aout','septembre'))
    
    num_mois = corresp$numero_saison[which(corresp$periode == nom_mois)]
    
    # Conditions spécifiques pour dim CA
    if(lettre_dimension=="CA"){
      if(ACPclimat == "oui"){
        path_table = paste0(gitCaractMilieu,"/output/ACP/",lettre_dimension,"/avec_ACP_clim/",
                            num_mois,nom_mois,"/axesFAMD_",lettre_dimension,"_usages_",nom_mois,".csv")
      }else{
        path_table = paste0(gitCaractMilieu,"/output/ACP/",lettre_dimension,"/sans_ACP_clim/",
                            num_mois,nom_mois,"/axesFAMD_",lettre_dimension,"_usages_",nom_mois,".csv")
      }
    } else {
      path_table = paste0(gitCaractMilieu,"/output/ACP/",lettre_dimension,"/",
                          num_mois,nom_mois,"/axesFAMD_",lettre_dimension,"_usages_",nom_mois,".csv")
    }
    
    table_vars_usa <- fread(path_table, drop="V1")
    return(table_vars_usa)
  }
  
  # Obtenir une table avec en colonnes 
  liste_Tbl_MOIS = lapply(liste.dim,  MergeTbl_for1Month_accrossDim)
  # Identifier les usages présents (pour remplir le by de merge)
  index_axe3 = grep("axe3",names(liste_Tbl_MOIS[[1]])) +1
  usages_presents = names(liste_Tbl_MOIS[[1]])[index_axe3:length(names(liste_Tbl_MOIS[[1]]))]
  # appliquer merge à travers la liste de 6 tables
  Tbl_MOIS = as.data.frame(Reduce(function(x, y) merge(x, y, by=c("x","y",usages_presents)), liste_Tbl_MOIS))
  
  return(Tbl_MOIS)
  
}
# Sort liste de 5 tables, 1/mois
data_vars_usages <- lapply(liste.mois, MergeTbl_for1Month)

# Test de modèle pour le mois de juillet
# Sur l'usage Couchade (Co)
data_glm = data_vars_usages[[3]]

# Rapide visualisation sur quelques axes
ggplot(data_glm, aes(x=axe1_PV, y=axe2_PV, color=as.factor(Co)))+
  geom_point(alpha=0.2)
ggplot(data_glm, aes(x=axe1_B, y=axe2_B, color=as.factor(Co)))+
  geom_point(alpha=0.2)

# Fit modèle
glm.Co <- glm(Co ~ axe1_CA + axe2_CA + axe3_CA +
      axe1_B + axe2_B + axe3_B +
      axe1_PV + axe2_PV + axe3_PV +
      axe1_CS + axe2_CS + axe3_CS + 
      axe1_D + axe2_D + 
      axe1_I + axe2_I + axe3_I,
    family=binomial, 
    data=data_glm)
summary(glm.Co) # display results

# Sélection des variables
glm.Co.step <- stepAIC(glm.Co)
glm.Co.step$anova
# best model sans axe1 et axe3 de I, et sans axe1 de D


summary(glm.Co.step)

confint(glm.Co.step) # 95% CI for the coefficients

# Odds ratio
exp(coef(glm.Co.step)) # exponentiated coefficients
exp(confint(glm.Co)) # 95% CI for exponentiated coefficients
pred = predict(glm.Co.step, type="response") # predicted values
residuals(glm.Co.step, type="deviance") # residuals 

#null.model <- glm(Co ~ 1, family = binomial,data=data_glm)
pseudoR2 <- (glm.Co.step$null.deviance - glm.Co.step$deviance)/glm.Co.step$null.deviance
pseudoR2


# Etape AUC/ROC tout pour trouver le seuil de proba où ça bascule en présence
#https://stats.stackexchange.com/questions/172945/rmse-root-mean-squared-error-for-logistic-models

#### GLM spatial #####
#https://cran.r-project.org/web/packages/glmmfields/vignettes/spatial-glms.html
m_spatial <- glmmfields(Co ~ axe1_CA + axe2_CA + axe3_CA +
                          axe1_B + axe2_B + axe3_B +
                          axe1_PV + axe2_PV + axe3_PV +
                          axe1_CS + axe2_CS + axe3_CS + 
                          axe1_D + axe2_D + 
                          axe1_I + axe2_I + axe3_I,
                        data = data_glm, 
                        family = binomial(link = "logit"),
                        lat = "y", lon = "x", 
                        nknots = 5, iter = 500, chains = 4,
                        seed = 123 # passed to rstan::sampling()
                        )
m_spatial
plot(m_spatial, type = "spatial-residual", link = TRUE) +
  geom_point(size = 3)
# Residual vs fitted
plot(m_spatial, type = "residual-vs-fitted")
# link scale:
p <- predict(m_spatial)
head(p)
# response scale:
p <- predict(m_spatial, type = "response")
head(p)
# get our parameter estimates as a data frame:
head(tidy(m_spatial, conf.int = TRUE, conf.method = "HPDinterval"))






# RF

