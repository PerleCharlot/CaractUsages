### Titre -------------------------------------
# Nom : Modélisation des usages
# Auteure : Perle Charlot
# Date de création : 09-09-2022
# Dates de modification : 29-09-2022

### Librairies -------------------------------------

library(glmmfields)
library(MASS)
library(raster)
library(data.table)
library(sf)
# library(magick)
# library(corrplot)
library(ggplot2)
# library(dplyr)
# library(tidyverse)
# library(fasterize)

### Fonctions -------------------------------------

# Fonction qui sort les valeurs des var env associées à un usage donné
ExtractData1Use <- function(usage
                            ,fenetre_temporelle = liste.mois
){
  # # TEST
  # usage = "couchade"
  # fenetre_temporelle = liste.mois
  
  DataUsageMois <- function(mois){
    # # TEST
    # mois = fenetre_temporelle[2]
    
    raster.usage.mois.path = list.files(paste0(wd,"/output/par_periode/",mois),usage,full.names = T)
    if(length(raster.usage.mois.path)!= 0){
      raster.usage.mois = raster(raster.usage.mois.path)
      # RAsters environnement
      all.tif = list.files(paste0(gitCaractMilieu,"/output/ACP/"),pattern = ".tif",recursive = T,full.names = T)
      all.tif = all.tif[grep(mois,all.tif)]
      # Retirer dossiers : toutes, ACP_FAMD, CA/sans_ACP_clim
      all.tif = all.tif[!grepl("toutes",all.tif)]
      all.tif = all.tif[!grepl("ACP_FAMD",all.tif)]
      all.tif = all.tif[!grepl("CA/sans_ACP_clim",all.tif)]
      stack.env = stack(all.tif)
      
      # Combiner env + usage
      stack.env.usage = stack(raster.usage.mois, stack.env)
      # enlève les NA des bords
      stack.env.usage <- raster::crop(stack.env.usage, limiteN2000)
      # stack -> df
      df.env.usage = as.data.frame(as.data.table(stack.env.usage[]))
      df.env.usage = na.omit(df.env.usage)
      # Remove nom mois dans colonnes
      to.modif = names(df.env.usage)[grep(mois,names(df.env.usage))]
      new.names = unlist(lapply(to.modif, function(x) substring(x,1,nchar(x) - nchar(mois) - 1)))
      names(df.env.usage)[grep(mois,names(df.env.usage))] = new.names
      return(df.env.usage)
      
    } else{
      cat("\n L'usage",usage,"n'est pas présent au mois de",mois,".")
      return(NA)}
  }
  df.e.u = do.call(rbind,lapply(fenetre_temporelle, DataUsageMois))
  # Remove NA (from absence months)
  df.e.u = na.omit(df.e.u)
  return(df.e.u)
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
gitCaractMilieu = "C:/Users/perle.charlot/Documents/PhD/DATA/R_git/CaractMilieu"

#### Autre ####
liste.mois = c("mai","juin","juillet","aout","septembre")
# Liste dimensions
liste.dim =  c("CA","B","PV","CS","D","I")

### Programme -------------------------------------

# Présence/absence usage ~ variables du milieu
# Quelles variables du milieu prendre ? 
# A - toutes (~ 40)
# B - les axes des AFDM par dimension (~20)

# GLM : binomiale
# Option B = axes des AFDM par dimension

limiteN2000 <- st_read(limiteN2000)

df.Co = ExtractData1Use("couchade")
df.Ni = ExtractData1Use("nidification")


# Test de modèle pour le mois de juillet
# Sur l'usage Couchade (Co)
data_glm = df.Co

# Rapide visualisation sur quelques axes
ggplot(data_glm, aes(x=axe1_PV, y=axe2_PV, color=as.factor(couchade)))+
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

