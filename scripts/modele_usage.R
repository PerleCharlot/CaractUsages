### Titre -------------------------------------
# Nom : Modélisation des usages
# Auteure : Perle Charlot
# Date de création : 09-09-2022
# Dates de modification : 12-10-2022

### Librairies -------------------------------------
library(glmmfields)
library(MASS)
library(raster)
library(data.table)
library(sf)
library(caret)
library(DescTools)
library(pROC)
# library(magick)
# library(corrplot)
library(ggplot2)
# library(dplyr)
library(tidyverse)
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
      
      # ! ne sert à rien de le recalculer à chaque fois si existe déjà
      # Sauvegarder table des valeurs des rasters par mois
      if(!dir.exists(paste0(input_path,"/vars_env/",mois))){
        dir.create(paste0(input_path,"/vars_env/",mois))
      }
      if(length(list.files(paste0(input_path,"/vars_env/",mois))) == 0){
        df.env = as.data.frame(as.data.table(stack.env[]))
        df.env = cbind(coordinates(stack.env), df.env)
        # Remove nom mois dans colonnes
        to.modif = names(df.env)[grep(mois,names(df.env))]
        new.names = unlist(lapply(to.modif, function(x) substring(x,1,nchar(x) - nchar(mois) - 1)))
        names(df.env)[grep(mois,names(df.env))] = new.names
        write_csv2(df.env, paste0(input_path,"/vars_env/",mois,"/dt_vars_FAMD.csv"))
      } else{cat(paste0("\n Le tableau des variables FAMD du mois de ",mois," a déjà été calculé."))}
      
      # Combiner env + usage
      stack.env.usage = stack(raster.usage.mois, stack.env)
      # enlève les NA des bords
      stack.env.usage <- raster::crop(stack.env.usage, limiteN2000.shp)
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
  #return(df.e.u)
  write.csv(df.e.u, paste0(wd,"/output/niches/df_niche_",usage,".csv"))
}

# Crée un ENM (pour 1 usage donné), puis prédit tous les mois
CreateModelUsage <- function(nom_court_usage){
  
  # # TEST
  # nom_court_usage = "Ni"
  
  if(!dir.exists(paste0(output_path,"/niches/",nom_court_usage))){
    dir.create(paste0(output_path,"/niches/",nom_court_usage))
  }
  
  # Import données
  corresp_nom_us = data.frame(Nom=c("Nidification","Lek","Couchade","Pâturage","Randonnée","VTT"),
                              Nom_court = c("Ni","Lk","Co","Pa","Rp","Vt"),
                              nom_long = c("nidification","parade","couchade","paturage","randonnee_pedestre","VTT"))
  nom_beau = corresp_nom_us$Nom[corresp_nom_us$Nom_court == nom_court_usage]
  nom_lg = corresp_nom_us$nom_long[corresp_nom_us$Nom_court == nom_court_usage]
  data_glm = fread(paste0(output_path,"/niches/df_niche_",nom_lg,".csv"), drop="V1")
  names(data_glm)[1] = "usage"
  # Visualisation distribution sur axes FAMD
  plot = data_glm %>% 
    pivot_longer(!usage, names_to = "variables", values_to = "valeurs") %>%
    ggplot(aes(x=valeurs, y=variables, colour=as.factor(usage))) +
    geom_boxplot()+ 
    scale_colour_discrete(name = nom_beau, labels = c("Absence", "Présence"))+
    theme(text = element_text(size=14))
  # Sauvegarde
  png(file=paste0(output_path,"/niches/",nom_court_usage,"/viz_distrib_famd.png"), 
      width=1400, height=800)
  print(plot)
  dev.off()
  
  # Modèle
  sample <- sample(c(TRUE, FALSE), nrow(data_glm), replace=TRUE, prob=c(0.7,0.3))
  train <- data_glm[sample, ]
  test <- data_glm[!sample, ] 
  model.glm <- glm(usage ~ ., family=binomial, data=train)
  model.glm.step <- stepAIC(model.glm)
  # Enregistrer le modèle pour pouvoir ensuite echo = F
  save(model.glm.step, file = paste0(output_path,"/niches/",nom_court_usage,"/modele.rdata"))
  
  # Coefficients
  # exp(coef(glm.Co.step))
  
  # Tester fiabilité modèle : score de Brier, matrice de confusion et AUC/ROC
  brier_score = BrierScore(model.glm.step)
  # Prédire sur jeu test pour tester accuracy + AUC
  test$pred_prob <- predict(model.glm.step, test, type="response")
  
  test$pred_resp <- ifelse(test$pred_prob > 0.50, 1, 0)
  test$usage = as.factor(test$usage)
  test$pred_resp = as.factor(test$pred_resp)
  # Sauvegarde AUC/ROC
  png(file=paste0(output_path,"/niches/",nom_court_usage,"/roc.png"), 
      width=800, height=800)
  print(roc(test$usage ~ test$pred_prob, plot = TRUE, print.auc = TRUE))
  dev.off()
  # TODO : trouver meilleur seuil pour convertir proba en P/A
  # library(ROCR)
  # ROCR_pred_test <- prediction(test$pred_prob,test$usage)
  # ROCR_perf_test <- performance(ROCR_pred_test,'sens',"spec")
  # ROCR_perf_test2 <- performance(ROCR_pred_test,'tpr',"fpr")
  # plot(ROCR_perf_test2)
  # cost_perf = performance(ROCR_pred_test, "cost") 
  # ROCR_pred_test@cutoffs[[1]][which.min(cost_perf@y.values[[1]])]
  
  # Sauvegarde matrice de confusion
  CM = confusionMatrix(test$pred_resp, test$usage)
  save(CM, file=paste0(output_path,"/niches/",nom_court_usage,"/CM.rdata"))
  
  # Prédiction spatialisée
  predUsageMois <- function(mois){
    # # TEST
    # mois = "mai"
    
    if(!dir.exists(paste0(output_path,"/niches/",nom_court_usage,"/predictions/"))){
      dir.create(paste0(output_path,"/niches/",nom_court_usage,"/predictions/"))
    }
    
    t = try(raster(paste0(output_path,"/par_periode/",mois,"/",nom_lg,".tif")),
            silent = TRUE)
    if(inherits(t, "try-error")) {
      # Si l'usage n'est pas présent le mois considéré
      cat(paste0("\nL'usage ", nom_beau," n'est pas présent le mois ",mois,"."))
    } else{
      
      df.env <- fread(paste0(input_path,"/vars_env/",mois,"/dt_vars_FAMD.csv"),dec=",")
      df.env$prob <- predict(model.glm.step, df.env, type="response")
      df.env$pred <- ifelse(df.env$prob > 0.50, 1, 0)
      
      raster_prob = rasterFromXYZ(data.frame(df.env$x, df.env$y, df.env$prob), crs=EPSG_2154)
      raster_obs = raster(paste0(output_path,"/par_periode/",mois,"/",nom_lg,".tif"))
      # mask (car 0 confondu avec NA hors N2000)
      raster_obs <- mask(raster_obs, limiteN2000.shp)
      
      all_rasters = stack(raster_obs, raster_prob)
      names(all_rasters) = c("observation","probabilite_prediction")
      # plot(all_rasters, colNA='black')
      
      writeRaster(all_rasters, 
                  paste0(output_path,"/niches/",nom_court_usage,
                         "/predictions/rasters_pred_",mois,".tif"),
                  overwrite=TRUE)
      cat(paste0("\nL'usage ", nom_beau," a été prédit pour le mois de ",mois,"."))
    }
  }
  
  lapply(liste.mois, predUsageMois)
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
# Liste usages
liste.usages = c("Ni","Lk","Co","Pa","Rp","Vt")

### Programme -------------------------------------

# Présence/absence usage ~ variables du milieu
# Quelles variables du milieu prendre ? 
# A - toutes (~ 40)
# B - les axes des AFDM par dimension (~20)

#### GLM : binomiale ####
# Option B = axes des AFDM par dimension

limiteN2000.shp <- st_read(limiteN2000)

# Création d'un df par usage : utilisation fonction ExtractData1Use
lapply(c("nidification",
         "couchade",
         "paturage",
         "randonnee_pedestre",
         "VTT",
         "parade"),ExtractData1Use)

# Exploitation des données : création modèle
set.seed(1)
lapply(liste.usages, CreateModelUsage)

# raster_pred = rasterFromXYZ(data.frame(df.env$x, df.env$y, df.env$pred), crs=EPSG_2154)
# # - raster accord entre obs et pred
# somme_rast = raster_pred + raster_obs
# raster_accord <- somme_rast ==2
# # - raster prédit mais non obs
# diff_rast = raster_pred - raster_obs
# raster_predNobs <- diff_rast == 1
# # - raster obs mais non prédit
# diff_rast2 = raster_obs - raster_pred
# raster_obsNpred <- diff_rast2 == 1

# library(mapview)
# mapviewOptions(na.color = "transparent",
#                basemaps="OpenStreetMap",
#                viewer.suppress=TRUE,
#                trim=TRUE)
# view(raster_obs)

# à mettre dans le rmd
load(file = paste0(output_path,"/niches/",nom_court_usage,"/modele.rdata"))
summary(model.glm.step) # display results
load(file=paste0(output_path,"/niches/",nom_court_usage,"/CM.rdata"))
CM


# # Test de modèle sur l'usage Couchade (Co)
# data_glm = df.Ni
# names(data_glm)[1] = "usage"

# #make this example reproducible
# set.seed(1)
# #Use 70% of dataset as training set and remaining 30% as testing set
# sample <- sample(c(TRUE, FALSE), nrow(data_glm), replace=TRUE, prob=c(0.7,0.3))
# train <- data_glm[sample, ]
# test <- data_glm[!sample, ] 
# # Fit modèle
# model.glm <- glm(usage ~ .,
#     family=binomial, 
#     data=train)
# summary(model.glm) # display results
# # Sélection des variables
# model.glm.step <- stepAIC(model.glm)
# # glm.Co.step$anova
# # summary(glm.Co.step)
# #confint(glm.Co.step) # 95% CI for the coefficients
# # Odds ratio
# # exp(coef(glm.Co.step)) # exponentiated coefficients
# # exp(confint(glm.Co)) # 95% CI for exponentiated coefficients
# # pred = predict(glm.Co.step, type="response") # predicted values
# # residuals(glm.Co.step, type="deviance") # residuals 
# 
# # #null.model <- glm(Co ~ 1, family = binomial,data=data_glm)
# # pseudoR2 <- (glm.Co.step$null.deviance - glm.Co.step$deviance)/glm.Co.step$null.deviance
# # pseudoR2
# 
# ## S3 method for class 'glm'
# library(DescTools)
# BrierScore(model.glm.step)
# 
# # Etape AUC/ROC tout pour trouver le seuil de proba où ça bascule en présence
# #https://stats.stackexchange.com/questions/172945/rmse-root-mean-squared-error-for-logistic-models
# 
# library(pROC)
# test$pred_prob <- predict(model.glm.step, test, type="response")
# test_roc = roc(test$usage ~ test$pred_prob, plot = TRUE, print.auc = TRUE)
# 
# # exp(coef(glm.Co.step))
# 
# test$pred_resp <- ifelse(test$pred_prob > 0.50, 1, 0)
# 
# test$usage = as.factor(test$usage)
# test$pred_resp = as.factor(test$pred_resp)
# 
# confusionMatrix(test$pred_resp, test$usage)


#### GLM spatial #####
# #https://cran.r-project.org/web/packages/glmmfields/vignettes/spatial-glms.html
# m_spatial <- glmmfields(Co ~ axe1_CA + axe2_CA + axe3_CA +
#                           axe1_B + axe2_B + axe3_B +
#                           axe1_PV + axe2_PV + axe3_PV +
#                           axe1_CS + axe2_CS + axe3_CS + 
#                           axe1_D + axe2_D + 
#                           axe1_I + axe2_I + axe3_I,
#                         data = data_glm, 
#                         family = binomial(link = "logit"),
#                         lat = "y", lon = "x", 
#                         nknots = 5, iter = 500, chains = 4,
#                         seed = 123 # passed to rstan::sampling()
#                         )
# m_spatial
# plot(m_spatial, type = "spatial-residual", link = TRUE) +
#   geom_point(size = 3)
# # Residual vs fitted
# plot(m_spatial, type = "residual-vs-fitted")
# # link scale:
# p <- predict(m_spatial)
# head(p)
# # response scale:
# p <- predict(m_spatial, type = "response")
# head(p)
# # get our parameter estimates as a data frame:
# head(tidy(m_spatial, conf.int = TRUE, conf.method = "HPDinterval"))
# 



#### RF ####

