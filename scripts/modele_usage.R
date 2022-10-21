### Titre -------------------------------------
# Nom : Modélisation des usages
# Auteure : Perle Charlot
# Date de création : 09-09-2022
# Dates de modification : 21-10-2022

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
library(FactoMineR)
library(factoextra)
library(ggtext)
library(glue)
### Fonctions -------------------------------------

# Fonction qui vérifie que le raster ait le bon CRS et extent, et le modifie si besoin
AjustExtCRS <- function(path.raster.to.check, path.raster.ref=chemin_mnt){
  
  # # TEST
  # path.raster.to.check = r_hiv[2]
  # raster.ref = MNT
  
  raster.to.check <- raster(path.raster.to.check)
  raster.ref <- raster(path.raster.ref)
  
  ext.to.check <- extent(raster.to.check)
  bon.extent <- extent( raster.ref)
  
  sameCRS <- compareCRS(raster.to.check,EPSG_2154)
  sameExtent <- (ext.to.check == bon.extent)
  
  if(any(!sameCRS,!sameExtent)) {
    raster.to.check <- projectRaster(raster.to.check, raster.ref)
    writeRaster(raster.to.check, path.raster.to.check, overwrite=TRUE)
    cat("\nRaster ", names(raster.to.check)," a été modifié et sauvegarde.")
  }
  
}

# Fonction qui sort les valeurs des var env associées à un usage donné
ExtractData1Use <- function(usage,
                            fenetre_temporelle = liste.mois,
                            type_donnees # "brute", "ACP"
){
  # # TEST
  # usage = "couchade"
  # fenetre_temporelle = liste.mois
  # type_donnees = "brute"
  
  DataUsageMois <- function(mois){
    # # TEST
    # mois = fenetre_temporelle[2]
    
    # Raster de l'usage considéré
    raster.usage.mois.path = list.files(paste0(wd,"/output/par_periode/",mois),
                                        usage,full.names = T)
    if(length(raster.usage.mois.path)!= 0){
      raster.usage.mois = terra::rast(raster.usage.mois.path)
      
      # Rasters environnement : sur le mois considéré
      if(type_donnees == "ACP"){
        all.tif = list.files(paste0(gitCaractMilieu,"/output/ACP/"), 
                             pattern = ".tif",recursive = T,full.names = T)
        all.tif = all.tif[grep(mois,all.tif)]
        # Retirer dossiers : toutes, ACP_FAMD, CA/sans_ACP_clim
        all.tif = all.tif[!grepl("toutes",all.tif)]
        all.tif = all.tif[!grepl("ACP_FAMD",all.tif)]
        all.tif = all.tif[!grepl("CA/sans_ACP_clim",all.tif)]
        # Sauvegarder table des valeurs des rasters par mois
        if(!dir.exists(paste0(input_path,"/vars_env/ACP/",mois))){
          dir.create(paste0(input_path,"/vars_env/ACP/",mois),recursive=T)
        }
        DOS = paste0(input_path,"/vars_env/ACP/",mois)
      } 
      if(type_donnees == "brute"){
        all.tif = list.files(paste0(gitCaractMilieu,"/output/stack_dim/"),
                             pattern = ".tif",recursive = T,full.names = T)
        all.tif = all.tif[grep(mois,all.tif)]
        # Sauvegarder table des valeurs des rasters par mois
        if(!dir.exists(paste0(input_path,"/vars_env/brute/",mois))){
          dir.create(paste0(input_path,"/vars_env/brute/",mois),recursive=T)
        }
        DOS = paste0(input_path,"/vars_env/brute/",mois)
      }
      # Stacker les variables de l'environnement
      stack.env2 = stack(all.tif)
      stack.env = terra::rast(all.tif)
      
      # Correction noms des vars en brute
      findCorrectName <- function(y){
        index_noms_faux = grep(pattern = y , x = noms_faux)
        if(length(index_noms_faux)!=0){
          return(data.frame(old_name =noms_faux[index_noms_faux],new_name=y))
        } else{return(c(old_name =NA,new_name=NA))}
      }
      if(type_donnees == "brute"){
        noms_faux = names(stack.env)[!names(stack.env) %in% table_variables$Nom]
        table_correctName = do.call(rbind,lapply(table_variables$Nom, findCorrectName))
        table_correctName = na.omit(table_correctName)
        # matcher les noms
        table_correctName = table_correctName[order(match(table_correctName$old_name, noms_faux)),]
        names(stack.env)[!names(stack.env) %in% table_variables$Nom] <- table_correctName$new_name
      }
      
      FILS = list.files(DOS)
      
      # ! ne sert à rien de le recalculer à chaque fois si existe déjà
      if(length(FILS) == 0){
        df.env = as.data.frame(as.data.table(stack.env[]))
        # Rendre noms variables propres
        if(type_donnees == "ACP"){
          # Remove nom mois dans colonnes
          to.modif = names(df.env)[grep(mois,names(df.env))]
          new.names = unlist(lapply(to.modif, function(x) substring(x,1,nchar(x) - nchar(mois) - 1)))
          names(df.env)[grep(mois,names(df.env))] = new.names
        }
        df.env = cbind(coordinates(stack.env2), df.env)
        write_csv2(df.env, paste0(DOS,"/dt_vars.csv"))
      } else{cat(paste0("\n Le tableau des variables ",
                        type_donnees," du mois de ",mois," a déjà été calculé."))}
      
      # Combiner env + usage
      stack.env.usage = c(raster.usage.mois, stack.env)
      # # Combiner env + usage
      # stack.env.usage = stack(raster.usage.mois, stack.env)
      # enlève les NA des bords
      stack.env.usage <- raster::crop(stack.env.usage, limiteN2000.shp)
      # stack -> df
      df.env.usage = as.data.frame(as.data.table(stack.env.usage[]))
      df.env.usage = na.omit(df.env.usage)
     
      # Arranger les noms, again
      if(type_donnees == "ACP"){
        # Remove nom mois dans colonnes
        to.modif = names(df.env.usage)[grep(mois,names(df.env.usage))]
        new.names = unlist(lapply(to.modif, function(x) substring(x,1,nchar(x) - nchar(mois) - 1)))
        names(df.env.usage)[grep(mois,names(df.env.usage))] = new.names
      }
      # TODO: ajouter les coordonnées x y ???
      return(df.env.usage)
    } else{
      cat("\n L'usage",usage,"n'est pas présent au mois de",mois,".")
      return(NA)}
  }
  
  df.e.u = do.call(rbind,lapply(fenetre_temporelle, DataUsageMois))
  # Remove NA (from absence months)
  df.e.u = na.omit(df.e.u)
  
  if(!dir.exists(paste0(wd,"/output/niches/",type_donnees,"/"))){
    dir.create(paste0(wd,"/output/niches/",type_donnees,"/"),recursive=T)
  }
  write.csv(df.e.u, paste0(wd,"/output/niches/",type_donnees,"/df_niche_",usage,".csv"))
}

# Crée un ENM (pour 1 usage donné), puis prédit tous les mois
CreateModelUsage <- function(nom_court_usage, type_donnees){
  
  # # TEST
  # nom_court_usage = "Ni"
  # type_donnees = "brute" # "ACP"
  
  
  if(!dir.exists(paste0(output_path,"/niches/",type_donnees,"/",nom_court_usage))){
    dir.create(paste0(output_path,"/niches/",type_donnees,"/",nom_court_usage))
  }
  
  # Import données
  corresp_nom_us = data.frame(Nom=c("Nidification","Lek","Couchade","Pâturage","Randonnée","VTT"),
                              Nom_court = c("Ni","Lk","Co","Pa","Rp","Vt"),
                              nom_long = c("nidification","parade","couchade","paturage","randonnee_pedestre","VTT"))
  nom_beau = corresp_nom_us$Nom[corresp_nom_us$Nom_court == nom_court_usage]
  nom_lg = corresp_nom_us$nom_long[corresp_nom_us$Nom_court == nom_court_usage]
  data_glm = fread(paste0(output_path,"/niches/",type_donnees,"/df_niche_",nom_lg,".csv"), drop="V1")
  names(data_glm)[1] = "usage"
  # Visualisation distribution le long des variables environnementales
  if(type_donnees == "brute"){

    
    # Nécessité de scale les données, sinon les valeurs extrêmes rendent le plot illisible
    data_glm_scaled = as.data.frame(scale(data_glm[,2:dim(data_glm)[2]]))
    data_glm_scaled = cbind(usage = data_glm$usage, data_glm_scaled)

    # Colorer les variables par couleur de dimension
    data_glm_s_p = data_glm_scaled %>% 
      pivot_longer(!usage, names_to = "variables", values_to = "valeurs")

    col_vars_sub = col_vars %>% subset(select=c(colour_dim,Nom,dim))
    
    #data_glm_s_p2 = merge(data_glm_s_p, col_vars_sub, by.x="variables", by.y="Nom")
    
    # labs = col_vars$Nom[col_vars$Nom %in% names(data_glm)]
    # labs = na.omit(labs[match(names(data_glm_scaled),labs)])
    # color = col_vars$colour_dim[col_vars$Nom %in% labs]
    # name <- glue("<i style='color:{color}'>{labs}")
    # # TESt initial
    # plot =  data_glm_scaled %>% 
    #   pivot_longer(!usage, names_to = "variables", values_to = "valeurs") %>%
    #   ggplot(aes(x=valeurs, y=variables, colour=as.factor(usage))) +
    #   geom_boxplot()+ 
    #   scale_y_discrete(labels = name)+
    #   scale_colour_discrete(name = nom_beau, labels = c("Absence", "Présence"))+
    #   theme(text = element_text(size=14),axis.text.y = element_markdown())

    a = sort(unique(data_glm_s_p$variables), decreasing=TRUE)
    a = data.frame(Nom=a)
    b = merge(a, col_vars_sub, by="Nom",sort=FALSE)
    b2 = b$colour_dim
    
    plot = data_glm_s_p %>%
      ggplot(aes(x=valeurs, y=variables, colour=as.factor(usage))) +
      geom_boxplot()+ 
      scale_colour_discrete(name = nom_beau, labels = c("Absence", "Présence"))+
      theme(text = element_text(size=14), axis.text.y = element_text(colour = rev(b2)))
  }
  if(type_donnees == "ACP"){
    plot =  data_glm %>% 
      pivot_longer(!usage, names_to = "variables", values_to = "valeurs") %>%
      ggplot(aes(x=valeurs, y=variables, colour=as.factor(usage))) +
      geom_boxplot()+ 
      scale_colour_discrete(name = nom_beau, labels = c("Absence", "Présence"))+
      theme(text = element_text(size=14))
  }

  # Sauvegarde
  png(file=paste0(output_path,"/niches/",type_donnees,"/",nom_court_usage,"/viz_distrib_famd.png"), 
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
  save(model.glm.step, file = paste0(output_path,"/niches/",type_donnees,
                                     "/",nom_court_usage,"/modele.rdata"))
  
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
  png(file=paste0(output_path,"/niches/",type_donnees,
                  "/",nom_court_usage,"/roc.png"), 
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
  save(CM, file=paste0(output_path,"/niches/",type_donnees,
                       "/",nom_court_usage,"/CM.rdata"))
  
  # Prédiction spatialisée
  predUsageMois <- function(mois){
    # # TEST
    # mois = "juin"
    
    if(!dir.exists(paste0(output_path,"/niches/",type_donnees,"/",
                          nom_court_usage,"/predictions/"))){
      dir.create(paste0(output_path,"/niches/",type_donnees,"/",
                        nom_court_usage,"/predictions/"))
    }
    
    t = try(raster(paste0(output_path,"/par_periode/",mois,"/",nom_lg,".tif")),
            silent = TRUE)
    if(inherits(t, "try-error")) {
      # Si l'usage n'est pas présent le mois considéré
      cat(paste0("\nL'usage ", nom_beau," n'est pas présent le mois ",mois,"."))
    } else{
      
      df.env <- fread(paste0(input_path,"/vars_env/",type_donnees,"/",
                             mois,"/dt_vars.csv"),dec=",")
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
                  paste0(output_path,"/niches/",type_donnees,"/",
                         nom_court_usage,
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
path_table_variables <- paste0(gitCaractMilieu,"/input/liste_variables.csv")

#### Autre ####
liste.mois = c("mai","juin","juillet","aout","septembre")
df.mois = data.frame(nom_mois = liste.mois, numero_mois = c("05","06","07","08","09"))
# Liste dimensions
liste.dim =  c("CA","B","PV","CS","D","I")
# Liste usages
liste.usages = c("Ni","Lk","Co","Pa","Rp","Vt")

### Programme -------------------------------------

# Lignes à valider
limiteN2000.shp <- st_read(limiteN2000)
# Colorer les variables en fonction des dimensions
corresp_col = data.frame(dim = liste.dim,
                         colour_dim = c("dodgerblue","darkgoldenrod1","darkgreen",
                                        "brown","blueviolet","darkgray"))
dim_vars = data.frame(vars_noms= apply(expand.grid(paste0("axe",1:3,"_"), liste.dim), 1, paste, collapse=""),
                      dim = unname(unlist(as.list(data.frame(t(replicate(3,liste.dim)))))))
dim_col = merge(dim_vars, corresp_col, by='dim',all=T)


table_variables <- fread(path_table_variables, header=T)

col_vars <- merge(corresp_col, table_variables, by.x="dim",by.y="Dimension")

#### GLM : binomiale ####
# axes des AFDM par dimension

# Création d'un df par usage : utilisation fonction ExtractData1Use
lapply(c("nidification",
         "couchade",
         "paturage",
         "randonnee_pedestre",
         "VTT",
         "parade"),function(x) ExtractData1Use(usage=x, type_donnees = "brute"))

# TODO : enlever ACP1_clim et ACP2_clim qui ne doivent pas être dans les vars

# Exploitation des données : création modèle
set.seed(1)
lapply(liste.usages, function(x) CreateModelUsage(nom_court_usage=x,type_donnees = "brute"))

# TODO : modif fct CreateModelUsage pour considérer vars brute

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

#### Visualisation espace ACP ####


VizIndDim <- function(usage){
  # TEST
  mois = "juillet"
  usage = "Ni"
  
  for(mois in liste.mois){
    
    data_mois <- fread(paste0(gitCaractMilieu,"/output/ACP/ACP_FAMD/",
                              df.mois$numero_mois[which(df.mois$nom_mois == mois)],mois,
                              "/tblFAMD_ACP_FAMD_",mois,".csv"),drop="V1")
    tbl_data = data_mois[,1:17]
    
    pca1 <- PCA(tbl_data, graph = FALSE)
    #PCA -> ggplot
    # https://tem11010.github.io/Plotting-PCAs/
    tbl_data$pc1 <- pca1$ind$coord[, 1]
    tbl_data$pc2 <- pca1$ind$coord[, 2]  
    pca.vars <- pca1$var$coord %>% data.frame
    pca.vars$vars <- rownames(pca.vars)
    
    #https://stackoverflow.com/questions/51219267/pca-scaling-not-applied-to-individuals
    # #https://rdrr.io/cran/factoextra/src/R/fviz_pca.R#sym-fviz_pca_biplot
    # r_scale = (max(pca1$ind$coord) - min(pca1$ind$coord))/
    #   (max(pca1$var$coord) - min(pca1$var$coord))
    
    N = nchar(mois) + 1
    pca.vars$vars_noms = substr(pca.vars$vars,1, nchar(pca.vars$vars)-N)
    
    pca.vars.m <- melt(pca.vars, id.vars = "vars",warning =F)
    circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
      r = diameter / 2
      tt <- seq(0,2*pi,length.out = npoints)
      xx <- center[1] + r * cos(tt)
      yy <- center[2] + r * sin(tt)
      return(data.frame(x = xx, y = yy))
    }
    circ <- circleFun(c(0,0),2,npoints = 500)
    
    seuil = 0.2
    pca.vars.2 = pca.vars %>%
      filter(Dim.1 > seuil | Dim.1 < -seuil | Dim.2 > seuil | Dim.2 < -seuil) 
    nb_var_drop = dim(pca.vars)[1] -dim(pca.vars.2)[1]
    cat(paste0(nb_var_drop," variables ont été retirées pour faciliter la visualisation."))
    
    # Valeurs des usages
    # USAGE
    liste_rast = list.files(paste0(output_path,"/niches/",usage,"/predictions/"),".tif",
                            full.names = T)
    if(length(liste_rast[grep(mois, liste_rast)])==0){
      cat("\nPas d'usage ",usage,"au mois de ",mois,".")
    }else{
    r = stack(liste_rast[grep(mois, liste_rast)])
    names(r)=c("observation", "probabilite_prediction")
    obs_sp = rasterToPolygons(r$observation,fun=function(x)x==1, dissolve=TRUE)
    # AXES ACP
    r.AFDM <- stack(list.files(paste0(gitCaractMilieu,"/output/ACP/ACP_FAMD/",
                                      df.mois$numero_mois[which(df.mois$nom_mois == mois)],mois,
                                      "/"),".tif", full.names = T))
    # EXTRACTION VALEURS
    library(exactextractr)
    values_usages <- as.data.frame(exact_extract(r.AFDM, obs_sp, include_xy=T))
    head(values_usages)
    values_usages$axe1_ACP_FAMD = values_usages[,1]
    values_usages$axe2_ACP_FAMD = values_usages[,2]
    
    pca.vars.2 = merge(pca.vars.2, dim_col, by="vars_noms")
    # Densité individus + flèches variables, en couleur de dimension
    P3 = ggplot() +
      geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
      geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_hex(data=values_usages, aes(x="axe1_ACP_FAMD", 
                                       y="axe2_ACP_FAMD")) +
      geom_segment(data = pca.vars.2, aes(x = 0, xend = Dim.1*4, y = 0, yend = Dim.2*4),
                   arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
                   lwd = 0.6, colour='black') + 
      geom_text(data = pca.vars.2,size = 10,
                aes(x = Dim.1*3.5, y =  Dim.2*3.5,
                    label = vars_noms),colour=pca.vars.2$colour_dim,
                check_overlap = F) +
      labs(x="ACP 1",y="ACP2", fill="Densité\n(en pixels)",
           title= paste0("Usage ",usage," - ",mois))+
      coord_equal() +
      scale_fill_continuous(type = "viridis",direction=-1)+
      theme_bw() +
      theme(panel.grid = element_blank(), text = element_text(size=15),
            panel.border = element_rect(fill= "transparent"))
    
    png(file=paste0(output_path,"/niches/",usage,"/densite_dim_",mois,".png"), 
        width=1400, height=800)
    print(P3)
    dev.off()
    
    # Graph du cercle de corrélation avec toutes les variables
    pca.vars
    pca.vars = merge(pca.vars, dim_col, by="vars_noms")
    P1 = ggplot() +
      geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
      geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_segment(data = pca.vars, aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),
                   arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
                   lwd = 0.5) + 
      geom_text(data = pca.vars,size=10,
                aes(x = Dim.1*1.15, y =  Dim.2*1.15,
                    label = vars_noms),colour=pca.vars$colour_dim,
                check_overlap = F) +
      labs(x="ACP 1",y="ACP2", title= mois)+
      coord_equal() +
      theme_minimal() +
      theme(panel.grid = element_blank(), 
            text = element_text(size=15),
            panel.border = element_rect(fill= "transparent"))
    
    png(file=paste0(gitCaractMilieu,"/output/ACP/ACP_FAMD/",
                    df.mois$numero_mois[which(df.mois$nom_mois == mois)],mois,
                    "/plot_rmd/cercle_correlation_",mois,".png"), 
        width=1400, height=800)
    print(P1)
    dev.off()
    
    # Graph du cercle de corrélation avec les variables (> seuil)
    P2 = ggplot() +
      geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
      geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_segment(data = pca.vars.2, aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),
                   arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
                   lwd = 0.5) + 
      geom_text(data = pca.vars.2,size=10,
                aes(x = Dim.1*1.15, y =  Dim.2*1.15,
                    label = vars_noms),colour=pca.vars.2$colour_dim,
                check_overlap = F) +
      labs(x="ACP 1",y="ACP2", title= mois)+
      coord_equal() +
      theme_minimal() +
      theme(panel.grid = element_blank(), text = element_text(size=15),
            panel.border = element_rect(fill= "transparent"))
    
    png(file=paste0(gitCaractMilieu,"/output/ACP/ACP_FAMD/",
                    df.mois$numero_mois[which(df.mois$nom_mois == mois)],mois,
                    "/plot_rmd/cercle_correlation_seuil_",mois,".png"), 
        width=1400, height=800)
    print(P2)
    dev.off()}
  }
}

VizIndDim("Ni")





