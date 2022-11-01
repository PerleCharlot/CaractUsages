### Titre -------------------------------------
# Nom : Modélisation des usages
# Auteure : Perle Charlot
# Date de création : 09-09-2022
# Dates de modification : 01-11-2022

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
library(exactextractr)
### Fonctions -------------------------------------

# Fonction qui encode des colonnes d'un df en factoriel, si leur nature
# est notée "qualitative" dans une autre tbale de référence
RecodeFacto <- function(tableau_donnees_a_recoder, tableau_nature_donnees){
  # # TEST
  # tableau_donnees_a_recoder = tbl_data
  # tableau_nature_donnees = col_vars_formate # Avec 2 colonnes "ID" et "nature"
  
  tableau_donnees_a_recoder = as.data.frame(tableau_donnees_a_recoder)
  
  vars_quali_all = tableau_nature_donnees$ID[tableau_nature_donnees$nature == "qualitative"] 
  ind_vars_quali_tbl = which(names(tableau_donnees_a_recoder) %in% vars_quali_all)
  if(length(ind_vars_quali_tbl) == 0){
    return(tableau_donnees_a_recoder)
  }else{
    ind_vars_quanti_tbl = which(!names(tableau_donnees_a_recoder) %in% vars_quali_all)
    
    vars_quali_fact = as.data.frame(lapply(ind_vars_quali_tbl, 
                                           function(x) as.factor(tableau_donnees_a_recoder[,x])))
    names(vars_quali_fact) = names(tableau_donnees_a_recoder)[ind_vars_quali_tbl]
    
    tbl_quanti = data.frame(tableau_donnees_a_recoder[,ind_vars_quanti_tbl])
    names(tbl_quanti) = names(tableau_donnees_a_recoder)[ind_vars_quanti_tbl]
    
    tableau_donnees_recoded = cbind(vars_quali_fact, 
                                    tbl_quanti)

    return(tableau_donnees_recoded)
  }
}

# Fonction pour tracer cercle corrélation
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

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
                            type_donnees # "brute", "ACP", "axes_AFDM"
){
  # # TEST
  # usage = "paturage"
  # fenetre_temporelle = liste.mois
  # type_donnees = "ACP"
  
  DataUsageMois <- function(mois){
    # # TEST
    # mois = fenetre_temporelle[2]
    
    # Raster de l'usage considéré
    raster.usage.mois.path = list.files(paste0(wd,"/output/par_periode/",mois),
                                        pattern= paste0(usage,".tif$"),
                                        full.names = T)
    if(length(raster.usage.mois.path)!= 0){
      raster.usage.mois = terra::rast(raster.usage.mois.path)
      
      # Rasters environnement : sur le mois considéré
      if(type_donnees == "ACP"){
        all.tif = list.files(paste0(gitCaractMilieu,"/output/ACP/toutes/"), 
                             pattern = ".tif",recursive = T,full.names = T)
        all.tif = all.tif[grep(mois,all.tif)]
        all.tif = all.tif[grep("axe",all.tif)]
      }
      if(type_donnees == "axes_AFDM"){
        all.tif = list.files(paste0(gitCaractMilieu,"/output/ACP/"), 
                             pattern = ".tif",recursive = T,full.names = T)
        all.tif = all.tif[grep(mois,all.tif)]
        # Retirer dossiers : toutes, ACP_FAMD, CA/sans_ACP_clim
        all.tif = all.tif[!grepl("toutes",all.tif)]
        all.tif = all.tif[!grepl("ACP_FAMD",all.tif)]
        all.tif = all.tif[!grepl("CA/sans_ACP_clim",all.tif)]
        all.tif = all.tif[grep("axe",all.tif)]
      } 
      if(type_donnees == "brute"){
        all.tif = list.files(paste0(gitCaractMilieu,"/output/stack_dim/"),
                             pattern = ".tif$",recursive = T,full.names = T)
        all.tif = all.tif[grep(mois,all.tif)]
      }
      # Sauvegarder table des valeurs des rasters par mois
      if(!dir.exists(paste0(input_path,"/vars_env/",type_donnees,"/",mois))){
        dir.create(paste0(input_path,"/vars_env/",type_donnees,"/",mois),recursive=T)
      }
      DOS = paste0(input_path,"/vars_env/",type_donnees,"/",mois)
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
        if(type_donnees == "ACP" | type_donnees == "axes_FAMD"){
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
      # ajout coordonnées x y 
      df.crds <- terra::crds(stack.env.usage,df=TRUE)
      df.env.usage = cbind(df.env.usage,df.crds)
      # Arranger les noms, again
      if(type_donnees == "ACP"| type_donnees == "axes_FAMD"){
        # Remove nom mois dans colonnes
        to.modif = names(df.env.usage)[grep(mois,names(df.env.usage))]
        new.names = unlist(lapply(to.modif, function(x) substring(x,1,nchar(x) - nchar(mois) - 1)))
        names(df.env.usage)[grep(mois,names(df.env.usage))] = new.names
      }
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
  # type_donnees = "ACP" # "ACP" "axes_AFDM" "brute"
  
  
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
  
  # Make factor variable (with absence/presence name)
  data_glm$usage <- as.factor(data_glm$usage)
  data_glm$usage <- fct_recode(data_glm$usage,
             "presence" = "1", 
             "absence" = "0")
  str(data_glm)
  
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
  if(type_donnees == "axes_AFDM" | type_donnees == "ACP"){
    plot =  data_glm %>% subset(select=-c(x,y) ) %>%
      pivot_longer(!usage, names_to = "variables", values_to = "valeurs") %>%
      ggplot(aes(x=valeurs, y=variables, colour=as.factor(usage))) +
      geom_boxplot()+ 
      scale_colour_discrete(name = nom_beau)+
      theme(text = element_text(size=14))
  }
  # Sauvegarde
  png(file=paste0(output_path,"/niches/",type_donnees,"/",nom_court_usage,"/viz_distrib_famd.png"), 
      width=1400, height=800)
  print(plot)
  dev.off()
  
  # TODO : incorporer workflow BIOMOD2
  # - construire les modèles avec GLM et RF
  # - ensemble modeling : comparer chaque modèle + regroupement modèles
  # - prédictions : prédire avec le meilleur modèle (seul ou regroupement)
  
  # Modèle
  sample <- sample(c(TRUE, FALSE), nrow(data_glm), replace=TRUE, prob=c(0.7,0.3))
  train <- data_glm[sample, ]
  test <- data_glm[!sample, ] 
  
  # Weight absence obs (to give the same weight
  # to all absence and all presence observations)
  n_abs <- sum(train$usage == 0)
  n_pre <- sum(train$usage == 1)
  # n_abs+ n_pre == dim(train)[1]
  WEIGHT <-  n_pre / n_abs
  w.vect <- ifelse(train$usage == 0,WEIGHT,1)
  
  # Conditions de contrôle
  fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = 10,
    repeats = 10,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    savePredictions=TRUE)
  
  # Fit
  if(type_donnees == "brute" | type_donnees == "axes_AFDM"){
    
    # model.glm <- glm(usage ~ ., family=binomial, data=train, 
    #                  weights = w.vect)
    # model.glm <- stepAIC(model.glm)
    
    train <- susbset(train, select=-c("x","y"))
    
    # GLM
    model.glm <- caret::train(usage ~ .,
                              train,
                              method = "glm",
                              family = "binomial",
                              trControl = fitControl,
                              metric = 'ROC',
                              weights = w.vect)
  }
  if(type_donnees == "ACP"){
    
    # Formule
    formula.usage = as.formula(usage ~ axe1_toutes + axe2_toutes + 
                                 I(axe1_toutes^2)+ I(axe2_toutes^2) + #termes quadratiques
                                 axe1_toutes*axe2_toutes)

    # GLM
    model.glm <- caret::train(formula.usage,
                             train,
                             method = "glm",
                             family = "binomial",
                             trControl = fitControl,
                             metric = 'ROC',
                             weights = w.vect)
      
    #densityplot(model.glm, pch = "|")
    # # estimate of the uncertainty in our accuracy estimate
    # model.glm$results
    
    # summary(model.glm$finalModel)
    # model.glm.f <- model.glm$finalModel
    
    # TODO : inclure RF + comparaison modèle
    # # RF
    # model.rf <- caret::train(formula.usage,
    #                              train,
    #                              method = "ranger",
    #                              trControl = fitControl,
    #                              metric= 'ROC',
    #                              weights = w.vect)
    # 
    # # Comparaison entre modèles
    # resamps <- resamples(list(RF = model.rf,
    #                           GLM = model.glm))
    # summary(resamps)
    # trellis.par.set(caretTheme())
    # dotplot(resamps, metric = "ROC")
    

  }

  # Enregistrer le modèle pour pouvoir ensuite echo = F
  save(model.glm, file = paste0(output_path,"/niches/",type_donnees,
                                     "/",nom_court_usage,"/modele_glm.rdata"))
  
  # Find best threshold
  probs <- seq(.1, 0.9, by = 0.02)
  ths <- thresholder(model.glm,
                     threshold = probs,
                     final = TRUE,
                     statistics = "all")
  # ind_thr = which(ths[,"Balanced Accuracy"] == max(ths[,"Balanced Accuracy"]))
  # ths[ind_thr,"prob_threshold"]
  ind_thr = which(ths[,"J"] == max(ths[,"J"]))
  proba_threshold = ths[ind_thr,"prob_threshold"]
  
  # # Tester fiabilité modèle : score de Brier, matrice de confusion et AUC/ROC
  # brier_score = BrierScore(model.glm$finalModel)

  # Prédire sur jeu test pour tester accuracy + AUC
  
  # test$pred_prob <- predict(model.glm, test, type="response")
  # test$pred_resp <- ifelse(test$pred_resp_presence > 0.5, 1, 0)
  # test$pred_resp = as.factor(test$pred_resp)
  
  pred_prob <- predict(model.glm, test, type="prob")
  test <- cbind(test, pred_prob)
  test$pred_resp <- ifelse(test$presence > proba_threshold, 1, 0)
  test$pred_resp <- as.factor(test$pred_resp)
  test$pred_resp <- fct_recode(test$pred_resp,
                               "presence" = "1", 
                               "absence" = "0")
  
  # Sauvegarde fichier test
  write.csv(test,paste0(output_path,"/niches/",type_donnees,
                        "/",nom_court_usage,"/table_test.csv") )
  
  # Sauvegarde AUC/ROC
  png(file=paste0(output_path,"/niches/",type_donnees,
                  "/",nom_court_usage,"/roc.png"), 
      width=800, height=800)
  print(roc(test$usage ~ test$presence, plot = TRUE, print.auc = TRUE))
  #print(roc(test$usage ~ test$pred_resp, plot = TRUE, print.auc = TRUE))
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
      df.env2 = na.omit(df.env)
      #df.env$prob <- predict(model.glm, df.env, type="prob")
      
      prob_spatial <- predict(model.glm, df.env2, type="prob")
      
      #df.env$pred <- ifelse(df.env$prob > 0.50, 1, 0)
      
      raster_prob = rasterFromXYZ(data.frame(df.env2$x, df.env2$y, prob_spatial$presence), crs=EPSG_2154)
      
      #raster_prob = rasterFromXYZ(data.frame(df.env$x, df.env$y, df.env$prob), crs=EPSG_2154)
      
      raster_obs = raster(paste0(output_path,"/par_periode/",mois,"/",nom_lg,".tif"))
      # mask (car 0 confondu avec NA hors N2000)
      raster_obs <- mask(raster_obs, limiteN2000.shp)
      
      
      raster_obs <- projectRaster(raster_obs, raster_prob)
      
      all_rasters = stack(raster_obs, raster_prob)
      names(all_rasters) = c("observation","probabilite_presence")
      #plot(all_rasters, colNA='black')
      
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

# Fonction qui retourne des graphiques présentant la distribution 
# des usages (en densité de présence ou proba) dans un espace écologique
# (soit sur une AFDM d'une dimension, soit sur ACP de toutes les dimensions),
# pour une période donnée
EspEco <- function(usage, type_donnees, liste_mois, espace){
  # # TEST
  # liste_mois = "juillet" # Mois considéré (ou autre fenêtre temporelle)
  # usage =c("Rp") # Nom court de l'usage considéré, ex "Ni"
  # type_donnees = "brute" # "brute" ou "ACP", relatif aux variables utilisées
  # #pour construire les niches d'usages
  # espace = 'dimension' # "global" ou "dimension", projection dans un espace global ou par dimension
  
  # S'assurer ordre alphabétique (important pour la suite)
  usage = usage[order(usage)]
  
  if(!dir.exists(paste0(output_path,"/niches/",type_donnees,"/",paste0(usage,collapse="_"),"/esp_eco/"))){
    dir.create(paste0(output_path,"/niches/",type_donnees,"/",paste0(usage,collapse="_"),"/esp_eco/"), recursive = T)
  }
  
  outpath2 = paste0(output_path,"/niches/",type_donnees,"/",paste0(usage,collapse="_"),"/esp_eco/")
  stock_outpath2 = outpath2
  
  cat("\nUsage ",usage," :")
  
  PlotEspEco <- function(tableau_donnees_usages,
                         tableau_ACP,type){
    # # TEST
    # tableau_donnees_usages = values_usages
    # tableau_ACP = pca.vars.quanti
    # type = "quanti" # "quanti" "quali" "all"
    
    
    tableau_donnees_usages = na.omit(tableau_donnees_usages)
    
    tabeau_donnees_cercle <- circleFun(c(0,0),2,npoints = 500)
    
    if(type=="all"){
      tableau_ACP = tableau_ACP %>% select(-vars)
      names(tableau_ACP)[grep("vars_noms",names(tableau_ACP))] = "vars"
      tableau_for_plot = merge(tableau_ACP, dim_col, by.x="vars", by.y="vars_noms")
    }
    if(type=="quali"){
      tableau_ACP = tableau_ACP %>% select(-vars)
      names(tableau_ACP)[grep("vars_noms",names(tableau_ACP))] = "vars"
      tableau_for_plot = tableau_ACP
      col_of_dim = unique(col_vars$colour_dim[col_vars$dim == dim])
      
      tableau_for_plot$colour_dim = rep(col_of_dim,dim(tableau_for_plot)[1])
    }
    if(type=="quanti"){
      tableau_for_plot = merge(col_vars, tableau_ACP, by.x="Nom", by.y="vars_noms")
    }
    
    # Seuil pour garder variables sur graph (pour visibilité)
    seuil = 0.2
    tableau_for_plot.filtred = tableau_for_plot %>%
      filter(Dim.1 > seuil | Dim.1 < -seuil | Dim.2 > seuil | Dim.2 < -seuil)
    
    nb_var_drop = dim(tableau_for_plot)[1] -dim(tableau_for_plot.filtred)[1]
    cat(paste0("\n",nb_var_drop," variables ont été retirées pour faciliter la visualisation."))
    
    # Calculs coefficients pour afficher sur même plan indiv et variables
    Mu = max(max(tableau_donnees_usages$axe1), max(tableau_donnees_usages$axe2))
    Mv = max(max(tableau_for_plot.filtred$Dim.1),max(tableau_for_plot.filtred$Dim.2))
    RATIO = ifelse(Mu>Mv,round(Mu/Mv,2),1)
    RA2 = 0.8 * RATIO
    
    P3 = ggplot() +
      geom_path(data = tabeau_donnees_cercle,aes(x,y), 
                lty = 2, color = "grey", alpha = 0.7) +
      geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_hex(data=tableau_donnees_usages,
               aes(x=axe1, 
                   y=axe2)) +
      geom_segment(data = tableau_for_plot.filtred, 
                   aes(x = 0, xend = Dim.1*RATIO, y = 0, yend = Dim.2*RATIO),
                   arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
                   lwd = 0.6, colour="black") + 
      geom_text(data = tableau_for_plot.filtred, size = 10,
                aes(x = Dim.1*RA2, y =  Dim.2*RA2,
                    label = vars), 
                colour= tableau_for_plot.filtred$colour_dim,
                check_overlap = T) +
      labs(x="ACP 1",y="ACP2", fill="Densité\n(en pixels)",
           title= paste0("Usage ",paste0(usage,collapse="_")," - ",mois))+
      coord_equal() +
      scale_fill_continuous(type = "viridis",direction=-1)+
      theme_bw() +
      theme(panel.grid = element_blank(), text = element_text(size=15),
            panel.border = element_rect(fill= "transparent"))
    
    png(file=paste0(outpath2, "/densite_dim_",dim,"_",type,"_",
                    mois,".png"), 
        width=1400, height=800)
    print(P3)
    dev.off()
    
    # Graph du cercle de corrélation avec toutes les variables
    P1 = ggplot() +
      geom_path(data = tabeau_donnees_cercle,aes(x,y), 
                lty = 2, color = "grey", alpha = 0.7) +
      geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
      geom_segment(data = tableau_for_plot, 
                   aes(x = 0, xend = Dim.1*RATIO, y = 0, yend = Dim.2*RATIO),
                   arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
                   lwd = 0.6) + 
      geom_text(data = tableau_for_plot,size=10,
                aes(x = Dim.1*RA2, y =  Dim.2*RA2,
                    label = vars),
                colour=tableau_for_plot$colour_dim,
                check_overlap = F) +
      labs(x="ACP 1",y="ACP2", title= mois)+
      coord_equal() +
      theme_minimal() +
      theme(panel.grid = element_blank(), 
            text = element_text(size=15),
            panel.border = element_rect(fill= "transparent"))
    
    png(file=paste0(outpath2, "/cercle_correlation_",dim,"_",type,"_",
                    mois,".png"), 
        width=1400, height=800)
    print(P1)
    dev.off()
    
    if(nb_var_drop != 0){
      # Graph du cercle de corrélation avec les variables (> seuil)
      P2 = ggplot() +
        geom_path(data = tabeau_donnees_cercle,aes(x,y), 
                  lty = 2, color = "grey", alpha = 0.7) +
        geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
        geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
        geom_segment(data = tableau_for_plot.filtred, 
                     aes(x = 0, xend = Dim.1*RATIO, y = 0, yend = Dim.2*RATIO),
                     arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
                     lwd = 0.6) + 
        geom_text(data = tableau_for_plot.filtred,size=10,
                  aes(x = Dim.1*RA2, y =  Dim.2*RA2,
                      label = vars),
                  colour=tableau_for_plot.filtred$colour_dim,
                  check_overlap = F) +
        labs(x="ACP 1",y="ACP2", title= mois)+
        coord_equal() +
        theme_minimal() +
        theme(panel.grid = element_blank(), 
              text = element_text(size=15),
              panel.border = element_rect(fill= "transparent"))
      
      png(file=paste0(outpath2, "/cercle_correlation_seuil_",dim,"_",type,"_",
                      mois,".png"), 
          width=1400, height=800)
      print(P2)
      dev.off()}
  }
  
  for(mois in liste_mois){
    
    cat("\nMois de ",mois," en cours.")
    
    if(espace == "global"){
      ### Lecture des valeurs des variables (= axes AFDM ou variables initiales)
      data_mois <- fread(paste0(gitCaractMilieu,"/output/ACP/ACP_FAMD/",
                                df.mois$numero_mois[which(df.mois$nom_mois == mois)],mois,
                                "/tblFAMD_ACP_FAMD_",mois,".csv"),drop="V1")
      tbl_data = data_mois[,1:17]
      ### Calcul d'un analyse factorielle sur ces variables
      ana.fact <- PCA(tbl_data, graph = FALSE)
      ### Mise en forme dt de l'analyse factorielle
      # https://tem11010.github.io/Plotting-PCAs/
      tbl_data$pc1 <- ana.fact$ind$coord[, 1]
      tbl_data$pc2 <- ana.fact$ind$coord[, 2]  
      pca.vars <- ana.fact$var$coord %>% data.frame
      pca.vars$vars <- rownames(pca.vars)
      # https://stackoverflow.com/questions/51219267/pca-scaling-not-applied-to-individuals
      # https://rdrr.io/cran/factoextra/src/R/fviz_pca.R#sym-fviz_pca_biplot
      N = nchar(mois) + 1
      pca.vars$vars_noms = substr(pca.vars$vars,1, nchar(pca.vars$vars)-N)
      
      #pca.vars.m <- melt(pca.vars, id.vars = "vars",warning =F)
      
      ### Lecture des usages (probabilité ou présence/absence)
      liste_rast = list.files(paste0(output_path,"/niches/",type_donnees,"/",
                                     usage,"/predictions/"),".tif$",
                              full.names = T)
      if(length(liste_rast[grep(mois, liste_rast)])==0){
        cat("\nPas d'usage ",usage,"au mois de ",mois,".")
      }else{
        r = stack(liste_rast[grep(mois, liste_rast)])
        names(r)=c("observation", "probabilite_prediction")
        obs_sp = rasterToPolygons(r$observation,fun=function(x)x==1, dissolve=TRUE)
        # AXES ACP
        r.vars <- stack(list.files(paste0(gitCaractMilieu,"/output/ACP/ACP_FAMD/",
                                          df.mois$numero_mois[which(df.mois$nom_mois == mois)],mois,
                                          "/"),".tif", full.names = T))
        # EXTRACTION VALEURS
        values_usages <- as.data.frame(exact_extract(r.vars, obs_sp, include_xy=T))
        values_usages$axe1 = values_usages[,1]
        values_usages$axe2 = values_usages[,2]
        
        # pca.vars.2 = merge(pca.vars.2, dim_col, by="vars_noms")
        
        ### Construction des figures
        dim= 'AFDM'
        
        PlotEspEco(values_usages,pca.vars, "all")
        
      }
    }
    
    if(espace == "dimension"){
      ### Lecture des valeurs des variables (= axes AFDM ou variables initiales)
      pot_fil = list.files(paste0(gitCaractMilieu, "/output/ACP/"),mois,recursive = T, full.names = T)
      pot_fil = pot_fil[grep("tblFAMD_",pot_fil)]
      # à virer : toutes ACP_FAMD TEST/CA_avec_ACP_clim
      path_files = pot_fil[!grepl("TEST|ACP_FAMD|toutes",pot_fil)]
      # Liste de df, 1 par dimension
      liste_tbl_data <- lapply(path_files, function(x) fread(x,drop="V1"))
      ### Calcul d'un analyse factorielle sur ces variables
      for(i in 1:length(liste_tbl_data)){
        
        dim = liste.dim[order(liste.dim)][i]
        cat("\n Dimension",dim ,"en cours.")
        
        outpath2 = paste0(stock_outpath2,dim)
        
        if(!dir.exists(outpath2)){
          dir.create(outpath2, recursive = T)
        }
        
        data_mois = liste_tbl_data[[i]]
        ind_x = grep("^x$",names(data_mois))
        tbl_data = data_mois[,1:(ind_x-1)]
        
        # Recoder vars factorielles
        tbl_data = RecodeFacto(tableau_donnees_a_recoder = tbl_data,
                               tableau_nature_donnees=col_vars_formate)
        t = try(FAMD(tbl_data, graph = FALSE))
        if(inherits(t, "try-error")) {
          # PCA si seulement quanti
          ana.fact <- PCA(tbl_data , graph = FALSE)
          
          pca.vars <- ana.fact$var$coord %>% data.frame
          pca.vars$vars <- rownames(pca.vars)
          pca.vars$vars_noms = pca.vars$vars
        } else{
          #FAMD si miste quali/quanti
          ana.fact <- FAMD(tbl_data , graph = FALSE)
          
          pca.vars.quanti <- ana.fact$quanti.var$coord %>% data.frame
          pca.vars.quanti$vars <- rownames(pca.vars.quanti)
          pca.vars.quanti$vars_noms = pca.vars.quanti$vars
          
          pca.vars.quali <- ana.fact$quali.var$coord %>% data.frame
          # Mais quelle bidouille de l'enfer pour récupérer les noms
          # des variables quaitatives dans l'AFDM .....
          Y = data.frame(ana.fact$call$quali.sup$quali.sup)
          vars_nom = c()
          for(j in 1:dim(Y)[2]){
            vars_nom = c(vars_nom,paste0(names(Y)[j],levels(Y[,j])))
          }
          pca.vars.quali$vars <- vars_nom
          pca.vars.quali$vars_noms = pca.vars.quali$vars
        }
        
        ### Mise en forme dt de l'analyse factorielle  
        tbl_data$pc1 <- ana.fact$ind$coord[, 1]
        tbl_data$pc2 <- ana.fact$ind$coord[, 2]
        
        ### Lecture des usages (probabilité ou présence/absence)
        liste_rast = list.files(paste0(output_path,"/niches/",type_donnees,"/",
                                       usage,"/predictions/"),".tif$",
                                full.names = T)
        
        if(length(liste_rast[grep(mois, liste_rast)])==0){
          cat("\nPas d'usage ",usage,"au mois de ",mois,".")
        }else{
          r = stack(liste_rast[grep(mois, liste_rast)])
          names(r) = unlist(lapply(usage, function(x) paste0(x,c("_observation", "_probabilite_prediction"))))
          
          # Calcul espace éco en prédit et en obs
          r_obs <- raster::subset(r, grep('observation', names(r)))
          r_pred <- raster::subset(r, grep('prediction', names(r)))
          
          # si travail sur un seul usage, ne sert à rien mais n'affecte pas
          r_obs <- sum(r_obs)
          names(r_obs) = "observation"
          # # a voir comment on traite l'overlay de proba
          # r_pred <- prod(r_pred)
          
          # Checker si multiusage ou pas
          N = ifelse(length(usage)>1,length(usage),1)
          obs_sp = rasterToPolygons(r_obs,fun=function(x)x==N, dissolve=TRUE)
          
          # AXES des analyses factorielles (ACP ou AFDM)
          A = list.files(paste0(gitCaractMilieu,"/output/ACP/",dim,"/",
                                paste0(df.mois$numero_mois[which(df.mois$nom_mois == mois)],mois),"/"),
                         ".tif$", full.names = T)
          A = A[grepl("axe",A)]
          r.vars <- stack(A)
          # Valeurs des obs/pred dans les axes factoriels
          values_usages <- as.data.frame(exact_extract(r.vars, obs_sp, include_xy=T))
          
          values_usages$axe1 = values_usages[,1]
          values_usages$axe2 = values_usages[,2]
          ### Construction des figures
          
          if(inherits(t, "try-error")) {
            PlotEspEco(values_usages,pca.vars,"quanti")
          } else  {          
            PlotEspEco(values_usages,pca.vars.quanti,"quanti")
            PlotEspEco(values_usages,pca.vars.quali,"quali")
          }
        }
      }
    }
  }
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
# formatage pour que cette table soit prise dans fonction RecodeFacto()
col_vars_formate =  col_vars %>%  
  dplyr::select(Nom,Nature) %>%
  transmute(ID=Nom, nature=Nature)

#### GLM : binomiale ####
# axes des AFDM par dimension

# Création d'un df par usage : utilisation fonction ExtractData1Use
lapply(c("nidification",
         "couchade",
         "paturage",
         "randonnee_pedestre",
         "VTT",
         "parade"),function(x) ExtractData1Use(usage=x, type_donnees = "brute"))

lapply(c("nidification",
         "couchade",
         "paturage",
         "randonnee_pedestre",
         "VTT",
         "parade"),function(x) ExtractData1Use(usage=x, type_donnees = "ACP"))

lapply(c("nidification",
         "couchade",
         "paturage",
         "randonnee_pedestre",
         "VTT",
         "parade"),function(x) ExtractData1Use(usage=x, type_donnees = "axes_FAMD"))


# Exploitation des données : création modèle
set.seed(1)
lapply(liste.usages, function(x) CreateModelUsage(nom_court_usage=x,type_donnees = "brute"))

lapply(c("Rp","Pa"), function(x) CreateModelUsage(nom_court_usage=x,type_donnees = "ACP"))

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

# Tests de la fonction :
 # - pour 1 usage, pour un mois --> OK
EspEco(usage="Rp",type_donnees = "brute",liste_mois="juillet",espace="global")
 # - pour 1 usage, sur une liste de mois-->OK
EspEco(usage="Rp",type_donnees = "brute",liste_mois=liste.mois,espace="global")
 # - pour tous les usages, sur une liste de mois --> OK
lapply(liste.usages, function(x) EspEco(usage = x, 
                                        liste_mois = liste.mois,
                                           type_donnees = "brute",
                                           espace="global"))
 # - pour un multi-usage, sur un mois --> OK
EspEco(usage=c("Pa","Vt","Ni"),type_donnees = "brute",liste_mois="juillet",espace="dimension")

 # - pour un multi-usage, sur une liste de mois 
# --> ne marche pas quand un des 3 pas présents sur un mois
EspEco(usage=c("Pa","Vt","Ni"),type_donnees = "brute",liste_mois=liste.mois,espace="dimension")


# pour 1 usage, sur un mois, en dimension --> OK
EspEco(usage="Rp",type_donnees = "brute",liste_mois= "juillet",espace="dimension")

 # pour 1 usage, sur une liste de mois, en dimension --> OK
EspEco(usage="Pa",type_donnees = "brute",liste_mois=liste.mois,espace="dimension")
