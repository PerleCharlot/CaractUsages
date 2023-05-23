### Titre -------------------------------------
# Nom : exploration
# Auteure : Perle Charlot
# Date de création : 27-04-2023
# Dates de modification : 09-05-2023

### Librairies -------------------------------------
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(reshape2)
library(viridis)
library(tidyverse)
library(ggh4x)
### Fonctions -------------------------------------

# Fonction qui calcule le % de pixels (géographiques) d'un usage (u1) partagés avec chaque autre usage
# avec en référence l'extent de la zone d'étude ou l'extent de l'usage (u1)
Comp_p_tens <- function(data,use1,use2, seq){
  # # TEST
  # data <-  df_proba_PCA_glob_sub2
  # use1 <- "Rp_proba"
  # use2 <- "Vt_proba"
  # seq <- seq(0,1,0.01)
  
  df_test <- data %>%
    filter(Use == use1 | Use == use2) %>%
    pivot_wider(names_from = Use,values_from = Proba)
  ComputePourcentTension <- function(s){
    # # TEST
    # s <- seq_seuil
    
    n_pix_u1 <- dim(df_test[df_test[use1] > s,])[1]
    
    n_pix_tension <- dim(df_test[df_test[use1] > s & df_test[use2] > s,])[1]
    n_pix_studysite <- dim(df_test)[1] # pourcentage par rapport à la taille de la zone d'étude
    # et non pas par rapport à l'extent de l'usage
    p_tens_studysite <- (n_pix_tension/n_pix_studysite)*100
    p_tens_u1 <- (n_pix_tension/n_pix_u1)*100
    
    surf_1pix_km2 <- (25*25) / 1000000
    
    return(data.frame("p_studysite" = p_tens_studysite,
                      "p_u1" = p_tens_u1,
                      "u1"=use1,
                      "u2"=use2,
                      "seuil"= s,
                      "area_u1_km2"= n_pix_u1 *surf_1pix_km2))
  }
  
  df_p <- do.call(rbind,(lapply(seq, ComputePourcentTension))) 
  
  return(df_p)
}

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

### ESPACE GEOGRAPHIQUE ####

# Pour moi, Schoener dans l'espace géographique ça calcule la similarité de distribution pourça :
df_proba_PCA_glob_sub2 <- df_proba_PCA_glob_sub %>%
  pivot_wider(names_from = Use,values_from = Proba) %>%
  group_by(x,y,Month) %>%
  mutate(num_pixel = 1)
df_proba_PCA_glob_sub2$num_pixel <- 1:dim(df_proba_PCA_glob_sub2)[1]
df_proba_PCA_glob_sub2 <- df_proba_PCA_glob_sub2 %>%
  pivot_longer(cols=ends_with("proba"),
               names_to = "Use" ,values_to = "Proba")
# Pour moi, Schoener dans l'espace géographique ça calcule la similarité de distribution pour ça :
# où chaque pixel (x) représente un "site" avec une certaine "abondance" (= probabilité ici)
P <- ggplot(data=df_proba_PCA_glob_sub2, aes(x=num_pixel))+
  geom_point(aes(y=Proba, col=Use),alpha=ifelse(df_proba_PCA_glob_sub2$Proba>0.75,1,0.1))+
  facet_grid(Use~.,labeller = labeller(Use = label_test))+
  labs(y="Predicted Probability", x="Pixel Number",subtitle = "All pixels")+
  geom_hline(yintercept = 0.75)+
  guides(col = FALSE)+
  theme_bw()
# Zoom sur les mille premiers pixels
Pzoom <- df_proba_PCA_glob_sub2 %>%
  filter(num_pixel<1000) %>%
ggplot( aes(x=num_pixel))+
  geom_hline(yintercept = 0.75)+
  geom_point(aes(y=Proba, col=Use), alpha=0.5)+
  facet_grid(Use~.,labeller = labeller(Use = label_test))+
  guides(col = FALSE)+
  labs(y="Predicted Probability",x="Pixel Number", subtitle = "Sample: pixel 1 to 1000")+
  theme_bw()
PP <- P + Pzoom
PP
png(file = paste0(Gspace_path,"/probability_by_site_all_uses_",example_month,".png"),width=1400, height=800)
plot(PP)
dev.off()

# Décomposé par latitude et longitude (préférence d'habitat + disponibilité environnement en fond)
Px <- ggplot(data=df_proba_PCA_glob_sub, aes(x=x))+
  geom_histogram(fill="grey", bins=500)+
  geom_point(aes(y=Proba*200, col=Use), alpha=0.05)+
  #geom_line(aes(y=Proba*200, color=Use), alpha=0.5)+
  facet_grid(~Use,labeller = labeller(Use = label_test))+
  scale_y_continuous(name = "Pixel Count",
                     sec.axis = sec_axis(~./200, name="Predicted Probability"))+
  guides(col = FALSE)+
  theme_bw()
Py <- ggplot(data=df_proba_PCA_glob_sub, aes(x=y))+
  geom_histogram(fill="grey", bins=500)+
  #geom_line(aes(y=Proba*200, color=Use), alpha=0.5)+
  geom_point(aes(y=Proba*200, color=Use), alpha=0.05)+
  facet_grid(~Use,labeller = labeller(Use = label_test))+
  scale_y_continuous(name = "Pixel Count",
                     sec.axis = sec_axis(~./200, name="Predicted Probability"))+
  guides(col = FALSE)+
  theme_bw()+coord_flip()
map <-  ggplot(data=df_proba_PCA_glob_sub) +
  geom_raster(aes(x = x, y = y, fill = Proba)) +
  scale_fill_viridis(limits=c(0,1)) +
  facet_grid(.~Use,labeller = labeller(Use = label_test),drop=T)+
  guides(fill = FALSE)+
  #coord_equal()+
  theme_bw()
plot_explo_geo <- map + Px  + Py + plot_layout(nrow=3,ncol=1)

png(file = paste0(Gspace_path,"/probability_by_xy_all_uses_",example_month,".png"),width=1400, height=800)
plot(plot_explo_geo)
dev.off()


# Identification des pixels à forte tension spatiale = CHEVAUCHEMENT GEOGRAPHIQUE
# = proba forte + chevauchement avec autre usage
seq_seuil <- seq(0,1,0.01)

liste_use1 <- unique(df_proba_PCA_glob_sub2$Use)
df0 <- data.frame()
for(i in 1:length(liste_use1)){
  use1 <- liste_use1[i]
  names(use1) <- names(uses.proba[uses.proba == use1])
  liste_use2 <- liste_use1[-i]
  
  dd <- do.call(rbind,(lapply(liste_use2,
                          function(x) Comp_p_tens(data = df_proba_PCA_glob_sub2,
                                                  use1 = use1,
                                                  use2 = x,
                                                  seq = seq_seuil))))
  
  #names(dd)  <- liste_use2 
  #df <- cbind(seuil=seq_seuil, dd)
  # df <- df %>%
  #   pivot_longer(cols=ends_with("proba"),
  #                names_to = "Use2" ,values_to = "p_tens")
  # df$Use1 <- use1
  
  df0 <- rbind(df0,dd)
}

# Sur 1 même graphique, extent de chaque usage + % pixel overlap (par rap extent zone d'étude et extent u1)
df0bis <- df0 %>%
  pivot_longer(cols=c("p_studysite","p_u1"),
               names_to = "y" , values_to = "y_values")
label_y<- c("% Overlapping Pixels - Study Area","% Overlapping Pixels - Use")
names(label_y) <- c("p_studysite","p_u1")

df0bis$u2 <- as.factor(df0bis$u2)

pchv <- df0bis %>% ggplot(aes(x=seuil, y=y_values, color=u2)) +
  geom_line()+geom_point(size=1)+
  geom_hline(yintercept = 50, linetype='dashed')+
  geom_vline(xintercept = 0.5, linetype='dashed')+
  labs( x="Probability Threshold",
       title="Geographical Overlapping Intensity",
       subtitle=example_month)+
  facet_grid(y~u1,labeller = labeller(u1 = label_test, y=label_y), scales = "free")+
  #force_panelsizes(rows = c(0.2,0.5,0.8))+
  theme(axis.title.y = element_blank())+
scale_color_discrete(name = "Use", labels = label_test)

# représente l'extent de u1 (~ prévalence ??)
pchv2 <- df0 %>% ggplot(aes(x=seuil, y=area_u1_km2)) +
  geom_line() +
  geom_ribbon(aes(ymax=area_u1_km2),ymin=0,alpha=0.3) +
  facet_wrap(.~u1,labeller = labeller(u1 = label_test),nrow = 1)+
  labs(x="Probability Threshold",
       y="Area (km²)") +
  theme_bw()

design <- "
  111111
  111111
  111111
  111111
  222222
"

P_chvch <- pchv + pchv2 + plot_layout(design = design)
P_chvch 

png(file = paste0(Gspace_path,"/overlapping_intensity_",example_month,".png"),width=2200, height=1200)
plot(P_chvch)
dev.off()


# analyser similarité de niche éco


### ESPACE ECOLOGIQUE ####

design <- "
  111111
  222222
"
Ppca1 + Ppca2 + plot_layout(design = design)

df_proba_PCA_glob_sub %>%
  pivot_longer(cols=c("PCA1","PCA2"), names_to = "Axe",values_to = "axe_value") %>%
  ggplot(aes(x=axe_value))+
  geom_histogram(fill="grey", bins=500)+
  geom_point(aes(y=Proba*400, col=Use), alpha=0.05)+
  facet_grid(Axe~Use,labeller = labeller(Use = label_test), scales = "free")+
  labs(x="Environmental Axe")+
  scale_y_continuous(name = "Pixel Count",
                     sec.axis = sec_axis(~./400, name="Predicted Probability"))+
  guides(col = FALSE)+
  theme_bw()


p_med_pca1 <- df_t2 %>% as.data.frame() %>%
  ggplot(aes(x=cut_x,y=median_Proba))+
  geom_point()+
  facet_grid(~Use,labeller = labeller(Use = label_test))

P_chelou <- Ppca1 + p_med_pca1 + plot_layout(design = design)
png(file = paste0(Espace_path,"/niche_similarit_test_",example_month,".png"),width=2200, height=1200)
plot(P_chelou)
dev.off()

Ppca1 <- ggplot(data=df_proba_PCA_glob_sub2, aes(x=PCA1))+
  geom_histogram(fill="grey", bins=500)+
  geom_point(aes(y=Proba*400, col=Use), alpha=0.05)+
  #geom_line(aes(y=Proba*200, color=Use), alpha=0.5)+
  facet_grid(~Use,labeller = labeller(Use = label_test))+
  scale_y_continuous(name = "Pixel Count",
                     sec.axis = sec_axis(~./400, name="Predicted Probability"))+
  guides(col = FALSE)+
  theme_bw()
Ppca2 <- ggplot(data=df_proba_PCA_glob_sub, aes(x=PCA2))+
  geom_histogram(fill="grey", bins=500)+
  #geom_line(aes(y=Proba*200, color=Use), alpha=0.5)+
  geom_point(aes(y=Proba*400, color=Use), alpha=0.05)+
  facet_grid(~Use,labeller = labeller(Use = label_test))+
  scale_y_continuous(name = "Pixel Count",
                     sec.axis = sec_axis(~./400, name="Predicted Probability"))+
  guides(col = FALSE)+
  theme_bw() #+coord_flip()
map2 <- ggplot(data=df_proba_PCA_glob_sub) +
  stat_summary_hex(aes(x=PCA1, y=PCA2, z= Proba),
                   fun = function(x) median(x), bins=75,colour='transparent')+
  scale_fill_viridis(na.value = "transparent",
                     limits=c(0,1))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_vline(xintercept = 0,linetype="dashed")+
  facet_grid(.~Use,labeller = labeller(Use = uses.labs),drop=T)+
  guides(fill = FALSE)+
  theme_bw()
plot_explo_eco <- map2 + Ppca1  + Ppca2 + plot_layout(nrow=3,ncol=1)
plot_explo_eco
png(file = paste0(Espace_path,"/probability_by_pca1pca2_all_uses_July.png"),width=1400, height=800)
plot(plot_explo_eco)
dev.off()

# Grille de 150 * 150 (résolution pixels écologiques)
taille_maille <- 10
df_t <- df_proba_PCA_glob_sub %>%
  dplyr::select(Use, Proba, PCA1, PCA2) %>%
  mutate(cut_x = cut(PCA1, breaks = round(seq(from = min(PCA1, na.rm = T), 
                                              to = max(PCA1, na.rm = T), 
                                              length.out = taille_maille),2),
                     include.lowest = T),
         cut_y = cut(PCA2, breaks = round(seq(from = min(PCA2, na.rm = T),
                                              to = max(PCA2, na.rm = T),
                                              length.out = taille_maille),2),
                     include.lowest = T)) %>%
  group_by(cut_x, cut_y,Use, .drop = F) %>% 
  summarise(n_pixels = n(), 
            med_proba = median(Proba,na.rm=T),
            Use=Use,
            num_combi_env = 1)%>%
  distinct() %>%
  pivot_wider(names_from = Use, values_from = med_proba) %>%
  mutate_all( ~replace_na(.,0))

df_t$num_combi_env <- 1:dim(df_t)[1]
df_t2 <- df_t %>%
  pivot_longer(cols=ends_with("proba"),
               names_to = "Use" ,values_to = "median_Proba")

Pe <- ggplot(data=df_t2, aes(x=num_combi_env))+
  geom_point(aes(y=Proba, col=Use),alpha=ifelse(df_t2$Proba>0.75,1,0.1))+
  facet_grid(Use~.,labeller = labeller(Use = uses.labs))+
  labs(y="Predicted Probability", x="Environmental Pixel Number",subtitle = "All pixels")+
  geom_hline(yintercept = 0.75)+
  guides(col = FALSE)+
  theme_bw()
Pzoome <- df_t2 %>%
  filter(num_combi_env<1000) %>%
  ggplot( aes(x=num_combi_env))+
  geom_hline(yintercept = 0.75)+
  geom_point(aes(y=Proba, col=Use), alpha=0.5)+
  facet_grid(Use~.,labeller = labeller(Use = uses.labs))+
  guides(col = FALSE)+
  labs(y="Predicted Probability",x="Environmental Pixel Number", subtitle = "Sample: pixel 1 to 1000")+
  theme_bw()

# TODO : remove x axis text

coef <- max(log(df_t2$n_pixels))
Pe2 <- ggplot(data=df_t2, aes(x=reorder(num_combi_env, -n_pixels)))+
  geom_col(aes(y=log(n_pixels)),col="grey")+
  geom_point(aes(y=Proba*coef, col=Use), alpha=ifelse(df_t2$Proba>0.75,1,0.3)) +
  facet_grid(Use~.,labeller = labeller(Use = uses.labs), scales = "free")+
  geom_hline(yintercept = 0.75*coef)+
  scale_y_continuous(name = "Log(Geographic Pixels Count)",
                     sec.axis = sec_axis(~./coef, name="Predicted Probability"))+
  labs(x="Environmental Pixel Number",subtitle = "All pixels (sorted)")+
  guides(col = FALSE)+
  theme_bw()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank())

# Garder les 1000 lignes
test <- df_t2 %>%
  group_by(Use) %>%
  # arrange(-n_pixels) %>%
  slice_max(n_pixels, n=1000)

Pzoome2 <- test %>% 
  ggplot(aes(x=reorder(num_combi_env, -n_pixels)))+
  geom_col(aes(y=log(n_pixels)),fill="grey")+
  geom_point(aes(y=Proba*coef, col=Use), alpha=0.5) +
  facet_grid(Use~.,labeller = labeller(Use = uses.labs), scales = "free")+
  geom_hline(yintercept = 0.75*coef)+
  scale_y_continuous(name = "Log(Geographic Pixels Count)",
                     sec.axis = sec_axis(~./coef, name="Predicted Probability"))+
  labs(x="Environmental Pixel Number",subtitle = "Sample: first 1000 pixels")+
  guides(col = FALSE)+
 theme_bw()+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank())


PPe <- Pe + Pzoome
PPe2 <- Pe2 + Pzoome2

png(file = paste0(Espace_path,"/probability_by_site_all_uses_July.png"),width=1400, height=800)
plot(PPe)
dev.off()

png(file = paste0(Espace_path,"/probability_by_site_all_uses_disponib_July.png"),width=1400, height=800)
plot(PPe2)
dev.off()


#count pixls géo dans espace éco définit seulement par PCA1
ggplot(data=df_t2, aes(x=cut_x))+
  geom_point(aes(y=log(n_pixels)),col="grey")+
  facet_grid(Use~.,labeller = labeller(Use = uses.labs), scales = "free")

ggplot(data=df_t2, aes(x=cut_x))+
  geom_col(aes(y=log(n_pixels)),col="grey")+
  geom_point(aes(y=Proba*coef, col=Use), alpha=ifelse(df_t2$Proba>0.75,1,0.3)) +
  facet_grid(Use~.,labeller = labeller(Use = uses.labs), scales = "free")+
  geom_hline(yintercept = 0.75*coef)+
  scale_y_continuous(name = "Log(Geographic Pixels Count)",
                     sec.axis = sec_axis(~./coef, name="Predicted Probability"))+
  labs(x="Environmental Pixel Number",subtitle = "All pixels (sorted)")+
  guides(col = FALSE)+
  theme_bw()+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank())