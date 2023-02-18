#Packages----
library(raster)
library(terra)
library(sf)
library(magrittr)
library(purrr)
library(stars)
library(rmapshaper)
source("./R/funciones.R")

#Seleccionar Mallas ----
malla_wgs <- raster("./input/Malla_neotropical_wgs.tif")#wgs84
malla_lcc <- raster("./input/Malla_neotropical.tif")#LCC

#Open SDM ----
list_mdp <-list.files("./input/SDM_Ureta_et_al_2022", pattern = ".tif$", full.names = TRUE)

### Tapir sdm
Tapir <- raster(list_mdp[grep("Tapirus", basename(list_mdp))])
Tapir <- crop(Tapir, malla_wgs, datatype='INT1S') %>% mask(., malla_wgs, datatype='INT1S') %>% projectRaster(., malla_lcc, method="ngb")
dataType(Tapir) <- 'INT1S'

### Remove tapir in the list
list_mdp <- list_mdp[!grepl("Tapirus", basename(list_mdp))]

## Overlap Tapir vs sdm----
pb <- txtProgressBar(0, length(list_mdp), style = 3)
result_overlap <- lapply(1:length(list_mdp), function(x){
  setTxtProgressBar(pb, x)
  spp <- crop(raster(list_mdp[x]), malla_wgs, datatype='INT1S') %>% mask(., malla_wgs, datatype='INT1S') %>% projectRaster(., malla_lcc, method="ngb")
  spp <- Tapir + spp; dataType(spp) <- 'INT1S'; return(spp)
  })
saveRDS(result_overlap, file = "./output/result_overlap_Tapir_vs_spp.Rds")

## Distr. Tapir en PA ----
PA <- read_sf("./input//PAs_Mexico_todas.shp")
#% inside PAs
Tapir_PA <- mask(Tapir, PA)
#Under protection
freq(Tapir_PA)
#1 = 633
#Total Tapir distribution area
freq(Tapir)
#1 = 1928
#Percentage protected
round((633*100)/1928)
#33%

## No especies con intersect----
base_area <- (res(result_overlap[[1]])[1]^2)* 1e-6 #km2
pb <- txtProgressBar(0, length(result_overlap), style = 3)
Interseccion_M1 <- map_dfr(1:length(result_overlap), function(x){
  setTxtProgressBar(pb, x)
  x.1 <- result_overlap[[x]]
  if(!is.null(x.1)){
    if(2 %in% unique(raster::values(x.1))){
      x.1 <- freq(x.1) %>% as.data.frame()
      x.2 <- data.frame(especie = basename(list_mdp[x]) %>% gsub(".tif", "", .) %>% gsub("\\.", " ", .),
                        intersect = "Yes",
                        "Overlap (km2)" = round(x.1[which(x.1[[1]] == 2), 2]*base_area, 2),
                        "% overlap" = (round(x.1[which(x.1[[1]] == 2),2]*base_area,2) * 100)/(sum(round(x.1[which(x.1[[1]] == 2 | x.1[[1]] == 1), 2] * base_area))), 
                        "No overlap (km2)" = round(x.1[which(x.1[[1]] == 1),2]*base_area, 2),
                        check.names = FALSE)
    } else {
      x.1 <- freq(x.1) %>% as.data.frame()
      x.2 <- data.frame(especie = basename(list_mdp[x]) %>% gsub(".tif", "", .) %>% gsub("\\.", " ", .),
                        intersect = "No",
                        "Overlap (km2)" = 0,
                        "% overlap" = 0,
                        "No overlap (km2)" = round(x.1[which(x.1[[1]] == 1), 2] * base_area, 2),
                        check.names = FALSE)
    }  
  } else {
    x.1 <- freq(x.1) %>% as.data.frame()
    x.2 <- data.frame(especie = basename(list_mdp[x]) %>% gsub(".tif", "", .) %>% gsub("\\.", " ", .),
                      intersect = "No",
                      "Overlap (km2)" = 0,
                      "% overlap" = 0,
                      "No overlap (km2)" = round(x.1[which(x.1[[1]] == 1), 2] * base_area, 2),
                      check.names = FALSE)
  }
  return(x.2)
})
head(Interseccion_M1)
mean(Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap` > 0)])
length(Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap` > 0)])

###No species intersect >10%
length(Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap` >= 10)])
Interseccion_M1$especie[Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap` >= 10)]]

boxplot(Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap`>= 10)])
saveRDS(Interseccion_M1, file = "./output/table_overlap_10Perc_Tapir_vs_spp.Rds")
Interseccion_M1[which(Interseccion_M1$`% overlap` >= 10)]
length(which(Interseccion_M1$`% overlap` >= 10))


## Potential species richness----
pb <- txtProgressBar(0, length(result_overlap), style = 3)
result_S <- lapply(which(Interseccion_M1$`% overlap` >= 10), function(x){
  setTxtProgressBar(pb, x)
  if(2 %in% unique(result_overlap[[x]][])){
    x<- result_overlap[[x]]; x[x != 2] <- NA;x[x == 2] <- 1 
  } else {
    x <- NULL
  }
  return(x)
})
result_S <- compact(result_S); result_S$fun <- sum; result_S <- do.call(mosaic, result_S)
result_S[result_S == 0] <- NA; dataType(result_S) <- 'INT2U'
writeRaster(result_S, "./output/Species_TapirDist.Tif", 
            overwrite = TRUE, datatype='INT2U')

#Grafica de especies en IUCN ----
#This data table was taked from the supplementary material of Ureta et al. 2022
IUCN <- read.csv("./input/Generalidades_estadoconservacion.csv")
Interseccion_M1$IUCN <- map_chr(Interseccion_M1$especie, function(x){
  x.1 <- IUCN$Risk.category..IUCN.[which(IUCN$Scientific.name == x)]
  
  if(length(x.1) == 0){
    x.1 <- "Unclassified"
  } else {
    if(x.1 == ""){
      x.1 <- "Unclassified"
    }
  }
  
  return(x.1)
})

Interseccion_M1$`intersect_10%`<- "No"
Interseccion_M1$`intersect_10%`[which(Interseccion_M1$`% overlap` >=10)]<- "Yes"

Interseccion_M1_IUCN <- Interseccion_M1[which(Interseccion_M1$`intersect_10%` == "Yes"),]
Interseccion_M1_IUCN <- merge(Interseccion_M1_IUCN, IUCN, by.x = "especie", by.y = "Scientific.name")
names(Interseccion_M1_IUCN)
Interseccion_M1_IUCN$especie

IUCN <- table(Interseccion_M1_IUCN$IUCN) %>% as.data.frame()
IUCN$Perc <- round(IUCN$Freq*100/ 43)
IUCN

Risk <- table(Interseccion_M1_IUCN$Population.IUCN.SEMARNAT.NatureServe) %>% as.data.frame()
Risk$Perc <- round(Risk$Freq*100/ 43)
Risk
saveRDS(Interseccion_M1_IUCN, "./output/IUCN_Tapir_10p.Rds")

#Accumulation curve----
Tapir_shp <- Tapir; Tapir_shp[Tapir_shp == 0] <- NA
plot(Tapir_shp)
Tapir_shp <- terra::as.polygons(terra::rast(Tapir_shp), dissolve=FALSE) %>% as(.,"Spatial") %>% st_as_sf()
Tapir_shp$id <- 1:nrow(Tapir_shp)

id_azar <-  sample(factor(rep(1:100, length.out= nrow(Tapir_shp))))
lista_mdp2 <- st_as_stars(stack(result_overlap))

spaccum_Tapir <- sp_accum(id_azar_n = id_azar, Tapir_shp = Tapir_shp, 
                          models_list = lista_mdp2, 
                          names_list = Interseccion_M1$especie)
saveRDS(spaccum_Tapir, file = "./output/spaccum_Tapir_1.rds")

