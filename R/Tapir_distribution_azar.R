#Paquetes----
library(raster)
library(terra)
library(sf)
library(magrittr)
library(purrr)
library(stars)
library(landscapeR)
library(rmapshaper)
source("./R/funciones.R")

# Modelo distribución azar ----
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

##Tapir azar ----
raster_tapir <- malla_lcc; raster_tapir <- raster_tapir * 0; raster_tapir <- raster_tapir + Tapir
raster_tapir[raster_tapir == 1] <- NA

# num <- 20
# size <- rand_vect(N = num, M = freq(Tapir)[which(freq(Tapir)[,1] == 1),2], sd = 20, pos.only = TRUE, seed = 1)
# Tapir_azar <- makeClass(raster_tapir, num, size); plot(Tapir_azar)
#saveRDS(Tapir_azar, "./output/Sel_area_randomly.Rds")
Tapir_azar <- readRDS("./output/Sel_area_randomly.Rds")

## Overlap random patches vs sdm----
# pb <- txtProgressBar(0, length(list_mdp), style = 3)
# result_overlap_random <- lapply(1:length(list_mdp), function(x){
#   setTxtProgressBar(pb, x)
#   spp <- crop(raster(list_mdp[x]), malla_wgs, datatype='INT1S') %>% mask(., malla_wgs, datatype='INT1S') %>% projectRaster(., malla_lcc, method="ngb")
#   spp <- Tapir_azar + spp; dataType(spp) <- 'INT1S'; return(spp)
# })
# saveRDS(result_overlap, file = "./output/result_overlap_random_vs_spp.Rds")
result_overlap <- readRDS("./output/result_overlap_random_vs_spp.Rds")

##No especies con intersect azar ----
base_area <- (res(result_overlap_random[[1]])[1]^2)* 1e-6#km2
pb <- txtProgressBar(0, length(result_overlap_random), style = 3)
Interseccion_M2 <- map_dfr(1:length(result_overlap_random), function(x){
  setTxtProgressBar(pb, x)
  x.1 <- result_overlap_random[[x]]
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
hist(Interseccion_M2$`% overlap`)
mean(Interseccion_M2$`% overlap`[which(Interseccion_M2$`% overlap` > 0)])
length(Interseccion_M2$`% overlap`[which(Interseccion_M2$`% overlap` >= 10)])
length(which(Interseccion_M2$intersect == "Yes"))
boxplot(Interseccion_M2$`% overlap`[which(Interseccion_M2$`% overlap`>0)])
saveRDS(Interseccion_M2, file = "./output/table_overlap_random_vs_spp.Rds")

## Riqueza azar ----
pb <- txtProgressBar(0, length(result_overlap_random), style = 3)
result_random_S <- lapply(which(Interseccion_M2$`% overlap` >= 10), function(x){
  setTxtProgressBar(pb, x)
  if(2 %in% unique(result_overlap_random[[x]][])){
    x<- result_overlap_random[[x]]; x[x != 2] <- NA; x[x == 2] <- 1
  } else {
    x <- NULL
  }
  return(x)
})
result_random_S <- compact(result_random_S); result_random_S$fun <- sum; result_random_S <- do.call(mosaic, result_random_S)
result_random_S[result_random_S == 0] <- NA; dataType(result_random_S) <- 'INT2U'
writeRaster(result_random_S, "./output/Species_Random_patches.Tif", overwrite = TRUE, datatype='INT2U')

#Grafica de especies en IUCN ----
IUCN <- read.csv("./input/Generalidades_estadoconservacion.csv")
Interseccion_M2$IUCN <- map_chr(Interseccion_M2$especie, function(x){
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
Interseccion_M2$`intersect_10%`<- "No"
Interseccion_M2$`intersect_10%`[which(Interseccion_M2$`% overlap` >= 11)]<- "Yes"
Interseccion_M2_IUCN <- Interseccion_M2[which(Interseccion_M2$`intersect_10%` == "Yes"),]
Interseccion_M2_IUCN <- merge(Interseccion_M2_IUCN, IUCN, by.x = "especie", by.y = "Scientific.name")
names(Interseccion_M2_IUCN)
Interseccion_M2_IUCN$especie
Interseccion_M2_IUCN$IUCN
Interseccion_M2_IUCN$Population.IUCN.SEMARNAT.NatureServe

#Accumulation curve----
Tapir_shp_azar <- Tapir_azar; Tapir_shp_azar[Tapir_shp_azar == 0] <- NA
Tapir_shp_azar <- terra::as.polygons(terra::rast(Tapir_shp_azar), dissolve=FALSE) %>% as(.,"Spatial") %>% st_as_sf()
Tapir_shp_azar <- Tapir_shp_azar[1:1992,]
Tapir_shp_azar$id <- 1:nrow(Tapir_shp_azar)

id_azar2 <- sample(factor(rep(1:100, length.out= nrow(Tapir_shp_azar))))

spaccum_Azar <- sp_accum(id_azar_n = id_azar2, 
                         Tapir_shp = Tapir_shp_azar, 
                         models_list = st_as_stars(stack(result_overlap_random)), 
                         names_list = Interseccion_M2$especie)
saveRDS(spaccum_Azar, file = "./output/spaccum_azar_1.rds")





