#Packages----
library(raster)
library(terra)
library(sf)
library(magrittr)
library(purrr)
library(stars)
library(rmapshaper)

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


## Potential species richness----
pb <- txtProgressBar(0, length(result_overlap), style = 3)
result_S <- lapply(1:length(result_overlap) , function(x){
  setTxtProgressBar(pb, x)
  if(2 %in% unique(result_overlap[[x]][])){
    x<- result_overlap[[x]]; x[x != 2] <- NA  
  } else {
    x <- NULL
  }
  return(x)
})
result_S <- compact(result_S); result_S$fun <- sum; result_S <- do.call(mosaic, result_S)
result_S[result_S == 0] <- NA; dataType(result_S) <- 'INT2U'
writeRaster(result_S, "./output/Species_TapirDist.Tif", overwrite = TRUE, datatype='INT2U')
plot(trim(result_S))


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
base_area <- (res(result_overlap[[1]])[1]^2)* 1e-6#km2
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
Interseccion_M1$`% overlap`
mean(Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap` > 0)])
length(Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap` > 5)])
length(which(Interseccion_M1$intersect == "Yes"))
boxplot(Interseccion_M1$`% overlap`[which(Interseccion_M1$`% overlap`>0)])
saveRDS(Interseccion_M1, file = "./output/table_overlap_Tapir_vs_spp.Rds")

