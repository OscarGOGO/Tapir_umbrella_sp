# Funciones ----
rand_vect <- function(N, M, sd = 1, pos.only = TRUE, seed = NULL) {
  if(!is.null(seed)){
    set.seed(seed)  
  }
  
  vec <- rnorm(N, M/N, sd)
  if(abs(sum(vec)) < 0.01){vec <- vec + 1} 
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  
  for(. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  
  if(pos.only){
    while (any(vec < 0)) {
      negs <- vec < 0
      pos  <- vec > 0
      vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 2
      vec[pos][i]  <- vec[pos][i <- sample(sum(pos), 1)] - 2
    }
  } 
  
  return(vec)
}

sp_accum <- function(id_azar_n, Tapir_shp, models_list, names_list){
  '%!in%' <- function(x,y)!('%in%'(x,y))
  id_azar_n <- as.numeric(id_azar_n)
  x.1 <- list()
  pb <- txtProgressBar(0, 100, style = 3)
  
  for(i in 1:100){
    setTxtProgressBar(pb, i)
    i.0 <- which(id_azar_n %in% seq(1:i))
    i.1 <- models_list[ms_dissolve(Tapir_shp[i.0,]) %>% st_buffer(0)] %>% as(., "Raster")
    names(i.1)[1:length(names_list)] <- names_list
    i.1 <- as(models_list, "Raster") %>% mask(., as(Tapir_shp[i.0,], "Spatial"))
    i.2 <- map_dbl(1:dim(i.1)[3], function(j){
      return(if(2 %in% unique(i.1[[j]][])){1}else{0})})
    
    if(i > 1){
      i.3a <- names(i.1)[which(i.2 == 1)]  
      i.3a <- i.3a[which(i.3a %!in% i.3)]
      x.1[[i]] <- length(i.3a)
      i.3 <- c(i.3, i.3a)
    } else {
      i.3 <- names(i.1)[which(i.2 ==1)]  
      x.1[[i]] <- length(i.3)
    }
  }
  x.1 <- do.call(rbind, x.1) %>% as.data.frame(); names(x.1) <- "S"
  x.1$Cumsum <- cumsum(x.1$S)
  return(x.1)
}
