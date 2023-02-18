# Plot accumulation species ----
library(ggplot2)
spaccum_Tapir <- readRDS("./output/spaccum_Tapir_1.rds")
spaccum_Azar <- readRDS("./output/spaccum_azar_1.rds")
spaccum_Azar$Group <- "Random selection"
spaccum_Tapir$Group <- "Tapir distribution model"
spaccum_1 <- rbind(spaccum_Tapir, spaccum_Azar)
spaccum_1$id <- c(1:100, 1:100)
spaccum_1$area <- (spaccum_1$id*(4720.418^2))*1e-6
spaccum_1$PorcentajeArea <- 1:100

p0 <- ggplot(data=spaccum_1, aes(x=PorcentajeArea, y=Cumsum, group=Group, col = Group)) +
  geom_line(size = 2)+
  labs(x = "Distribution area (%)",
       y = "Cumulative species richness", size = 10)+
  theme_minimal(base_size = 18)+
  theme(text = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(hjust = 1))+
  scale_color_manual(values = c("#FF7F0E", "#1F77B4"))
p0
ggsave("./output/Figure1a.tif", plot = p0, device = "tiff", width = 15,
       height = 11, compression = "lzw")

# SP vs traslape---
library(ggcharts)
Interseccion_M1 <- readRDS("./output/table_overlap_10Perc_Tapir_vs_spp.Rds")
Interseccion_M2 <- readRDS("./output/table_overlap_random_vs_spp.Rds")
graph_area_intersection <- map_dfr(1:100, function(x){
  x.1 <- rbind(data.frame("Type" = "Baird's tapir distribution",
                          "Pintersect" = x,
                          "NoSpecies" = length(which(Interseccion_M1$`% overlap`>= x)), 
                          check.names = FALSE),
               data.frame("Type" = "Random selection",
                          "Pintersect" = x,
                          "NoSpecies" = length(which(Interseccion_M2$`% overlap`>= x)),
                          check.names = FALSE)
  )
  return(x.1)})
graph_area_intersection$Type <- factor(graph_area_intersection$Type, 
                                       levels = c("Baird's tapir distribution", "Random selection"))
graph_area_intersection$Pintersect <- factor(graph_area_intersection$Pintersect, levels = 1:100)
graph_area_intersection <- graph_area_intersection[c(which(graph_area_intersection$Type == "Baird's tapir distribution"),
                                                     which(graph_area_intersection$Type == "Random selection")),]
row.names(graph_area_intersection) <- NULL

p1<-pyramid_chart(graph_area_intersection[c(1:30, 101:130),], Pintersect, NoSpecies, Type,
                  xlab = "Number of species")
p1
ggsave("./output/Figure1b.tif", plot = p1, device = "tiff", width = 15,
       height = 11, compression = "lzw")
