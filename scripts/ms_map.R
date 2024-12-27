library(marmap);library(dplyr);library(readxl);library(sf);library(ggrepel)


shapefile_path<-"C:/Users/leahm/OneDrive - University of Otago/Documents/git-otago/Fiordland_reporting/shapefiles"

# depth ----

fiordland_base <- marmap::getNOAA.bathy(lon1 = 165, lon2 = 169,
                        lat1 = -43.5, lat2 = -47.5, resolution = 0.5, keep = TRUE)

fiordland_base_raster<-marmap::as.raster(fiordland_base)

test_spdf <- as(fiordland_base_raster, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)%>%filter(layer < 0)
colnames(test_df) <- c("Depth (m)", "x", "y")

nrow(test_df)
test_df<-test_df%>%
  filter(`Depth (m)` < 10)
nrow(test_df)

bathy<-sf::read_sf(shapefile_path, layer = "niwa-new-zealand-bathymetry-contours-2016")
mapview::mapview(bathy)
alliso50<-subset(bathy, ELEVATION == -50)
alliso50<-as.data.frame(st_coordinates(alliso50))
alliso200<-subset(bathy, ELEVATION == -200)
alliso200<-as.data.frame(st_coordinates(alliso200))

# other shapesfiles ----

NZ_coast<-sf::read_sf(shapefile_path, layer = "nz-coastlines-and-islands-polygons-topo-1500k") #https://data.linz.govt.nz/layer/51560-nz-coastlines-and-islands-polygons-topo-1500k/

NZ_lakes<-sf::read_sf(shapefile_path, layer = "nz-lake-polygons-topo-150k") #https://data.linz.govt.nz/layer/50293-nz-lake-polygons-topo-150k/
big_lakes<-subset(NZ_lakes, !is.na(name_ascii))

protected_areas<-sf::read_sf(shapefile_path, layer = "protected-areas") #https://data.linz.govt.nz/layer/53564-protected-areas/
natpark<-subset(protected_areas, (section == "s.4 - National Park"))
FNP<-subset(natpark, (name == "Fiordland National Park"))
mapview::mapview(FNP)
mpa<-subset(protected_areas, (section == "s.3 - Marine Reserve"))

MA<-sf::read_sf(shapefile_path, layer = "Marine_Protected_Areas_under_the_Marine_Management_Act") # https://catalogue.data.govt.nz/dataset/marine-protected-areas-under-the-marine-management-act/resource/ad9dec70-5163-4d21-ab1c-216d6ff314ac
FMA<-subset(MA, Name == "Fiordland Marine Area")

# fiord labels ----
fiord_labels<-data.frame(label = c("Piopiotahi/Milford Sound","Te H\u101pua/Sutherland Sound",
                                   "H\u101wea/Bligh Sound","Te Houhou/George Sound","Taitetimu/Caswell Sound",
                                   "Taiporoporo/Charles Sound","Hinenui/Nancy Sound","Te Awa-o-T\u16b/Thompson Sound",
                                   "Patea/Doubtful Sound","Te R\u101/Dagg Sound",
                                   "Te Puaitaha/\nBreaksea Sound","Tamatea/Dusky Sound","Taiari/Chalky Inlet",
                                   "Rakituma/Preservation Inlet", "Motup\u14dhue/Bluff", "Rakiura/Stewart Island"),
                         lat = c(-44.55, -44.72,
                                 -44.77, -44.85, -45.01,
                                 -45.05, -45.1, -45.15,
                                 -45.27, -45.39,
                                 -45.59, -45.75, -46.02,
                                 -46.1, -46.597, -47),
                         lon = c(167.8, 167.55,
                                 167.48, 167.35, 167.14,
                                 167.09, 167.03, 166.97,
                                 166.88, 166.78,
                                 166.67, 166.47, 166.51,
                                 166.6, 168.33, 167.75))


# map ----

bathy_FMA<-ggplot()+
  geom_raster(data=test_df, aes(x=x, y=y, fill=`Depth (m)`), alpha=0.8)+
  geom_sf(data = FMA, alpha = 0.05, color = "aquamarine", lwd = 0.4)+
  #Doubtful complex
  geom_polygon(mapping = aes(y=c(-45.13,-45.15,-45.23,-45.28,-45.43,-45.50,-45.44),
                             x=c(166.97,166.96,166.88,166.85,166.91,167.06,167.46)), fill = "gold", alpha = 0.9)+
  #Dusky complex
  geom_polygon(mapping = aes(y=c(-45.55,-45.61,-45.73,-45.83,-45.76,-45.44),
                             x=c(166.68,166.60,166.45,166.45,167.04,167.03)), fill = "purple", alpha = 0.9)+
  geom_sf(data = NZ_coast, alpha = 1, fill = "white", lwd = 0.1)+
  geom_path(alliso200, mapping = aes(X,Y,group = L2), color = "steelblue4", alpha = 0.7, linewidth = 0.1)+
  geom_path(alliso50, mapping = aes(X,Y,group = L2), color = "antiquewhite4", alpha = 0.9, linewidth = 0.1)+
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  xlab("Longitude")+
  ylab("Latitude")+
  geom_sf(data = FNP, fill = "darkgreen", alpha = 0.3, lwd = 0.1)+
  geom_sf(data = big_lakes, alpha = 0.7, fill = "steelblue2", lwd = 0.1)+
  coord_sf(xlim = c(165.2,168.5), ylim = c(-47.2,-44.05), crs = 4269)+
  theme(legend.position = c(0.90, 0.40),
        legend.title =  element_text(size = 6),
        legend.margin = margin(c(1, 1, 1, 1)),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size = 5),
        legend.box.background = element_rect(color = "white",fill = "white"),
        legend.key = element_rect(fill = NA),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+
  geom_text_repel(data = fiord_labels, aes(x = lon, y = lat, label = label), size = 2.3, min.segment.length = 0, force_pull = 2, box.padding = 0.1,
                  nudge_x = c(rep(-0.6,5),-0.7,-0.6,-0.8,-0.8,-0.6,-1,-0.6,-0.6,rep(-0.6,3)),
                  nudge_y = c(0.1,0.15,
                              0.1,0.1,0.15,
                              0.1,0.07,0.01,
                              0,0,
                              0.08,-0.01,-0.1,
                              -0.15,0.02,0))

bathy_FMA

###small NZ with box ----

dunedin<-data.frame(label = c("\u14ctepoti/Dunedin"),
                         lat = c(-45.880),
                         lon = c(170.501))

base<-ggplot()+
  geom_sf(data = NZ_coast, alpha = 0.9, fill = "white", lwd = 0.1)+
  theme(panel.background = element_rect(fill = "lightblue"),
        panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  xlab("Longitude")+
  ylab("Latitude")

NZ<-base+
  coord_sf(crs = 4269)+
  theme(panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  geom_rect(mapping = aes(xmin = 165.2, xmax = 168.7, ymin = -47.5, ymax = -44), fill = NA, color = "red")+
  theme_void()

NZ

#### map a all together ----

map_a<-cowplot::ggdraw() +
  cowplot::draw_plot(bathy_FMA) +
  cowplot::draw_plot(NZ, x = 0.2, y = 0.65, width = 0.2, height = 0.3)

ggsave("./figures/map.png", map_a, dpi = 700, height = 6, width = 4, units = 'in')
ggsave("./figures/map.svg", map_a, dpi = 700, height = 6, width = 4, units = 'in')

# #### map b close up of Dusky/Doubtful ----
# 
# complex_base <- marmap::getNOAA.bathy(lon1 = 166, lon2 = 167.5,
#                         lat1 = -45.9, lat2 = -45.0, resolution = 0.1)
# complex_base_raster<-marmap::as.raster(complex_base)
# complex_spdf <- as(complex_base_raster, "SpatialPixelsDataFrame")
# complex_df <- as.data.frame(complex_spdf)%>%filter(layer <= 0)
# colnames(complex_df) <- c("Depth (m)", "x", "y")
# 
# map_b<-ggplot()+
#   geom_raster(data=test_df, aes(x=x, y=y, fill=`Depth (m)`), alpha=0.8)+
#   geom_sf(data = FMA, alpha = 0.05, color = "aquamarine", lwd = 0.4)+
#   geom_sf(data = NZ_coast, alpha = 0.9, fill = "white", lwd = 0.1)+
#   geom_path(alliso200, mapping = aes(X,Y,group = L2), color = "steelblue4", alpha = 0.7, linewidth = 0.1)+
#   geom_path(alliso50, mapping = aes(X,Y,group = L2), color = "antiquewhite4", alpha = 0.9, linewidth = 0.1)+
#   theme(panel.background = element_rect(fill = "lightblue"),
#         panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "black"),
#         panel.border = element_rect(colour = "black", fill=NA, size=1))+
#   xlab("Longitude")+
#   ylab("Latitude")+
#   geom_sf(data = FNP, fill = "darkgreen", alpha = 0.3)+
#   geom_sf(data = big_lakes, alpha = 0.7, fill = "steelblue2")+
#   coord_sf(xlim = c(166.3,167.2), ylim = c(-45.85,-45.11), crs = 4269)+
#   theme(legend.position = "none",
#         axis.text = element_text(size = 8),
#         axis.title = element_text(size = 8))
# 
# map_b
# 
# ggsave("./figures/map_b.png", map_b, dpi = 700, height = 6, width = 4, units = 'in')
# ggsave("./figures/map_b.svg", map_b, dpi = 700, height = 6, width = 4, units = 'in')

