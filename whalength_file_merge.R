library(dplyr)
library(ggplot2)
library(ggpubr)

# trip and day ----

trip <- c('2022_05','2022_07', '2022_10', '2023_01', '2023_05', '2023_06','2023_07','2023_11')

#grab all individual csv outputs from Whalength and merge & write bh_tl file for each trip 
lapply(trip, function(x){
  #x<-"2023_07"
  trip<-x
  print(x)

# file paths ----

path <- paste0("./Images/", trip)
date_list <- grep(list.files(path), pattern='.xlsx', invert=TRUE, value=TRUE)

sg_path <- lapply(date_list, function(x) paste0(path,"/",x,"/screen grabs"))
mf_path <- lapply(date_list, function(x) paste0(path,"/",x,"/measure_files"))

whalength_list <- as.list(unlist(lapply(sg_path, function(x) list.files(x, ".CSV", ignore.case = TRUE, full.names = T))))

whalength_files <-
  lapply(whalength_list, function(x)
    read.csv(paste0(x), skip = 1)[, 1:50])

whalength_merge <- do.call(rbind, whalength_files)
head(whalength_merge)

ID_alt<-as.list(unlist(lapply(mf_path, function(x) list.files(x, glob2rx("altperimage*.xlsx"), full.names = T))))

ID_alt_read<-lapply(ID_alt, function(x) readxl::read_excel(x, range = cellranger::cell_cols("A:AI"), na = "NA")%>%filter(issue == "N"))
head(ID_alt_read)

ID_alt_merge <- do.call(rbind, ID_alt_read)

#names(ID_alt_read[[2]])
names(whalength_merge)
names(ID_alt_merge)

# merge measurements with ID data ----

length_ID_merge<-whalength_merge%>%
  filter(is.na(Label..for.mult.whales.per.image.))%>%
  left_join(ID_alt_merge, by = c("Filename" = "vlc_filename"))
names(length_ID_merge)

# add in data where more than one dolphin was measured in one photo ----

if(nrow(whalength_merge%>%
   filter(!is.na(Label..for.mult.whales.per.image.))) > 0){

length_ID_merge<-whalength_merge%>%
  filter(!is.na(Label..for.mult.whales.per.image.))%>%
  left_join(ID_alt_merge, by = c("Filename" = "vlc_filename", "Label..for.mult.whales.per.image." = "ID"))%>%
  bind_rows(length_ID_merge)

}

length_ID_merge<-length_ID_merge%>%
  filter(!is.na(ID))%>%
  arrange(Filename)%>%
  mutate(trip = trip)

# checks ----

length_ID_merge%>%
  mutate(length_diff = actual_length-Total.Length..m.)%>%
  filter(abs(length_diff)>0.00005)

length_ID_merge[length_ID_merge == 0] <- NA

#measurement per photo
write.csv(length_ID_merge, paste0(path,"./length_ID_merge_",trip,".csv"), row.names = F)

})


#####
# same measurement bh to df compared to total length

length_ID_merge_list<-list.files("./Images", "length_ID_merge_", ignore.case = TRUE, full.names = T, recursive = T)

length_ID_merge_read<-lapply(length_ID_merge_list, function(x) read.csv(x, header = T))

length_ID_merge <- do.call(rbind, length_ID_merge_read)
head(length_ID_merge)

length_ID_merge_Tt<-length_ID_merge%>%
  filter(ID != "calibration")%>%
  filter(ID != "noodle")%>%
  mutate(YEAR =  case_when(
    as.numeric(stringr::str_sub(trip,6,7)) < 9 ~ as.numeric(stringr::str_sub(trip,1,4)),
    as.numeric(stringr::str_sub(trip,6,7)) >= 9 ~ as.numeric(stringr::str_sub(trip,1,4))+1,
  ))

head(length_ID_merge_Tt)

# Fig. S4, noodle ----

noodle<-length_ID_merge%>%
  filter(ID == "calibration" | ID == "noodle")%>%
  mutate(Total.Length..m. = case_when(
    trip == "2023_07" ~ NA,
    TRUE ~ Total.Length..m.
  ))

noodle_diff<-noodle%>%
  ungroup()%>%
  mutate(diff_TL = Total.Length..m.-1.947, ## noodle total length
         diff_RBH = Rostrum.BH-0.462, ## noodle orange length
         diff_BHDF = BH.DF.insertion-1.485) ## noodle purple length

ggplot(noodle_diff)+
  geom_point(mapping = aes(y = diff_TL, x = 1.947, color = "Total length"))+
  geom_point(mapping = aes(y = diff_RBH, x = 0.462, color = "Orange"))+
  geom_point(mapping = aes(y = diff_BHDF, x = 1.485, color = "Purple"))+
  theme_bw()+
  xlab("Actual length (m)")+
  ylab("Difference between observed and actual length (m)")+
  facet_wrap(~trip, ncol = 4)+
  theme(legend.position = c(0.88,0.2))+
  guides(color=guide_legend(title='Segment measured'))+
  scale_color_manual(values = c("orange","purple","black"))

ggplot2::ggsave(paste0("./Figures/noodle.png"), device = "png", dpi = 700, height = 120, width = 250, units = 'mm')
