library(dplyr)
library(ggplot2)
library(ggpubr)

# regression plot function ----

tl_bhdf_regression<-function(x){
  ggplot(x, aes(x = mean_BHDF, y = mean_TL))+
    geom_point()+
    geom_smooth(method="lm", formula = y ~ x, se=FALSE)+
    stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
             r.accuracy = 0.01,
             label.x = 0.5, label.y = 2.9, size = 4) +
    stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                          label.x = 0.5, label.y = 3, size = 4)
}


# trip and day ----

trip <- c('2022_07', '2022_10', '2023_01')

lapply(trip, function(x){

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

names(whalength_merge)
names(ID_alt_merge)

# merge measurements with ID data ----

length_ID_merge<-whalength_merge%>%
  filter(is.na(Label..for.mult.whales.per.image.))%>%
  left_join(ID_alt_merge, by = c("Filename" = "vlc_filename"))
names(length_ID_merge)

# add in data where more than one dolphin was measured in one photo ----
length_ID_merge<-whalength_merge%>%
  filter(!is.na(Label..for.mult.whales.per.image.))%>%
  left_join(ID_alt_merge, by = c("Filename" = "vlc_filename", "Label..for.mult.whales.per.image." = "ID"))%>%
  bind_rows(length_ID_merge)%>%
  filter(!is.na(ID))%>%
  filter(ID != '-')%>%
  arrange(Filename)

# checks ----

length_ID_merge%>%
  mutate(length_diff = actual_length-Total.Length..m.)%>%
  filter(abs(length_diff)>0.00005)

# total length only for individuals where I entered TL measurement

tl_trip<-length_ID_merge%>%
  filter(!is.na(actual_length))%>%
  group_by(ID)%>%
  #do widths from these individuals as well
  dplyr::summarise(mean_TL = mean(Total.Length..m.), SD_TL = sd(Total.Length..m.))
  
ggplot(tl_trip)+
  geom_point(aes(x = mean_TL, y = SD_TL))  

tl_trip%>%filter(SD_TL > 0.1)

# rostrum to bh and bh to df for anyone who had it done
bh_trip<-length_ID_merge%>%
  group_by(ID)%>%
  dplyr::summarise(mean_RBH = mean(Rostrum.BH), SD_RBH = sd(Rostrum.BH),
                   mean_BHDF = mean(BH.DF.insertion), SD_BHDF = sd(BH.DF.insertion))%>%
  filter(!is.na(mean_BHDF))

bh_trip%>%filter(SD_RBH > 0.1 | SD_BHDF > 0.1)

bh_tl<-bh_trip%>%
  left_join(tl_trip, by = "ID")%>%
  filter(!is.na(mean_TL))%>%
  mutate(trip = trip)

write.csv(bh_tl, paste0(path,"./bh_tl_",trip,".csv"), row.names = F)

write.csv(length_ID_merge, paste0(path,"./length_ID_merge_",trip,".csv"), row.names = F)

#tl_bhdf_regression(bh_tl)

})


# merge all trip means ----
bh_tl_list<-list.files("./Images", "bh_tl_", ignore.case = TRUE, full.names = T, recursive = T)

bh_tl_read<-lapply(bh_tl_list, function(x) read.csv(x, header = T))

bh_tl_merge <- do.call(rbind, bh_tl_read)

tl_bhdf_regression(bh_tl_merge)

prop<-bh_tl_merge%>%
  mutate(TL_BHDF_prop = mean_TL/mean_BHDF)%>%
  dplyr::summarise(mean_prop = mean(TL_BHDF_prop))
  
est_calc<-bh_tl_merge%>%
  mutate(est_TL = mean_BHDF*prop$mean_prop,
         diff = est_TL-mean_TL)

ggplot(est_calc, aes(x = est_TL, y = diff))+
  geom_point()+
  geom_density_2d()

est_calc%>%filter(diff > 0.2)

# same measurement bh to df compared to total length

length_ID_merge_list<-list.files("./Images", "length_ID_merge_", ignore.case = TRUE, full.names = T, recursive = T)

length_ID_merge_read<-lapply(length_ID_merge_list, function(x) read.csv(x, header = T))

length_ID_merge <- do.call(rbind, length_ID_merge_read)

#only include REAL total lengths
length_ID_filter<-length_ID_merge%>%filter(!is.na(actual_length))


ggplot(length_ID_filter, aes(x = BH.DF.insertion, y = Total.Length..m.))+
  geom_point()+
  geom_smooth(method="lm", formula = y ~ x, se=FALSE)+
  stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
           r.accuracy = 0.01,
           label.x = 0.5, label.y = 2.9, size = 4) +
  stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                        label.x = 0.5, label.y = 3, size = 4)+
  theme_bw()+
  xlab("BH to DF length (m)")+
  ylab("Total length (m)")


# noodle ----

length_ID_merge%>%
  group_by(issue)%>%
  dplyr::summarise(mean_TL = mean(Total.Length..m.), SD_TL = sd(Total.Length..m.), 
                   mean_RBH = mean(Rostrum.BH), SD_RBH = sd(Rostrum.BH),
                   mean_BHF = mean(BH.DF.insertion), SD_BHF = sd(BH.DF.insertion))%>%
  mutate(diff_TL_per = (mean_TL/1.947), ## noodle total length
         diff_RBH_per = (mean_RBH/0.462), ## noodle orange length
         diff_BHF_per = (mean_BHF/1.485))%>% ## noodle purple length
  filter(!is.na(issue))
