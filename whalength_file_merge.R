library(dplyr)
library(ggplot2)
library(ggpubr)

# regression plot function ----

tl_bhdf_regression<-function(x){
  ggplot(x, aes(x = mean_BHDF, y = mean_TL))+
    geom_point()+
    geom_smooth(method="lm", formula = y ~ x, se=FALSE)+
    geom_smooth(method="glm", formula = y ~ x, se=FALSE)+
    stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
             r.accuracy = 0.01,
             label.x = 0.5, label.y = 2.9, size = 4) +
    stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                          label.x = 0.5, label.y = 3, size = 4)
}


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

# total length only for individuals where I entered TL measurement
tl_trip<-length_ID_merge%>%
  filter(!is.na(actual_length))%>%
  group_by(ID)%>%
  #do widths from these individuals as well
  dplyr::summarise(mean_TL = mean(Total.Length..m., na.rm=TRUE), SD_TL = sd(Total.Length..m., na.rm=TRUE),
                   mean_10W = mean(Width.at.10..TL, na.rm=TRUE), SD_10W = sd(Width.at.10..TL, na.rm=TRUE),
                   mean_20W = mean(Width.at.20..TL, na.rm=TRUE), SD_20W = sd(Width.at.20..TL, na.rm=TRUE),
                   mean_30W = mean(Width.at.30..TL, na.rm=TRUE), SD_30W = sd(Width.at.30..TL, na.rm=TRUE),
                   mean_40W = mean(Width.at.40..TL, na.rm=TRUE), SD_40W = sd(Width.at.40..TL, na.rm=TRUE),
                   mean_50W = mean(Width.at.50..TL, na.rm=TRUE), SD_50W = sd(Width.at.50..TL, na.rm=TRUE),
                   mean_60W = mean(Width.at.60..TL, na.rm=TRUE), SD_60W = sd(Width.at.60..TL, na.rm=TRUE),
                   mean_70W = mean(Width.at.70..TL, na.rm=TRUE), SD_70W = sd(Width.at.70..TL, na.rm=TRUE),
                   mean_80W = mean(Width.at.80..TL, na.rm=TRUE), SD_80W = sd(Width.at.80..TL, na.rm=TRUE),
                   mean_90W = mean(Width.at.90..TL, na.rm=TRUE), SD_90W = sd(Width.at.90..TL, na.rm=TRUE),
                   mean_FW = mean(Fluke.width, na.rm = TRUE), SD_FW = sd(Fluke.width, na.rm = TRUE))
  
ggplot(tl_trip)+
  geom_point(aes(x = mean_TL, y = SD_TL))  

tl_trip%>%filter(SD_TL > 0.1)

# rostrum to bh and bh to df for anyone who had it done
bh_trip<-length_ID_merge%>%
  group_by(ID)%>%
  dplyr::summarise(mean_RBH = mean(Rostrum.BH, na.rm=TRUE), SD_RBH = sd(Rostrum.BH, na.rm=TRUE),
                   mean_BHDF = mean(BH.DF.insertion, na.rm=TRUE), SD_BHDF = sd(BH.DF.insertion, na.rm=TRUE))%>%
  filter(!is.na(mean_BHDF))

bh_trip%>%filter(SD_RBH > 0.1 | SD_BHDF > 0.1)

bh_tl<-bh_trip%>%
  left_join(tl_trip, by = "ID")%>%
  filter(!is.na(mean_TL))%>%
  mutate(trip = trip)

#mean
write.csv(bh_tl, paste0(path,"./bh_tl_",trip,".csv"), row.names = F)

#measurement per photo
write.csv(length_ID_merge, paste0(path,"./length_ID_merge_",trip,".csv"), row.names = F)

})


# merge all mean trip data ----
bh_tl_list<-list.files("./Images", "bh_tl_", ignore.case = TRUE, full.names = T, recursive = T)

bh_tl_read<-lapply(bh_tl_list, function(x) read.csv(x, header = T))

bh_tl_merge <- do.call(rbind, bh_tl_read)

bh_tl_merge_Tt<-bh_tl_merge%>%
  filter(ID != "calibration" & ID != "noodle")

#need to replace ageclass_length_peryear for this to work Jan 2024

# length_mean_ID_demo<-bh_tl_merge_Tt%>%
#   mutate(YEAR =  case_when(
#     as.numeric(stringr::str_sub(trip,6,7)) < 9 ~ as.numeric(stringr::str_sub(trip,1,4)),
#     as.numeric(stringr::str_sub(trip,6,7)) >= 9 ~ as.numeric(stringr::str_sub(trip,1,4))+1,
#   ))%>%
#   left_join(ageclass_length_peryear, by = c("ID", "YEAR"))

# length_mean_ID_demo%>%
#   filter(mean_BHDF > 1 & mean_TL < 3.0)

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

######################
####### old mean approach below


#only include bhdf and total lengths measured at the same time
# length_ID_filter<-length_ID_merge_Tt%>%filter(!is.na(actual_length))
# 
# ggplot(length_ID_filter, aes(x = BH.DF.insertion, y = Total.Length..m.))+
#   geom_point(aes(color = AGECLASS), alpha = 0.5, size = 2)+
#   geom_smooth(method="lm", formula = y ~ x, se=FALSE)+
#   geom_smooth(method="auto", se=T)+
#   stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
#            r.accuracy = 0.01,
#            label.x = 0.5, label.y = 2.9, size = 4) +
#   stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
#                         label.x = 0.5, label.y = 3, size = 4)+
#   theme_bw()+
#   xlab("BH to DF length (m)")+
#   ylab("Total length (m)")+
#   scale_color_viridis_d(na.value="#000000")+
#   theme(legend.position = "bottom")
# 
# length_ID_filter%>%
#   filter(ID == "HIVE")%>%
#   mutate(tl_w40 = Width.at.40..TL/Total.Length..m.)
# 
# ## mean bhdf to mean length
# 
# bh_tl<-ggplot(length_mean_ID_demo, aes(x = mean_BHDF, y = mean_TL))+
#   geom_point(aes(color = AGECLASS), alpha = 0.5, size = 2)+
#   geom_smooth(method="lm", formula = y ~ x, se=FALSE)+
#   stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
#            r.accuracy = 0.01,
#            label.x = 0.5, label.y = 2.9, size = 4) +
#   stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
#                         label.x = 0.5, label.y = 3, size = 4)+
#   theme_bw()+
#   xlab("Mean BH to DF length (m)")+
#   ylab("Mean total length (m)")+
#   scale_color_viridis_d(na.value="#000000")+
#   theme(legend.position = "bottom")
# 
# ggsave("./Figures/bh_tl.png", bh_tl, dpi = 320, height = 100, width = 100, units = 'mm')
# 
# length_mean_ID_demo%>%filter(mean_BHDF > 1 & mean_TL < 2.5)
# 
# length_mean_ID_demo%>%filter(is.na(AGECLASS))
# 
# length_ID_merge_Tt%>%filter(ID == 'SHIVERS')%>%dplyr::select(BH.DF.insertion)%>%as.data.frame()
# 
# # width ----
# 
# #avg width at interval per individual per known ageclass
# width_mean_ID_demo<-length_mean_ID_demo%>%
#   ungroup()%>%
#   group_by(ID, SEX, BIRTH_YEAR, FIRST_YEAR, POD, AGECLASS)%>%
#   dplyr::summarise(across(mean_RBH:SD_90W, ~ mean(.x, na.rm = TRUE)))%>%
#   tidyr::pivot_longer(cols = ends_with("0W"),
#                       names_to = "widths",
#                       values_to = "width_m")%>%
#   filter(!grepl("SD", widths) & !is.na(AGECLASS))%>%
#   mutate(width = as.numeric(str_sub(widths,6,7)))
# 
# # avg width per ageclass per season
# width_mean_ID_demo_trip<-length_mean_ID_demo%>%
#   ungroup()%>%
#   group_by(ID, SEX, BIRTH_YEAR, FIRST_YEAR, POD, AGECLASS, trip)%>%
#   dplyr::summarise(across(mean_RBH:SD_90W, ~ mean(.x, na.rm = TRUE)))%>%
#   tidyr::pivot_longer(cols = ends_with("0W"),
#                       names_to = "widths",
#                       values_to = "width_m")%>%
#   filter(!grepl("SD", widths) & !is.na(AGECLASS))%>%
#   mutate(width = as.numeric(str_sub(widths,6,7)))
# 
# width_ratio_trip<-width_mean_ID_demo_trip%>%
#   ungroup()%>%
#   group_by(trip, AGECLASS, SEX, width)%>%
#   mutate(ratio = mean(width_m, na.rm = TRUE)/mean(mean_TL, na.rm = TRUE))
# 
# width_mean_ID_demo_trip%>%filter(SEX != 'X' & AGECLASS == 'A')%>%distinct(ID)
# 
# 
# width_sex<-ggplot()+
#   geom_point(width_mean_ID_demo_trip%>%filter(SEX != 'X' & AGECLASS == 'A'), mapping = aes(x = width, y = width_m/mean_TL, color = trip), alpha = 0.5)+
#   #mean across ageclass
#   geom_line(width_ratio_trip%>%filter(SEX != 'X' & AGECLASS == 'A'), mapping = aes(x = width, y = ratio, color = trip, linetype = SEX), size = 0.75)+
#   #geom_smooth(aes(color = AGECLASS))+
#   theme_bw()+
#   xlab("Total length intervals (%)")+
#   ylab("Width (m)")+
#   scale_color_viridis_d(na.value="#000000")+
#   theme(legend.position = "bottom")+
#   #facet_wrap(~AGECLASS*SEX)+
#   scale_x_continuous(breaks=seq(0, 100, 10))+
#   theme(panel.grid.minor = element_blank())
# 
# width_ratio<-width_mean_ID_demo%>%
#   ungroup()%>%
#   group_by(AGECLASS, width)%>%
#   mutate(ratio = mean(width_m, na.rm = TRUE)/mean(mean_TL, na.rm = TRUE))
# 
# width_ageclass<-ggplot()+
#   #mean individual points
#   geom_point(width_mean_ID_demo, mapping = aes(color = AGECLASS, x = width, y = width_m/mean_TL), alpha = 0.5)+
#   #mean across ageclass
#   geom_line(width_ratio, mapping = aes(x = width, y = ratio, color = AGECLASS), size = 0.75)+
#   #geom_smooth(aes(color = AGECLASS), se=F, method = 'loess')+
#   theme_bw()+
#   xlab("Width at length intervals (%)")+
#   ylab("Width-length ratio (m)")+
#   scale_color_viridis_d(na.value="#000000")+
#   theme(legend.position = "bottom")+
#   scale_x_continuous(breaks=seq(0, 100, 10))+
#   theme(panel.grid.minor = element_blank())
# 
# ggsave("./Figures/width_ageclass.png", width_ageclass, dpi = 320)
# # dev.off()
# WIDTH<-width_mean_ID_demo%>%
#   filter(AGECLASS == "J")%>%
#   mutate(y = width_m/mean_TL)
# 
# width_mean_ID_demo%>%filter(widths == "mean_40W" & AGECLASS == "J" & width_m/mean_TL >0.2)
# width_ageclass+
#   facet_wrap(~AGECLASS)
# #the below will need reworked to investigate individual level things
# 
# # ggplot(length_mean_ID_demo%>%filter(grepl("mean", widths) & 
# #                                   (ID == 'FRENZY' | ID == 'COMET' | ID == 'C3P0' | ID == 'ALYSA' | ID == 'GERBIL' | ID == 'TR120' | ID == 'TICK-TACK')),
# #        aes(x = width, y = width_m))+
# #   geom_point(aes(color = trip), alpha = 0.5)+
# #   geom_smooth(aes(color = trip))+
# #   theme_bw()+
# #   xlab("Total length intervals (m)")+
# #   ylab("Width (m)")+
# #   scale_color_viridis_d(na.value="#000000")+
# #   theme(legend.position = "bottom")+
# #   facet_wrap(~ID)
# # 
# # 
# length_mean_ID_demo%>%filter(SEX == 'F' & AGECLASS == 'S-A')%>%as.data.frame()
# length_mean_ID_demo%>%filter(ID == "SESAME" & YEAR == 2023)%>%as.data.frame()
# 
# female_width<-width_mean_ID_demo_trip%>%filter(SEX == "F")%>%group_by(ID)%>%mutate(distinct_trip = n_distinct(trip))#%>%filter(distinct_trip > 1)
# 
# ggplot(female_width)+
#   #mean individual points
#   geom_point(mapping = aes(color = trip, x = width, y = width_m), alpha = 0.5)+
#   geom_line(mapping = aes(x = width, y = width_m, color = trip), size = 0.75)+
#   facet_wrap(~ID)
# 
# unk_width<-width_mean_ID_demo_trip%>%filter(SEX == "X")%>%group_by(ID)%>%mutate(distinct_trip = n_distinct(trip))%>%filter(distinct_trip > 0)
# 
# ggplot(unk_width)+
#   #mean individual points
#   geom_point(mapping = aes(color = trip, x = width, y = width_m), alpha = 0.5)+
#   geom_line(mapping = aes(x = width, y = width_m, color = trip), size = 0.75)+
#   facet_wrap(~ID)
# 
# male_width<-width_mean_ID_demo_trip%>%filter(SEX == "M")%>%group_by(ID)%>%mutate(distinct_trip = n_distinct(trip))%>%filter(distinct_trip > 0)
# 
# ggplot(male_width)+
#   #mean individual points
#   geom_point(mapping = aes(color = trip, x = width, y = width_m), alpha = 0.5)+
#   geom_line(mapping = aes(x = width, y = width_m, color = trip), size = 0.75)+
#   facet_wrap(~ID)
# 
# ## fluke
# 
# ggplot()+
#   #mean individual points
#   geom_point(length_mean_ID_demo%>%filter(!is.na(AGECLASS) & !is.na(BIRTH_YEAR))%>%mutate(age = YEAR-BIRTH_YEAR), mapping = aes(color = AGECLASS, x = age, y = mean_FW), alpha = 0.5)+
#   #mean across ageclass
#   #geom_line(width_ratio, mapping = aes(x = width, y = ratio, color = AGECLASS), size = 0.75)+
#   #geom_smooth()+
#   theme_bw()+
#   #xlab("Width at length intervals (%)")+
#   #ylab("Width-length ratio (m)")+
#   scale_color_viridis_d(na.value="#000000")+
#   theme(legend.position = "bottom")
#   #scale_x_continuous(breaks=seq(0, 100, 10))+
#   #theme(panel.grid.minor = element_blank())
# 
# 
# length_mean_ID_demo%>%
#   filter(mean_FW > 0.9)
# 
# noodle ----

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

ggplot2::ggsave(paste0("./Figures/noodle.png"), device = "png", dpi = 700, height = 150, width = 300, units = 'mm')


noodlestats<-noodle%>%
  group_by(trip)%>%
  dplyr::summarise(mean_TL = mean(Total.Length..m.),
                   sd_TL = sd(Total.Length..m.),
                   mean_RBH = mean(Rostrum.BH),
                   sd_RBH = sd(Rostrum.BH),
                   mean_BHDF = mean(BH.DF.insertion),
                   sd_BHDF = sd(BH.DF.insertion)
                   )%>%
  mutate(CV_TL = sd_TL/mean_TL,
         CV_RBH = sd_RBH/mean_RBH,
         CV_BHDF = sd_BHDF/mean_BHDF)

noodle_plot<-ggplot(noodle)+
  geom_point(mapping = aes(y = Total.Length..m., x = 1.947, color = "Total length"))+
  geom_point(mapping = aes(y = Rostrum.BH, x = 0.462, color = "Orange"))+
  geom_point(mapping = aes(y = BH.DF.insertion, x = 1.485, color = "Purple"))+
  theme_bw()+
  xlab("Actual length (m)")+
  ylab("Observed length (m)")+
  facet_wrap(~trip)+
  theme(legend.position = c(0.8,0.15))+
  guides(color=guide_legend(title='Segment measured'))

ggplot2::ggsave(paste0("./Figures/noodle.png"), noodle_plot, device = "png", dpi = 700, height = 300, width = 200, units = 'mm')
