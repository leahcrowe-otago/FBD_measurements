library(dplyr)
library(lubridate)
library(stringr)
library(exifr)

# Date and path ----
#year <- '2022'
trip <- '2022_10'
date <- '20221025'
offset_s<-(seconds(10)) # how much faster is the metadata time than the GPS

path <- paste0("./Images/", trip, "/", date)
#path <- paste0("Z:/Fiordland bottlenose dolphin/Long Term Monitoring/Dusky Sound Dolphin Monitoring/", year,"/", trip,"/UAS/",date)
# LiDAR ----

lidar_list <- list.files(path, ".CSV", ignore.case = TRUE)

lidar_files <-
  lapply(lidar_list, function(x)
    read.csv(paste0(path, '/', x), sep = '\t', skip = 2))
lidar_merge <- do.call(rbind, lidar_files)

lidar_merge <- lidar_merge %>%
  dplyr::rename('gmt_date' = 'X.gmt_date') %>%
  mutate(gmt_datetime = ymd_hms(paste(gmt_date, gmt_time), tz = "GMT")) %>%
  mutate(nz_datetime = with_tz(gmt_datetime, tz = "Pacific/Auckland")) 

# .MOV start time ----

mov_list <- paste0(path, '/', list.files(path, ".MOV$"))
mov_metadata <- exifr::read_exif(unlist(mov_list))
mov_metadata <-
  exifr::read_exif(unlist(mov_list),
                   tags = c("FileName", "FileModifyDate", "GPSLatitude", "GPSLongitude", "GPSAltitude", "GPSCoordinates", "Duration"))
mov_metadata <- mov_metadata %>%
  mutate(offset_FileModifyDate =  ymd_hms(FileModifyDate) - offset_s)%>%
  mutate(mov_datetime_start_gmt = round_date(ymd_hms(offset_FileModifyDate, tz = "GMT") - seconds(Duration), unit = "seconds")) %>% #start time is the metadata time minus the duration of the video
  mutate(mov_datetime_start_nz = with_tz(mov_datetime_start_gmt, tz = "Pacific/Auckland"))

# .jpg from drone ----

jpg_list <- paste0(path, '/', list.files(path, ".JPG"))

if (length(list.files(path, ".JPG")) > 0) {
jpg_metadata <- exifr::read_exif(unlist(jpg_list))
jpg_metadata <-
  exifr::read_exif(unlist(jpg_list),
                   tags = c("FileName", "DateTimeOriginal", "GPSLatitude", "GPSLongitude", "GPSAltitude", "GPSCoordinates"))
jpg_metadata <- jpg_metadata %>%
  mutate(offset_DateTimeOriginal = ymd_hms(DateTimeOriginal) - offset_s)%>%
  mutate(jpg_datetime_nz = ymd_hms(offset_DateTimeOriginal, tz = "Pacific/Auckland")) %>% #
  mutate(jpg_datetime_gmt = with_tz(jpg_datetime_nz, tz = "GMT"))
} else {
  jpg_metadata = data.frame(SourceFile = NA, 
                              FileName = NA, 
                              GPSLatitude = NA, 
                              GPSLongitude = NA, 
                              jpg_datetime_nz = NA)
}

x<-jpg_metadata %>% 
  dplyr::select(SourceFile, FileName, GPSLatitude, GPSLongitude, jpg_datetime_nz)%>%
  dplyr::rename(nz_datetime = jpg_datetime_nz)
y<-mov_metadata %>% 
  dplyr::select(SourceFile, FileName, GPSLatitude, GPSLongitude, mov_datetime_start_nz)%>%
  dplyr::rename(nz_datetime = mov_datetime_start_nz)

DJI_lidar<-x%>%
  bind_rows(y)%>%
  left_join(lidar_merge, by = 'nz_datetime')%>%
  mutate(DJI_time_correction = offset_s)%>%
  filter(!is.na(SourceFile))

#lidar_merge%>%filter(nz_datetime == '2022-02-20 12:26:54')

if (file.exists(paste0(path,"/measure_files")) == FALSE) {
  dir.create(paste0(path,"/measure_files"))
}

##distance between points from LiDAR and points from GPS
m_nm<-1/1852

DJI_lidar<-DJI_lidar%>%
  mutate(dist_km = geosphere::distVincentyEllipsoid(matrix(c(GPSLongitude,GPSLatitude), ncol = 2),matrix(c(longitude, latitude), ncol =2),a=6378137, f=1/298.257222101)*m_nm)


write.csv(DJI_lidar, paste0(path,"/measure_files/DJI_lidar.csv"), row.names = F)

# .png snapshot from VLC time ----

png_df <- data.frame(vlc_filename = unlist(list.files(paste0(path,"/screen grabs"), ".png")))
png_df <- png_df %>%
  mutate(
    mov_filename = str_sub(vlc_filename, 1, 12),
    houroffset = as.numeric(str_sub(vlc_filename, 14, 15)),
    minoffset = as.numeric(str_sub(vlc_filename, 17, 18)),
    secoffset = as.numeric(str_sub(vlc_filename, 20, 21))
  )

png_datetime <- png_df %>%
  left_join(mov_metadata, by = c("mov_filename" = "FileName")) %>%
  mutate(png_datetime_gmt = mov_datetime_start_gmt + hours(houroffset) + minutes(minoffset) + seconds(secoffset)) %>%
  mutate(png_datetime_nz = mov_datetime_start_nz + hours(houroffset) + minutes(minoffset) + seconds(secoffset)) %>%
  dplyr::select(vlc_filename, png_datetime_gmt, png_datetime_nz)

# LiDAR data for PNGs ----

# time in seconds over which to average altitudes on either side of exact
# ie a value of 4 here will include 4 sec before and 4 sec after and average those 9 values

time_offset_alt_avg <- 2

png_datetime_split <- split(png_datetime, png_datetime$vlc_filename)

png_datetime_split<-lapply(png_datetime_split, function(x)
  x %>% 
    mutate(avg_start_nz = png_datetime_nz - seconds(time_offset_alt_avg),
           avg_end_nz = png_datetime_nz + seconds(time_offset_alt_avg)))


avg_alt<-lapply(png_datetime_split, function(x) 
{
  #x<-png_datetime_split[[1]]
  print(x)
  exact<-x %>%
    left_join(lidar_merge, by = c('png_datetime_nz' = 'nz_datetime'))
  
  y<-lidar_merge%>%
    filter(nz_datetime >= x$avg_start_nz & nz_datetime <= x$avg_end_nz)%>%
    mutate(num_avg = n())
print(y)
  exact%>%
    mutate(avg_alt_cm = mean(y$laser_altitude_cm),
           time_offset_alt_avg = time_offset_alt_avg,
           num_avg = nrow(y),
           range = max(y$laser_altitude_cm)-min(y$laser_altitude_cm))

  })

avg_alt_df<-do.call(rbind, avg_alt)

avg_alt_df<-avg_alt_df%>%
  dplyr::select(vlc_filename, png_datetime_gmt, png_datetime_nz, laser_altitude_cm, avg_alt_cm, tilt_deg, everything())%>%
  mutate(issue = case_when(
    laser_altitude_cm > 5000 ~ 'Y', #erroneously high alt value
    avg_alt_cm > 5000 ~ 'Y', #erroneously high average alt value
    time_offset_alt_avg == 2 & range > 80 ~ 'Y', #range of >80cm over 5secs
    time_offset_alt_avg == 1 & range > 40 ~ 'Y', #range of >40cm over 3secs
    TRUE ~ "N"
  ))

# READ IDs ----

IDs<-readxl::read_excel(paste0("./Images/",trip,"/",trip, " UAS IDs.xlsx"))

IDs<-IDs%>%
  filter(ymd(DATE) == ymd(date))%>%
  mutate(vlc_filename = case_when(
    !is.na(vlc_filename) ~ paste0(vlc_filename,".png"),
    is.na(vlc_filename) ~ vlc_filename
  ))

avg_alt_ID<-IDs%>%
  left_join(avg_alt_df, by = c("vlc_filename"))%>%
  mutate(laser_altitude_m = laser_altitude_cm*0.01,
         avg_alt_m = avg_alt_cm*0.01,
         actual_length = '',
         avg_length = '')%>%
  dplyr::select(-laser_altitude_cm, -avg_alt_cm)

write.csv(avg_alt_ID, paste0(path,'/measure_files/altperimage_',date,"-",as.character(Sys.Date()),'.csv'), row.names = F)

# Lidar choice ----

lidar_type = as.list(c('exact','average'))

# output for potential feed into whalength ----

lapply(lidar_type, function(x){
  
  if (x == 'exact'){
  avg_alt_ID<-avg_alt_ID%>%
    mutate(lidar_alt = laser_altitude_m)
} else if (x == 'average'){
  avg_alt_ID<-avg_alt_ID%>%
    mutate(lidar_alt = avg_alt_m)
}
  
whalelength<-avg_alt_ID%>%
  mutate(Folder = date,
         Content = 'sa',
         Notes = '',
         Best_image = vlc_filename,
         Blank1 = '',
         Blank2 = '',
         time = png_datetime_nz,
         Blank3 = '',
         tilt = tilt_deg, 
         Lidar = lidar_alt)%>% #cm to m for whalelength
  filter(issue == "N")%>% # no issues
  dplyr::select(ID,Folder, Content, Notes, Best_image, Blank1, Blank2,
                time, Blank3, tilt, Lidar, longitude, latitude) 
#movie_time is nonsense because it doesn't recognize the format of mm:ss as character to identify when in the vid the dolphin shows up
openxlsx::write.xlsx(whalelength, paste0(path,'/measure_files/whalelength_',date,'_',x,'_',as.character(Sys.Date()),'.xlsx'), rowNames = F, sheetName="Sheet1")
})

