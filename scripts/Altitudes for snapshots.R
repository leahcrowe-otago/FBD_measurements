library(dplyr)
library(lubridate)
library(stringr)
library(exifr)

# Date and path ----

date <- '20221010'
path <- paste0("./Images/", date)
offset_s<--(seconds(54)) # how much faster is the metadata time than the GPS
# LiDAR ----

lidar_list <- list.files(path, ".CSV")

lidar_files <-
  lapply(lidar_list, function(x)
    read.csv(paste0(path, '/', x), sep = '\t', skip = 2))
lidar_merge <- do.call(rbind, lidar_files)

lidar_merge <- lidar_merge %>%
  dplyr::rename('gmt_date' = 'X.gmt_date') %>%
  mutate(gmt_datetime = ymd_hms(paste(gmt_date, gmt_time), tz = "GMT")) %>%
  mutate(nz_datetime = ymd_hms(format(gmt_datetime, tz = "Pacific/Auckland"))) %>%
  arrange(nz_datetime)

# .MOV start time ----

mov_list <- paste0(path, '/', list.files(path, ".MOV$"))
mov_metadata <- exifr::read_exif(unlist(mov_list))
mov_metadata <-
  exifr::read_exif(unlist(mov_list),
                   tags = c("FileName", "FileModifyDate", "GPSLatitude", "GPSLongitude", "GPSAltitude", "GPSCoordinates", "Duration"))
mov_metadata <- mov_metadata %>%
  mutate(offset_FileModifyDate = ymd_hms(FileModifyDate) - offset_s)%>%
  mutate(mov_datetime_start_gmt = ymd_hms(offset_FileModifyDate, tz = "GMT") - seconds(Duration)) %>% #start time is the metadata time minus the duration of the video
  mutate(mov_datetime_start_nz = ymd_hms(format(mov_datetime_start_gmt, tz = "Pacific/Auckland")))

# .jpg from drone ----

jpg_list <- paste0(path, '/', list.files(path, ".JPG"))
jpg_metadata <- exifr::read_exif(unlist(jpg_list))
jpg_metadata <-
  exifr::read_exif(unlist(jpg_list),
                   tags = c("FileName", "DateTimeOriginal", "GPSLatitude", "GPSLongitude", "GPSAltitude", "GPSCoordinates"))
jpg_metadata <- jpg_metadata %>%
  mutate(offset_DateTimeOriginal = ymd_hms(DateTimeOriginal) - offset_s)%>%
  mutate(jpg_datetime_start_gmt = ymd_hms(offset_DateTimeOriginal, tz = "GMT")) %>% #
  mutate(jpg_datetime_start_nz = ymd_hms(format(jpg_datetime_start_gmt, tz = "Pacific/Auckland")))

# .png snapshot from VLC time ----

png_df <- data.frame(vlc_filename = unlist(list.files(path, ".png")))
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

write.csv(avg_alt_df, paste0(path,'/altperimage_',date,"-",as.character(Sys.Date()),'.csv'), row.names = F)

# Lidar choice ----

lidar_type = as.list(c('exact','average'))

#######

lapply(lidar_type, function(x){
  
  if (x == 'exact'){
  lidar_alt = avg_alt_df$laser_altitude_cm
} else if (x == 'average'){
  lidar_alt = avg_alt_df$avg_alt_cm
}
  
whalelength<-avg_alt_df%>%
  filter(issue == "N")%>% # no issues
  mutate(Folder = date,
         Content = 'sa',
         Notes = '',
         Best_image = vlc_filename,
         Blank1 = '',
         Blank2 = '',
         time = png_datetime_nz,
         Blank3 = '',
         tilt = tilt_deg,
         Lidar = lidar_alt*0.01)%>% #cm to m for whalelength
  dplyr::select(Folder, Content, Notes, Best_image, Blank1, Blank2,
                time, Blank3, tilt, Lidar, longitude, latitude)

openxlsx::write.xlsx(whalelength, paste0(path,'/whalelength_',date,'_',x,'_',as.character(Sys.Date()),'.xlsx'), rowNames = F, sheetName="Sheet1")
})
