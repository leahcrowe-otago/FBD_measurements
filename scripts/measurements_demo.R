library(dplyr)
library(lubridate)
library(stringr)
library(readxl)

# Date and path ----
calf_year <- '2022'
trip <- '2022_07'
#date <- '20220711'

#grabs all alts from altperimage files
path <- paste0("./Images/", trip)

measurement_list <- list.files(path, glob2rx("altperimage*.xlsx"), ignore.case = TRUE, recursive = T)

measurement_files <-
  lapply(measurement_list, function(x)
    read_excel(paste0(path, '/', x)))

measurement_files_filter<-lapply(measurement_files, function(x)
  x%>%
    filter(issue == 'N' & !is.na(actual_length))%>%
    distinct(DATE, vlc_filename, ID, actual_length) #avg_length_mean, avg_length_sd)
  
  )

measurement_merge <- do.call(rbind, measurement_files_filter)

measurement_ID<-measurement_merge%>%
  ungroup()%>%
  dplyr::select(ID, actual_length)%>%
  filter(ID != 'NA')%>%
  group_by(ID)%>%
  mutate(actual_length_mean = mean(actual_length), #mean of lengths measured using altitude at exact time
         actual_length_sd = sd(actual_length), #sd of lengths measured using altitude at exact time
         #avg_length_mean = mean(avg_length), #mean of lengths measured using averaged altitude +/- 2 sec from exact time
         #avg_length_sd = sd(avg_length), #sd of lengths measured using averaged altitude +/- 2 sec from exact time
         n = n(),
         trip = trip)%>%
  distinct(ID, actual_length_mean, actual_length_sd, n, trip)
  
library(ggplot2)
ggplot(measurement_ID)+
  geom_point(aes(x = actual_length_mean, y = actual_length_sd))
          
measurement_ID%>%
  filter(actual_length_sd > 0.1)%>%
  arrange(actual_length_sd)

measurement_merge%>%
  filter(ID == 'CHAZ')

ggplot(measurement_ID)+
  geom_histogram(aes(x = n), binwidth = 1, color = "black")+
  facet_wrap(~trip)

# life history update ----

source('C:/Users/leahm/Documents/git_otago/Fiordland_reporting/scripts/life_history_ageclass update.R', local = TRUE, verbose = F)$value
lifehist<-dbReadTable(con, "life_history_ageclass")

# merge lh with measurements ----

measurements_demo_ind<-measurement_ID%>%
  left_join(lifehist, by = c('ID' = 'NAME'))%>%
  dplyr::select(trip, ID, actual_length_mean, actual_length_sd, n, SEX, BIRTH_YEAR, FIRST_YEAR, ends_with(calf_year))%>%
  dplyr::rename(AGECLASS = ends_with(calf_year))%>%
  mutate(YEAR = calf_year)%>%
  mutate(AGECLASS = case_when(
    ID == "SUNSHINE" ~ "C",
    TRUE ~ AGECLASS
  ))

measurements_demo_ind$AGECLASS<-factor(measurements_demo_ind$AGECLASS, levels = c("C","J","S-A","A"))

# only age grouping ----

age_length<-measurements_demo_ind

age_length$BIRTH_YEAR<-as.numeric(age_length$BIRTH_YEAR)
age_length$FIRST_YEAR<-as.numeric(age_length$FIRST_YEAR)
age_length$YEAR<-as.numeric(age_length$YEAR)

age_calc<-age_length%>%
  ungroup()%>%
  mutate(age = case_when(
    !is.na(BIRTH_YEAR) ~ YEAR-BIRTH_YEAR,
    is.na(BIRTH_YEAR) ~ YEAR-FIRST_YEAR
  ))%>%
  filter(!is.na(BIRTH_YEAR))
  # mutate(age = case_when(
  #   is.na(BIRTH_YEAR) ~ paste0(age,"+"),
  #   !is.na(BIRTH_YEAR) ~ as.character(age)
  # ))
  
age_calc%>%
  filter(is.na(BIRTH_YEAR))

ggplot(age_calc, aes(x = age, y = log(actual_length_mean)))+
  geom_point()+
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(x, bs = "cs"))

age_length%>%
  group_by(trip, AGECLASS)%>%
  mutate(actual_length_mean_agg = mean(actual_length_mean), #mean of lengths measured using altitude at exact time
         actual_length_sd_agg = sd(actual_length_mean), #sd of lengths measured using altitude at exact time
         actual_length_range = paste0(round(min(actual_length_mean),2),'–',round(max(actual_length_mean), 2)),
         #avg_length_mean_agg = mean(avg_length_mean), #mean of lengths measured using averaged altitude +/- 2 sec from exact time
         #avg_length_sd_agg = sd(avg_length_mean), #sd of lengths measured using averaged altitude +/- 2 sec from exact time
         #avg_length_range = paste0(min(avg_length_mean),'–',max(avg_length_mean)),
         num_ind = n())%>%
  distinct(trip, AGECLASS, actual_length_mean_agg, actual_length_sd_agg, actual_length_range, num_ind)

ggplot(age_length)+
  geom_boxplot(aes(x=AGECLASS, y = actual_length_mean))

age_length%>%
  filter(actual_length_mean < 1.8 & AGECLASS == "C")

# age and sex grouping ----

agesex_length<-measurements_demo_ind%>%
  #mutate(YEAR = year(DATE))%>%
  group_by(YEAR, SEX, AGECLASS)%>%
  mutate(SEX = case_when(
    ID == "LOWSPOT" | ID == "KEEL" | ID == "DUB" ~ 'M',
    TRUE ~ SEX))

agesex_length%>%
  mutate(actual_length_mean_agg = mean(actual_length_mean), #mean of lengths measured using altitude at exact time
         actual_length_sd_agg = sd(actual_length_mean), #sd of lengths measured using altitude at exact time
         actual_length_range = paste0(round(min(actual_length_mean),2),'–',round(max(actual_length_mean), 2)),
         #avg_length_mean_agg = mean(avg_length_mean), #mean of lengths measured using averaged altitude +/- 2 sec from exact time
         #avg_length_sd_agg = sd(avg_length_mean), #sd of lengths measured using averaged altitude +/- 2 sec from exact time
         #avg_length_range = paste0(min(avg_length_mean),'–',max(avg_length_mean)),
         num_ind = n())%>%
  distinct(YEAR, SEX, AGECLASS, actual_length_mean_agg, actual_length_sd_agg, actual_length_range, num_ind)

ggplot(agesex_length)+
  geom_boxplot(aes(x=AGECLASS, y=actual_length_mean))+
  facet_wrap(~SEX)
