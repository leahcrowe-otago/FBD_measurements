library(readxl);library(dplyr);library(ggplot2)

TtUSPac<-read_excel('./data/2006-2022 Tt (US Pac).xlsx', sheet = "2006-2022 Tt (US Pac)")

str(TtUSPac)

unique(TtUSPac$Species)
unique(TtUSPac$State)
unique(TtUSPac$`Latitude Units`)
unique(TtUSPac$`Confidence Code`)
unique(TtUSPac$`Out of Habitat`)
unique(TtUSPac$`Fishery Interaction`)
unique(TtUSPac$`Age Class`)

TtUSPac%>%filter(State == 'WA')%>%as.data.frame()
TtUSPac%>%filter(`Out of Habitat` == "Y")%>%as.data.frame()
TtUSPac%>%filter(`Fishery Interaction` == "Y")%>%as.data.frame()

TtUSPac_confirmed<-TtUSPac%>%
  filter(!grepl("Unconfirmed", `Confidence Code`) & Genus == "Tursiops" & Species == "truncatus")%>%
  mutate(length_cm = case_when(
    `Length Units` == "in" ~ as.numeric(Length) * 2.54, 
    `Length Units` == "cm" ~ as.numeric(Length),
    is.na(`Length Units`) ~ as.numeric(Length)
    ),
    latitude_dd = case_when(
      `Latitude Units` == "deg/min/sec" ~ as.numeric(substr(Latitude, 1, 2)), #+ (as.numeric(substr(TtUSPac_confirmed$Latitude, 4, 5))/60)),
      `Latitude Units` == "decimal degrees" ~ as.numeric(substr(Latitude, 1, 2))
    )
    )%>%
  #filter(!is.na(latitude_dd))%>%
  filter(`Length actual/estimate` != "notMeasured")
  

TtUSPac_confirmed$latitude_dd

ggplot(TtUSPac_confirmed)+
  geom_point(aes(x = length_cm, y = latitude_dd, color = State))

ggplot(TtUSPac_confirmed)+
  geom_boxplot(aes(x = as.factor(latitude_dd), y = length_cm))+
  facet_wrap(~`Age Class`)

max(TtUSPac_confirmed$latitude_dd)
min(TtUSPac_confirmed$`Observation Date`)
max(TtUSPac_confirmed$`Observation Date`)


TtUSPac_confirmed%>%
  filter(length_cm > 300)%>%
  tally()

nrow(TtUSPac_confirmed)
