---
title: "Size estimates of bottlenose dolphins at the southern extreme of their range"
csl: aquatic-biology.csl
bibliography: references.bib
affiliations:
- name: Department of Marine Science, University of Otago, Dunedin, New Zealand
  index: 1
authors:
- name: Leah M. Crowe
  orcid: 0000-0001-9133-8522
  affiliation: 1
- name: Steve M. Dawson
  affiliation: 1
- name: William Rayment
  affiliation: 1
output:
  pdf_document: default
  html_document:
    df_print: paged
  word_document: default
---

Don't like this title because it doesn't encompass how they relate to Scottish dolphins and other Pacific

# ABSTRACT

The critically endangered bottlenose dolphins (*Tursiops truncatus*) inhabiting the fiords of southwestern Aotearoa-New Zealand are living at the southern edge of their range. Previous photogrammetric estimates of Fiordland bottlenose dolphins have included lasermetric estimates of dorsal fins and vessel-mounted systems which restrict length estimation to boat-positive individuals, often juveniles and sub-adults. Modern uncrewed aerial systems (i.e. drones) are relatively inexpensive, and provide more flexibility for collecting photogrammetric data from wildlife to aid in population monitoring. In this study, an off-the-shelf quadcopter equipped with a custom-built LIDAR datalogger was used to collect footage of individuals within the small bottlenose dolphin sub-populations in the Patea-Doubtful (approximately 60 individuals) and Tamatea-Dusky (approximately 120 individuals) Sound complexes. Photogrammetric measurements of individual total lengths were estimated using the software, *whalength*, and longitudinal life history data were applied to determine the range of total lengths for different age and sex-classes. Preliminary results indicate that the average total length for the different age-classes are as follows: adults were 2.88m (SD = 0.14, range = 2.68--3.12, n = 23), sub-adults (3--8 y) were 2.60m (SD = 0.13, n = 12), juveniles (1--2 y) were 2.07m (SD = 0.20, range = 1.86--2.26, n = 3), and calves of the year were 1.84m in July 2022 (SD = 0.14, range = 1.68--1.94, n = 3). Length measurements reported here are similar to estimates for bottlenose dolphins living at the northern extreme of their range, are similar in length to stranding records collected in the United States in the North Pacific and provide a baseline for future comparative studies. Our results further demonstrate the utility of incorporating photogrammetric methods into this monitoring programme to assess population health and resilience of Fiordland bottlenose dolphins.

# INTRODUCTION

## 1.1 Size as health indicator

Morphometric data provides insight into the well-being of wild populations, but can be challenging or impossible to collect directly. For marine species, data is often obtained through dead specimens, capture of individuals, or remote methods. The availability of commercial uncrewed aerial systems (UASs) has made collection of morphometric data of free-swimming marine species more accessible, and has revolutionized the field in the last decade. UASs derived photogrammetric measurement estimates have less error applied to larger objects <!--# back this up -->, and have widely been applied to studying large whale species (e.g. @christiansen2020; @durban2021); however, in recent years, this has become a more common technique for remotely studying the size of smaller marine animal species, including manta rays [@setyawan2022], pinnipeds [@krause2017], manatees [@ramos2022], and dolphins [@cheney2022; @christie2021; @deoliveira2023].

Body condition indices provide valuable insight into the health of populations, particularly relative to the health and resilience of endangered and threatened species. Previous studies have harnessed UAS data in comparative studies to investigate differences in body condition relative to seasonal and reproductive cycles [@bierlich2022; @christiansen2016], prey availability [@durban2021; @christiansen2021], and anthropogenic impacts [@christiansen2020; @stewart2021]. These studies can better inform management strategies by revealing processes that constrain population growth and recovery.

## 1.2 Bottlenose dolphins at the extremes and in between

Common bottlenose dolphins (*Tursiops truncatus)* exist in two predominant ecotypes (coastal and pelagic) inhabiting the worlds' oceans from temperate to tropical waters [@tezanos-pinto2008]. Previous studies have suggested that the larger morphs of this species exist at the extreme northern and southern extents of their range: the population in Scotland has measured larger than both captured individuals in Florida [@read1993] and stranded bottlenose found on the coast of the mid-Atlantic of the United States [@cheney2017], while bottlenose in Aotearoa-New Zealand have previously been estimated to be larger than stranded dolphins off the coast of Texas [@chong2001two]. Until recently, measurements of the bottlenose dolphin populations living at the extremes have been conducted through vessel-mounted stereo-photogrammetry, laser-metric methods, or from stranded specimens. UAS methodology has provided the ability to remotely measure many members of these populations. In this study, we employ this technology to better quantify size estimates of bottlenose dolphins at the southern extreme of their range, and compare them to their northern counterparts with data reported in @cheney2022. In addition, this study includes length measurement data of bottlenose dolphins stranded in the United States and its territories in the Pacific Ocean to further compare length variation by latitude for this species [@swfsc2023].

## 1.3 Fiordland bottlenose dolphins and previous attempts to quantify size

Fiordland bottlenose dolphins are comprised of four sub-populations living at the southern extreme of their range within the fourteen fiords of the Te Moana o Atawhenua-Fiordland Marine Area. The northern pod ranges between Piopiotahi-Milford Sound and Taiporoporo-Charles Sound [@lusseau2005, Crowe et al. in prep], the Doubtful pod ranges between Hinenui-Nancy and Te Ra-Dagg Sound with a few observed forays into Te Puaitaha-Breaksea Sound, the Dusky pod ranges between Te Puaitaha-Breaksea Sound and Rakituma-Preservation Inlet (Crowe et al. in prep), and the southern pod range between Taiari-Chalky Inlet and Rakiura-Stewart Island [@brough2015, Crowe et al. in prep] (Figure \##). Because of their small population sizes (2021 estimates: Doubtful pod = 6# +/- #, Dusky pod = 12# +/- #, @Crowe2022) , as well as their patterns of residency and geographic isolation, this population is collectively considered under the IUCN as 'Critically Endangered' [@currey2009; @FBDIUCN2011].

The first attempts at remotely measuring Fiordland bottlenose dolphins occurred in 1997 when the photo-identification study in the Patea-Doubtful Sound complex was in its seventh year [@chong2001two]. A stereo-photogrammetry system was mounted to the bow of a sailing yacht to obtain measurements of the total length, flippers, and flukes of individuals in the Doubtful pod, and results indicated that this population was larger than stranded bottlenose dolphins found off the coast of Texas, USA [@chong2001two]. The asymptotic length of the Doubtful pod was 3.20m, approximately 0.70m larger than the Texas dolphins; however, the growth curve was fit to measured lengths in order to estimate age, and was fixed between two points representing the estimated age of two stranded dolphins per teeth growth layers [@chong2001two].

Stereo-photogrammetry methods were again attempted in the 2010s to estimate total length as well as explore allometric relationship for Fiordland bottlenose dolphin as well as bottlenose across Aotearoa-New Zealand [@brough2013using]. The infrastructure for collecting these data proved challenging, but findings from seven Doubtful individuals measured, it was suggested that the length between rostrum and dorsal fin insertion could be used as a proxy for total length [@brough2013using], which is similar to other work that has investigated the allometric relationship between the blowhole to dorsal fin insertion length and total length (e.g. @cheney2017); in both comparisons, the proxy lengths are more common to capture above the surface and are more reliably straight based on the movement mechanics of dolphins (REFERENCE?). The measured dolphins from the Doubtful pod were similar in total length to the total length of stranded individuals across the country [@brough2013using].

In addition to efforts to measure total length, studies in the mid-2000s were conducted to investigate the presence of sexual dimorphism in Fiordland bottlenose dolphins using laser-metric methods [@rowe2009; @rowe2008]. These findings suggested dorsal fin size could be used as an indicator of sex [@rowe2009]; however, several individuals from the Doubtful pod were not assigned the correct sex in the life history data used at the time (unpublished data). This study's sex assignment method was additionally applied to the dolphins sighted in the Tamatea-Dusky complex in 2007/2008, the first few years of the monitoring study [@currey2008], and, similarly, many assignments have later proven incorrect as the long-term studies have unfolded (unpublished data).

Despite the previous efforts to measure Fiordland bottlenose dolphins, the methods used have only measured a small fraction of the Doubtful pod, and it is unclear to what extent there is sexual dimorphism in size for this population. In this study, an unoccupied aerial system was used to estimate length and width measurements of Fiordland bottlenose dolphins from overhead imagery. Flights were conducted at several sampling occasions over two years to capture growth (length) and seasonal changes in body condition (width). UASs allowed for greater coverage of the population across age and sex classes as previous vessel-mounted systems have been constrained to boat-positive individuals, which have primarily been juveniles (e.g. @chong2001two). Here we demonstrate that bottlenose dolphins living at their extreme southern range show similar growth and size trends as their northern counterparts, and suggest that there may not be drastic size differences of bottlenose dolphins across the Pacific Ocean.

# METHODS

## 2.1 Data collection

Overhead video footage of bottlenose dolphins was collected during regular survey efforts in the Patea-Doubtful and Tamatea-Dusky complex between February 2022 and November 2023. The Fiordland Marine Area is characterized by heavy rainfall (average ##m of rain each year), therefore, UAS flights were constrained to times when there was no rain, no impact from overhanging vegetation, and low Beaufort (\<= 2) during dolphin sightings. Surveys were conducted from a 5-m vessel powered by a 70-hp, 4-stroke outboard engine, and the UAS was launched and recovered by hand from the stern.

All UAS surveys were conducted using a DJI Inspire 1 Pro equipped with an Olympus 25mm f1.8 lens camera and a custom-added laser altimeter system consisting of a Light Detection and Ranging (LiDAR) and Global Positioning System (GPS) data logger set at a 1-sec sampling rate (see @dawson2017 for details on the configuration). LiDAR derived altitude measurements have previously been shown greatly reduce variability in body measurement calculations compared to barometric and GPS altimeter readings [@dawson2017; @ramos2022; @bierlich2021]. Flight altitude was targeted at 10--15m to comply with our permit issued by the New Zealand Department of Conservation (#87586-MAR), to obtain footage with a high enough resolution to identify individuals, and to also limit behavioural reactions to the platform by the animals. When over the animals, the camera was in a perdendicular position, and video footage collected was shot in 4K resolution (2160 x 3840 pixels) at 25 frames per second; the animals move too quickly for single imagery to be practical.

## 2.2 Video analysis

Video footage was reviewed using VLC player version 3.0.12 [@solutions2006vlc], and for application in total length and width measurements, snapshots of frames were taken as .png files when an individual was in the most elongated position at the surface. This position was most often characterized by the frame where there was a clear view of the rostrum and tail notch slightly subsurface and the blowhole just above the surface (FIGURE reference). In addition, snapshots were taken in instances where the blowhole and dorsal fin insertion were clearly seen for application in exploring the allometric relationship between this length and total length as has been explored in previous studies (e.g. @cheney2017 ). The time of the screen shot image was calculated by subtracting the time in the video from the timestamp of the video, and the altitude of the image was obtained by syncing this image timestamp to the feed from the LiDAR/GPS data logger.

Dolphins captured in image screen shots were scrutinized to identify the individual and associated life history information. Long-term photo-identification studies have occurred almost continuously since 1990 on the Patea-Doubtful Sound population and since 2006 on the Tamatea-Dusky population. Previously, individuals were primarily identified by the nicks, scratches, and scars on the dorsal fin [@currey2008;@j.c.currey2007]; for this study, the catalogue has been expanded to include oblique images of the entire dorsal body (e.g. @cheney2022 ). Demographic classification in terms of age and sex classes are possible for these two sub-population as new entrants have almost exclusively occurred through births to known individuals, and individuals have been sexed through visual observation of the genital slit or a close, consistent association with a calf for females.

## 2.3 Photogrammetry measurements

Field calibrations were done at least once per trip to ensure the systems were working as expected. We used a boat-hook wrapped in pool noodles as our known-length object as it was rigid with a low-profile, and floated just above the surface. Two different lengthed, contrasting coloured noodles on one fully extended boat-hook allowed for measurements of three varying lengths (total length and two shorter lengths of each of the coloured sections to represent smaller individuals) without adding too much cumbersome gear to the project. CALIBRATION RESULTS.

Altitude and tilt angle values from the LiDAR system were used within the software, *whalength*, to conduct photogrammetric measurements [@dawson2017]. An altitude filter was applied to control for errant altimeter readings where the range of values over a 5-second period must be less than 80 cm and also less than 40 cm within a 3-second period [@dickson2021]. Within the *whalength* software, settings for Inspire 1 Pro video stills was used, and the camera to LiDAR offset was set at -1.5cm. Target measurements included: total length (rostrum to fluke notch), rostrum to rostral edge of the blowhole, <!--# probs remove? --> caudal/anterior edge of the blowhole to dorsal fin insertion, width (at 10% increments along the total length), and fluke width<!--# probs remove? -->.

## 2.4 Summary/Statistics

Morphometric measurements were pooled by age and sex class according to the following: adults were at least nine years old or had been sighted over the course of at least eight years if birth year was unknown, sub-adults were less than nine and at least three years old, juveniles were than three and at least one year old, and calves included individuals in their birth year. Birth year was defined by the seasonality of Fiordland bottlenose dolphin calving and begins in the austral spring of the previous calendar year and ends at the end of the winter (e.g. 01 September 2022--31 August 2023); this is similar to definitions used previously [@henderson2014], but encompasses the entire spring and winter season as well as an earlier birth observed in recent years [@Crowe2022].

### 2.4.1 Length

Total length and blowhole to dorsal fin insertion length were analyzed in terms of age and sex classes.

## Growth curve (length)

Seasonal changes in width

## 2.5 West Coast US dolphins

# RESULTS

## 3.1 Data collection

Data were collected at least once per season in the Patea-Doubtful Sound complex, and were collected in the summer and winter of 2022 and the winter of 2023.

## 3.2 Video analysis

## 3.3 Photogrammetry measurements

## 3.4 Stats

## 3.5 West Coast US dolphins

Not coastal ecotype

Just as big

Similar latitude to Texas and mid-Atlantic

# DISCUSSION

1.  LENGTH
    1.  Similar size as northern counterparts and even other bottlenose along the US west coast.
    2.  Growth in terms of age as well as sex
2.  WIDTH
    1.  Seasonality
    2.  Reproductive
3.  Sexual dimorphism has been suggested to increase with age (Read et al. 1993), but this doesn't seem to hold for at least resident populations of the offshore ecotype (Cheney and this study).
4.  Actual altitude from LiDAR vs barometric altitude when flying

# ACKNOWLEDGEMENTS
