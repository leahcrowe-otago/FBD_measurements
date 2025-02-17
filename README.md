# FBD_measurements
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14563512.svg)](https://doi.org/10.5281/zenodo.14563512)

Poster presented at the 2024 Society for Marine Mammalogy conference:
[Crowe LM, Schofield MR, Dawson SM, Rayment WJ (Nov 2024). Growth estimates of bottlenose dolphins (_Tursiops truncatus_) at the southern extreme of their range. Poster presentation at the 25th Biennial Conference on the Biology of Marine Mammals, Perth, Western Australia, Australia.](https://github.com/leahcrowe-otago/FBD_measurements/blob/main/outputs/SMM%20abstract%202024/E%26E%2061_CROWE_LEAH_SMM2024.pdf)

# Workflow:

1. Identify when in video dolphin surfaces, ID that dolphin
2. Screen grab the best image of the dolphin above that shows rostrum to tail at/near the surface
  - used VLC player to go frame by frame
3. Add screen grab filename to the time and ID from #1
4. Run "./scripts/Altitudes for snapshots" for each day of droning
  - this adds a check to make sure GPS of the LiDAR and the metadata from the drone match up
  - finds altitude for each snapshot
  - matches up to IDs and outputs file to enter measurements
  - identifies if altimeter reading meets criteria per Dickson et al. 2020 ("issues" column)
5. Populate the "altperimage_*" file with measurements (.csvs from *Whalelength* are later merged in "./scripts/whalength_file_merge.R")
  - filter by issue == "N"
  - hide columns C, G, H, and J to AC
  - Run each resulting image through Whalength
    - I1P Video stills, -1.5 cm offset
6. Run "measurements_demo.Rmd" to wrangle data, incorporate life history info ("demo"graphy data), and create some supplementary figures
  - calls "whalength_file_merge.R"
    - merges all *Whalelength* outputs
    - creates Fig. S4, calibration measurements 
  - Supplementary Figures S5
  - data output: see data repo for metadata details
    -  individual *i* measurement data at each capture occasion *j* for model
    -  individual life history data
7. Model run in the "./scripts/stan/runstan_allo_mv_t0.R" file: 'run' the 'stan' code of 'allo'metric measurement data in a 'm'ulti'v'ariate model with fixed age at length zero (t0)
  - data formatting for model
  - calls the Stan model
    - main model: "./scripts/stan/vb_mod_all0_t0.stan"
    - for the supplementary sex/pod effect model: "./scripts/stan/sex_pod_effects.R" runs the model specified in "./scripts/stan/vb_mod_all0_t0_sexpod.stan"
  - initial values for von Bertalanffy model (init_vb)
  - fit the model (fit_vb)
  - save results
8. "./scripts/Results.R"
  - read data
  - summarise results
  - create table S2
  - create Figs. 2–4, S6, S7
9. "Supplementary.qmd" creates the suppmat file
  - creates supplementary tables
  - Section 1.3 created from './scripts/Results_age_adjusted.R'
  - Fig. S1 map created by calling "./scripts/ms_map.R"
10. Data: "./data/Measurements/Data for review"
  - Metadata-description.pdf: describes the data in the following files:
    - ij_1.csv
    - ij_2.csv
    - ij_3.csv
    - ij_ID.csv
