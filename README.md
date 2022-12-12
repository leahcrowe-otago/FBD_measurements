# FBD_measurements

# Workflow:

1. Identify when in video dolphin surfaces, ID that dolphin
2. Screen grab the best image of the dolphin above that shows rostrum to tail at/near the surface
  - used VLC player to go frame by frame
3. Add screen grab filename to the time and ID from #1
4. Run "Altitudes for snapshots" for each day of droning
  - this adds a check to make sure GPS of the LiDAR and the metadata from the drone match up
  - finds altitude for each snapshot
  - matches up to IDs and outputs file to enter measurements
  - identifies if altimeter reading meets criteria per Dickson et al. 2020 ("issues" column)
5. Populate the "altperimage_*" file with measurements (also consider instead merging output csvs from Whalelength)
  - filter by issue == "N"
  - hide columns C, G, H, and J to AC
  - Run each resulting image through Whalength
    - I1P Video stills, -1.5 cm offset
6. Run "measurements_demo.Rmd" to group by age and sex classes
