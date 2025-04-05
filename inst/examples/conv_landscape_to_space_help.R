dir_l <- "C:/temp/decompressed_spaces/globe_0-20Ma_1Myr/"

landscape_to_space(dir_l,
                   duration=list(from=-20, to=0,  by=1, unit="Ma"),
                   crs="+proj=longlat +datum=WGS84 +no_defs",
                   cost_function=list(xx=gcf$XXHarderNA_dist_Km)
                   )



