# TODO what to do with this file? Remove? Convert into vignette?

## METADATA ===============================================================
## Description: recreate gen3sis input space from previouly avaiable landscape
## landscape.
##
## R version: 4.2.2 for Windows
## Date: 2023-01-3 11:07:00
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##
dir_input <- "C:/temp/4ds/world60by10at4d"

############# RUN FOR static landscape (geodynamic=TRUE) -----------
dir_gen2 <- file.path(dir_input, "gen3sis2_test2")

library(spac3tools)

landscape_to_space(dir_input,
                   dir_output=dir_gen2,
                   duration=list(from=-60, to=0,  by=10, unit="Ma"),
                   crs="+proj=longlat +datum=WGS84 +no_defs",
                   cost_function=list(xx=gcf$XXHarderNA_dist_Km))

#unlink(file.path(dir_gen2, "decompressed"))

decompress_space(dir_input=dir_gen2,
                 cost_function_index=1,
                 dir_output=dir_gen2,
                 remove_temp_rasters=TRUE)

# run a simulation as raster
set.seed(102)
library(gen3sis2)
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape=file.path(dir_gen2, "decompressed"),
               output_directory=file.path(dir_input, "outputs"))

# convert raster to h3
space_raster_to_h3(dir_input=file.path(dir_gen2, "decompressed"),
                   dir_output=file.path(dir_input, "h3"),
                   res=2,
                   verbose=1,
                   compute_distances=TRUE)

set.seed(107)
# run a simulation as h3
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape=file.path(dir_input, "h3_v2"),
               output_directory=file.path(dir_input, "outputs"), verbose = 3)




############# RUN FOR static landscape (geodynamic=FALSE) -----------
dir_input_static <- "c:/temp/1d_rep10/"

dir_gen2 <- file.path(dir_input_static, "gen3sis2_test2")

library(spac3tools)


landscape_to_space(dir_input_static,
                   dir_output=dir_gen2,
                   duration=list(from=-9, to=0,  by=1, unit="Xa"),
                   crs="+proj=longlat +datum=WGS84 +no_defs",
                   cost_function=list(xx=gcf$XXHarderNA_dist_Km))

#unlink(file.path(dir_gen2, "decompressed"))

decompress_space(dir_input=dir_gen2,
                 cost_function_index=1,
                 dir_output=dir_gen2,
                 remove_temp_rasters=TRUE)


# run a simulation as raster
set.seed(102)
library(gen3sis2)
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape=file.path(dir_gen2, "decompressed"),
               output_directory=file.path(dir_input, "outputs"))

# convert raster to h3
space_raster_to_h3(dir_input=file.path(dir_gen2, "decompressed"),
                   dir_output=file.path(dir_gen2, "decompressed", "h3_3"),
                   res=3,
                   verbose=3,
                   compute_distances=TRUE)

set.seed(107)
# run a simulation as h3
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape="C:/Users/am92guke/Documents/iDiv/code/eco2evo/evo2eco_JOBZ/data/temp/h3_2", #file.path(dir_gen2, "decompressed", "h3_3")
               output_directory=file.path(dir_gen2, "outputs"), verbose = 3)

