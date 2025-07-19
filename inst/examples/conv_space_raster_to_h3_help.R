library(spac3tools)
\donttest{
# First load some raster spaces.rds input
dir_input <- system.file("extdata/TestSpaces/geodynamic_spaces/raster", package = "gen3sis2")
dir_output <- file.path(tempdir(), "h3")

raster_spaces <- readRDS(file.path(dir_input,"spaces.rds"))

# Input is of type "raster"
print(raster_spaces$meta$type)

# Convert it to h3 choosing appropriate h3 resolution (not the same as raster resolution)
space_raster_to_h3(dir_input = dir_input,
                   dir_output = dir_output,
                   res=2,
                   verbose=0,
                   compute_distances=TRUE)

h3_spaces <- readRDS(file.path(dir_output,"spaces.rds"))

# Now input is of type "h3"
print(h3_spaces$meta$type)
}

