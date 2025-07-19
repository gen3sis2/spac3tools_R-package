\donttest{
library(spac3tools)
# First load some raster spaces.rds input
dir_input <- system.file("extdata/TestSpaces/geodynamic_spaces/h3", package = "gen3sis2")
points_dir <- file.path(tempdir(),"points")
dir.create(points_dir)

file.copy(file.path(dir_input,"spaces.rds"), points_dir)

spaces_h3 <- readRDS(file.path(points_dir,"spaces.rds"))

# Input is of type "h3"
print(spaces_h3$meta$type)

# Conversion to points happens locally
space_h3_to_points(dir_input = points_dir)

spaces_points <- readRDS(file.path(points_dir,"spaces.rds"))

# Now input is of type "points"
print(spaces_points$meta$type)
}

