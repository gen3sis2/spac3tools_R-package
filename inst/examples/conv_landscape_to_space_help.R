# deprecated
# dir_l <- "C:/temp/decompressed_spaces/globe_0-20Ma_1Myr/"
#
# landscape_to_space(dir_l,
#                    duration=list(from=-20, to=0,  by=1, unit="Ma"),
#                    crs="+proj=longlat +datum=WGS84 +no_defs",
#                    cost_function=list(xx=gcf$XXHarderNA_dist_Km)
#                    )
\donttest{
library(spac3tools)

# Load gen3sis landscapes.rds
landscapes_path <- system.file("extdata/SouthAmerica/landscape/", package = "gen3sis2")
old_landcape <- readRDS(file.path(landscapes_path, "landscapes.rds"))

# It contains the old gen3sis format, with only environmental variables
names(old_landcape)

# It have 66 timesteps and x and y coordinates
ncol(old_landcape$temp)

# Convert the old landscapes.rds to the new spaces.rds
# Need to provide a cost function and duration (timesteps)
spaces_path <- file.path(tempdir(),"space")
dir.create(spaces_path)

landscape_to_space(
  dir_input = landscapes_path,
  dir_output = spaces_path,
  duration = list(from = 65, to = 0, by = -1, unit = "Ma"),
  cost_function = function(source, dest) {
    if(!all(source$habitable, dest$habitable)) {
      return(2/1000)
    } else {
      return(1/1000)
    }
  }
)

# Loading the converted spaces.rds
new_spaces <- readRDS(file.path(spaces_path, "spaces.rds"))

# The new spaces.rds from gen3sis2 is very different, containing more information
names(new_spaces)
names(new_spaces$env)
names(new_spaces$meta)
}
