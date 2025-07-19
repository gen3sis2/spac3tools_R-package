\dontrun{
library(spac3tools)

# The function reconstruct the distance matrices and rasters from the compressed spaces.rds
decompress_space(
  dir_input="your/compressed_spaces/directory", # just the directory, not the path to compressed spaces.rds
  dir_output="your/decompressed_spaces/directory",
  cost_function_index=1, # which cost function to use
  remove_temp_rasters=FALSE # remove the reconstructed raster or not
)
}
