\dontrun{
  ## This is a temperature (Chelsa v2 bio1) raster for South America
  ## Each layer must be one timestep, and they must be ordered
  temperature_raster <- terra::rast(system.file(
    "extdata/rasters/temperature_rasters_1dg.tif", package = "spac3tools")
  )
  # Names are important, as they denote which variable the raster represents
  raster_list <- list(
    temperature = temperature_raster
  )

  # Get a vector with H3 cells of interest
  south_america_cells <- c(
    "8166fffffffffff", "81663ffffffffff", "81677ffffffffff", "815f7ffffffffff",
    "81667ffffffffff", "818f7ffffffffff", "818a7ffffffffff", "818afffffffffff",
    "818b7ffffffffff", "818a3ffffffffff", "815f3ffffffffff", "818abffffffffff",
    "81807ffffffffff", "81803ffffffffff", "8181bffffffffff", "81813ffffffffff",
    "81a8bffffffffff", "81817ffffffffff", "818bbffffffffff", "81a8fffffffffff",
    "818b3ffffffffff", "81a87ffffffffff", "81a83ffffffffff", "81a93ffffffffff",
    "81c2fffffffffff", "81c23ffffffffff", "81c37ffffffffff", "81b2fffffffffff",
    "81cebffffffffff", "81b27ffffffffff", "81a97ffffffffff", "81b37ffffffffff",
    "81b23ffffffffff", "81b33ffffffffff", "818e7ffffffffff", "81c33ffffffffff",
    "81cfbffffffffff", "81df7ffffffffff", "818f3ffffffffff", "8166bffffffffff")

  # Construct and organize a dataset
  h3_list <- data_raster_to_h3(h3_address = south_america_cells, raster_list = raster_list)

  # Create the H3 spaces
  create_spaces_h3(
    h3_list = h3_list,
    cost_function = function(source, dest){
      return(1/1000)
    },
    output_directory = file.path(tempdir(),"h3"),
    full_dists = TRUE,
    overwrite_output = TRUE,
    verbose = TRUE,
    duration=list(from=-4, to=0, by=1, unit="Ma"),
    geodynamic = TRUE
  )
  space <- readRDS(file.path(tempdir(),"h3","spaces.rds"))
  gen3sis2::plot_space_overview(space)
}
