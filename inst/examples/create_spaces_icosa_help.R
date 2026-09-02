\dontrun{
  ## This is a temperature (Chelsa v2 bio1) raster for South America
  ## Each layer must be one timestep, and they must be ordered
  ## # COpied over from h3-help
  temperature_raster <- terra::rast(system.file("extdata/rasters/temperature.tiff",package = "spac3tools"))

  # Names are important, as they denote which variable the raster represents
  raster_list <- list(
    temperature = temperature_raster
  )

  # attach
  library(icosa)
  # construct grid
  hex <- hexagrid(spacing=4)
  # structure input
  icosa_list <- data_raster_to_icosa(icosa =hex, raster_list = raster_list)
  # needs to include tessellation vector so the grid can be reconstructed
  # for distance calculations

}
