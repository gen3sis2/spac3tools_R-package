# Test for conversion functions

# landscape_to_spaces ----
test_that("landscape_to_spaces works",{
  withr::with_tempdir({
    landscapes_path <- system.file("extdata/SouthAmerica/landscape/", package = "gen3sis2")
    old_landcape <- readRDS(file.path(landscapes_path, "landscapes.rds"))
    spaces_path <- file.path(getwd(),"spaces")
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

    new_spaces <- readRDS(file.path(spaces_path, "spaces.rds"))

    expect_contains(names(new_spaces),c("env","meta"))
    expect_equal(names(old_landcape), names(new_spaces$env))
  })
})

# space_raster_to_h3 ----
test_that("space_raster_to_h3 and space_h3_to_points works",{
  withr::with_tempdir({
    dir_input <- system.file("extdata/TestSpaces/geodynamic_spaces/raster", package = "gen3sis2")
    dir_output <- file.path(getwd(), "h3")

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

    expect_equal(h3_spaces$meta$type, "h3")
    expect_equal(h3_spaces$meta$type_spec$res, 2)

    space_h3_to_points(dir_input = dir_output)

    points_spaces <- readRDS(file.path(dir_output,"spaces.rds"))
    expect_equal(points_spaces$meta$type, "points")
  })
})
