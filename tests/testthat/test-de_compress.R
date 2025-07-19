# Test for compression and decompression functions

test_that("compress_space and decompress_space works", {
  withr::with_tempdir({
    dir_input <- system.file("extdata/TestSpaces/geodynamic_spaces/raster", package = "gen3sis2")
    compressed_output <- file.path(getwd(),"compressed")
    decompressed_output <- file.path(getwd(), "decompressed")

    compress_space(
      dir_input = dir_input, # the dir containing the spaces.rds
      dir_output = compressed_output # the dir to save the compressed spaces.rds
    )

    compressed_files <- list.files(compressed_output)
    expect_equal(compressed_files, "spaces.rds")

    decompress_space(
      dir_input=compressed_output, # just the directory, not the path to compressed spaces.rds
      dir_output=decompressed_output,
      cost_function_index=1, # which cost function to use
      remove_temp_rasters=TRUE # remove the reconstructed raster or not
    )

    expect_equal(list.files(decompressed_output), c("decompressed", "spaces.rds"))
    expect_contains(list.files(file.path(decompressed_output,"decompressed")), c("distances_full","distances_local","spaces.rds"))
  })
})
