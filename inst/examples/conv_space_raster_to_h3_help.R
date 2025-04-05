
dir_input="C:/temp/decompressed_spaces/world60by10at4d"
space_raster_to_h3(dir_input,
                   dir_output=file.path(dir_input, "h3_v2"),
                   res=0,
                   verbose=0,
                   compute_distances=TRUE)
