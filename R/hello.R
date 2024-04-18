# # Hello, world!
# #
# # This is an example function named 'hello'
# # which prints 'Hello, world!'.
# #
# # You can learn more about package authoring with RStudio at:
# #
# #   http://r-pkgs.had.co.nz/
# #
# # Some useful keyboard shortcuts for package authoring:
# #
# #   Install Package:           'Ctrl + Shift + B'
# #   Check Package:             'Ctrl + Shift + E'
# #   Test Package:              'Ctrl + Shift + T'
#
#
# dir_landis <- "C:/temp/4ds/world60by10at4d"
# dir_decompressed <- "C:/temp/decompressed_spaces/world60by10at4d"
#
#
# my_description <- "Temperatures on current koppen bands were extracted for each Koppen band (5Ma resolution)(Hagen et al. 2019)
# and had a focal mean applied at a strength that mimics empirical present temperature(WorldClim2 2018) spread for the koppen zones.
# These were first interpolated at a resolution of 1Ma. Lapse rates for each zones (Hagen et al. 2019) were applied to the
# respective elevation maps (Straume et al. 2020) also available at a resolution of 1Ma. Air surface temperatures were applied.
# LGM and LTG had strength corrected to match reference global air surface temperature maps (Westerhold et al. 2020).
# Global temperature differences were calculated using entire koppen band to account for sea surface temperature
# (Westerhold et al. 2020). For more details see (Hagen et al 2020/2021)."
#
# landscape_to_space(dir_input = dir_landis,
#                    duration=list(from=-60, to=0,  by=10, unit="Ma"),
#                    area.unit="km2",
#                    crs="+proj=longlat +datum=WGS84 +no_defs",
#                    cost_function=list(xx=gcf$XXHarderNA_dist_Km),
#                    author="Oskar Hagen",
#                    source="10.1371/journal.pbio.3001340",
#                    description=list(env="temperature in degree celcius;
#                                     humidity proxy from zero to one (respectively dry and arid)",
#                                     methods=my_description)
# )
#
# # decompress space
# decompress_space(dir_input=dir_landis,
#                  dir_output=dir_decompressed)
#
# # convert a space.rds type raster to h3
# space_raster_to_h3(dir_input=dir_decompressed, res=1)
