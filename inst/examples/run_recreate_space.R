## METADATA ===============================================================
## Description: recreate gen3sis input space from previouly avaiable landscape
## landscape.
##
## R version: 4.2.2 for Windows
## Date: 2023-01-3 11:07:00
## License: GPL3
## Author: Oskar Hagen (oskar@hagen.bio)
##=======================================================================##
dir_input <- "C:/temp/4ds/world60by10at4d"

############# RUN FOR static landscape (geodynamic=TRUE) -----------
dir_gen2 <- file.path(dir_input, "gen3sis2_test2")

library(spac3tools)

landscape_to_space(dir_input,
                   dir_output=dir_gen2,
                   duration=list(from=-60, to=0,  by=10, unit="Ma"),
                   crs="+proj=longlat +datum=WGS84 +no_defs",
                   cost_function=list(xx=gcf$XXHarderNA_dist_Km))

#unlink(file.path(dir_gen2, "decompressed"))

decompress_space(dir_input=dir_gen2,
                 cost_function_index=1,
                 dir_output=dir_gen2,
                 remove_temp_rasters=TRUE)

# run a simulation as raster
set.seed(102)
library(gen3sis2)
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape=file.path(dir_gen2, "decompressed"),
               output_directory=file.path(dir_input, "outputs"))

# convert raster to h3
space_raster_to_h3(dir_input=file.path(dir_gen2, "decompressed"),
                   dir_output=file.path(dir_input, "h3"),
                   res=2,
                   verbose=1,
                   compute_distances=TRUE)

set.seed(107)
# run a simulation as h3
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape=file.path(dir_input, "h3_v2"),
               output_directory=file.path(dir_input, "outputs"), verbose = 3)




############# RUN FOR static landscape (geodynamic=FALSE) -----------
dir_input_static <- "c:/temp/1d_rep10/"

dir_gen2 <- file.path(dir_input_static, "gen3sis2_test2")

library(spac3tools)


landscape_to_space(dir_input_static,
                   dir_output=dir_gen2,
                   duration=list(from=-9, to=0,  by=1, unit="Xa"),
                   crs="+proj=longlat +datum=WGS84 +no_defs",
                   cost_function=list(xx=gcf$XXHarderNA_dist_Km))

#unlink(file.path(dir_gen2, "decompressed"))

decompress_space(dir_input=dir_gen2,
                 cost_function_index=1,
                 dir_output=dir_gen2,
                 remove_temp_rasters=TRUE)


# run a simulation as raster
set.seed(102)
library(gen3sis2)
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape=file.path(dir_gen2, "decompressed"),
               output_directory=file.path(dir_input, "outputs"))

# convert raster to h3
space_raster_to_h3(dir_input=file.path(dir_gen2, "decompressed"),
                   dir_output=file.path(dir_gen2, "decompressed", "h3_3"),
                   res=3,
                   verbose=3,
                   compute_distances=TRUE)

set.seed(107)
# run a simulation as h3
run_simulation(config=file.path(dir_input, "configs", "config_M0HR.R"),
               landscape="C:/Users/am92guke/Documents/iDiv/code/eco2evo/evo2eco_JOBZ/data/temp/h3_2", #file.path(dir_gen2, "decompressed", "h3_3")
               output_directory=file.path(dir_gen2, "outputs"), verbose = 3)


## REMOVE AFTER NOT NEEDED
#
# ### [] create variables for directory relative location -------
# # dd= data directory
# # dd <- "C:/temp"
# dir.create(od, showWarnings = FALSE)
#
# # input_landscape <- file.path("C:/Dropbox/Angiosperms_gen3sis/input/", "w2004d/landscapes.rds")
# recreated_space_dir <- file.path(od,"recreated_ideal")
# h3_space_dir <- file.path(od,"h3_ideal_test")
#
# # read
# gls=readRDS(input_landscape)
#
# # lapply(gls, function(x){x[,1:4]})
# #### TODO FIX THIS HERE, ON HOW TO CALL PROPERLY FOR TEMPORARY NAMES...
#
# # n_lc <- lapply(gls, function(x){
# #   x[,c("x", "y", as.character(seq(0,390, 6)))]
# #   })
# # t_steps <- colnames(x[,-c(1,2)]))
#
#
# decompress_space(dir_input="C:/temp/rasterfree_input/new_test",
#                  cost_function=gcf$XXHarderNA_dist_Km,
#                  dir_output="C:/temp/decompressed_spaces",
#                  remove_temp_rasters=TRUE)
#
#
#
# ######### gen3sis reference landscape  ------------
# # rf_od <- file.path(od, "input_rasterfree")
# # lsname
# ## _load ----------
# # define rezz for hexagons..
# rezz <- 1
# compute_distances <- TRUE
# plot_stuff <- TRUE
# # orig_gen3sis_space
# # dir origin
# dir_o_gls <- recreated_space_dir
# # dir desteny
# dir_d_gls <- h3_space_dir
# o_gls_all <- readRDS(file.path(dir_o_gls, "landscapes.rds"))
# # origin full distances
# o_cd_fl  <- list.files(file.path(dir_o_gls, "distances_full"))
# # order o_cd_fl
# numb_ts <- as.numeric(gsub("[^0-9.-]+", "", o_cd_fl))
# o_cd_fl <- o_cd_fl[order(numb_ts)]
#
# # landscape output dir
#
# # create new input dir
# if (!dir.exists(dir_d_gls)){
#   dir.create(dir_d_gls, showWarnings = FALSE)
#   dir.create(file.path(dir_d_gls, "distances_full"), showWarnings = FALSE)
# }
# # dir.create("input_h3") # this is the new input
# # subset if necessary
# # ls_sub_time <- lapply(ls, function(x){
# #   nx <- x[,1:12]
# # })
# # ls <- ls_sub_time
#
# #### LOOP OVER TIME -------------------
# time_steps <- colnames(o_gls_all[[1]][-c(1,2)])
# n_ts <- length(time_steps)
#
# if (n_ts!=length(o_cd_fl)){
#   error("Mismatch on number of time steps")
# }
# # Get parent cell indexes
# h3_0 <- h3jsr::get_res0()
# # Get points of h3_0 at desired resolution
# c1 <- h3jsr::get_children(h3_address = h3_0, res = rezz)
# # Get hex axes
# h3c1 <- unlist(c1)
# p1 <- h3jsr::cell_to_point(h3_address = h3c1, simple = FALSE)
# poly1 <- cell_to_polygon(input= h3c1, simple=FALSE)
# # prepare variables
# all_pts <- o_gls_all[[1]][,c("x", "y")]
# # set as coordinates
# all_pts_sf <- st_as_sf(all_pts, coords = c("x","y"))
# sf::st_crs(all_pts_sf) <- crs_wgs84 # set coordinates
# # make H3 simple
# h3rezz_pts <- as.data.frame(st::st_coordinates(p1$geometry))
# cat("# cells is: ")
# (unlist(lapply(list("all_pts"=all_pts, "h3rezz_pts"=h3rezz_pts), nrow)))
# #####  run loop over each h3rezz point
# np <- length(p1$h3_address) # or nrow(h3rezz_pts) as above
# h3index <- p1$h3_address
# all_r_index <- rownames(gls[[1]])
# # closest raster index
# c_r_index <- rep(NA, np)
# # store error distance in meters
# #create numerical empty vector
# error_m <- c_r_index
#
#
# ## short initial loop to finds closest points... pnly do once
# for (h3pi in 1:np){ # loop for each point
#   # h3pi <- 12
#   print(paste0("h3pi or h3 index of points is: ", h3pi, "/", np))
#   #get closest point
#   distances <- distGeo(h3rezz_pts[h3pi,], all_pts)
#   m_dist <- min(distances)
#   closest_p_index <- all_r_index[distances==m_dist][1] #SHOULD ALLWAYS BE ONE TODO: implemente a check, TELL this info and sample 1
#   c_r_index[h3pi] <- closest_p_index
#   error_m[h3pi] <- m_dist
# }
# # prepare input variables...
# dummyM <- matrix(rep(NA, np*n_ts), ncol=n_ts)
# rownames(dummyM) <- c_r_index # TODO USE h3c1 index later here#
# colnames(dummyM) <- time_steps
# envs <- list("gtemp"=dummyM, "ghumd"=dummyM, "area"=dummyM)
# # gen_lc <- envs
#
# for (ti in 1:n_ts){ # LOOP OVER TIMESTEPS
#   # ti <- 1
#   print(paste("ti=", ti, ",   at time_step", time_steps[ti]))
#
#   # store gen3sis value
#   # get envs of closest index
#   envs$gtemp[,ti] <- o_gls_all$temp[c_r_index,time_steps[ti]]
#   envs$ghumd[,ti] <- o_gls_all$prec[c_r_index,time_steps[ti]]
#   # set the presences
#   lp1_m <- as.logical(as.numeric(envs[[1]][,time_steps[ti]])) # per time step
#   # set area of points falling into water to NA
#   envs$area[,time_steps[ti]] <- cell_area(h3c1, 'km2', simple=TRUE)[lp1_m]
#
#   for (h3pi in 1:np){ # loop for each point
#     # h3pi <- 1
#     print(paste0("h3pi or h3 index of points is: ", h3pi, "/", np))
#     # TODO make summary for hexagons...
#     # poly1_i <- poly1$geometry[h3pi]
#     # plot(poly1_i)
#
#     # WIP HERE.....
#     # sp::over(all_pts_sf, poly1_i, fn = NULL)
#   }
#   # get new data-----
#
#
#
#   # load distances and manipulate it...
#   if (compute_distances){
#     # lfd_ti = landscapes full distances at ti
#     lfd_ti <- readRDS(file.path(dir_o_gls, "distances_full", o_cd_fl[ti]))
#     print(paste("Loaded:", o_cd_fl[ti]))
#     # sub select cost distances
#     ib <- colnames(lfd_ti)%in%c_r_index[lp1_m] # TODO USE h3c1 index later here too, do match#
#     # ix <- which(c_r_index%in%colnames(cd))
#     # length(ix)
#     new_cost_dist_full <- lfd_ti[ib,ib]
#     # save cost function ti
#     saveRDS(new_cost_dist_full, file.path(dir_d_gls, "distances_full", o_cd_fl[ti]))
#     print(paste("Saved:", file.path(dir_d_gls, "distances_full", o_cd_fl[ti])))
#     # plot entire data with custom function
#   }
#
#   if (plot_stuff){
#     #possible points
#     # points(all_pts, pch=2, col=rgb(1,1,1,1,1))
#     # plot points in original raster
#
#     # plot all possible points
#     plot_points(xy=h3rezz_pts,
#                 vcol=1,
#                 cols=terrain.colors(2),
#                 vcex=1)
#
#     # plot only points with temperature on h3
#     plot_points(xy=h3rezz_pts[lp1_m, ],
#                 vcol=envs$gtemp[lp1_m],
#                 cols=gen3sis::color_richness(5),
#                 vcex=1.5)
#     # plot points with temperature on original raster input
#     plot_points(xy=all_pts_sf$geometry[!is.na(gls$temp[,3]) ],
#                 vcol=na.omit(gls$temp[,3]),
#                 cols=gen3sis::color_richness(30),
#                 vcex=1)
#
#     # plot cost distance values with the original location...
#     if (compute_distances) {
#       # points(all_pts[colnames(new_cost_dist_full),], pch=3)
#       # # plot cost distance location....
#       # points(all_pts[colnames(lfd_ti),], pch=3)
#       #
#       # # plot random distances
#       #
#       # # plot only points with temperature and plot all points with distances on top...
#       # df <- all_pts[colnames(new_cost_dist_full),]
#       # plot_points(xy=h3rezz_pts,
#       #             vcol=envs$gtemp, # use 1 for all points
#       #             cols=gen3sis::color_richness(5),
#       #             vcex=1.5)
#       #
#       # #  plot(df, col=rainbow(800)) # :) have a nice day
#       #
#       # # select random points...
#       # points <- sample(colnames(new_cost_dist_full), 2)
#       # points(df, pch=3)
#       # segments(x0=df[points[1], "x"],y0=df[points[1], "y"],x1=df[points[2], "x"],y1=df[points[2], "y"])
#       # round(new_cost_dist_full[points[1],points[2]], 2)
#       # text(50, 50, )
#       #
#       #
#       # df <- all_pts[colnames(new_cost_dist_full),]
#       # plot_points(xy=df,
#       #             vcol=1, # use 1 for all points
#       #             cols=gen3sis::color_richness(5),
#       #             vcex=1.5)
#     }
#   }
# } # END LOOP OVER TIME STEPS
# # save final landscape
# final_landscape <- lapply(envs, function(x){
#   h3pts <- h3rezz_pts
#   rownames(h3pts) <- c_r_index
#   colnames(h3pts) <- c("x", "y")
#   return(cbind(h3pts, x))
# })
# saveRDS(final_landscape, file.path(dir_d_gls,"landscapes.rds"))
# # TODO READ AND UPDATE METADATA.TXT
# print("Please update your meta data and store with with the landscapes.rds file")
