######################################
###            Daisyworld          ###
######################################
# Daisyworld minimal gen3sis2 configuration
#
# This configuration implements a minimal, documented version of the
# Daisyworld thought experiment (Watson & Lovelock, 1983). In Daisyworld,
# two daisy types (typically "black" and "white") differ in albedo. The
# area-weighted albedo of the planetary surface alters the planetary
# temperature, and local temperature in turn affects daisy growth. This
# feedback can lead to self-regulation of planetary temperature.
#
# This file is intentionally simple and designed to be read by
# gen3sis2::create_input_config(config_file = "path/to/this/file.R") and then
# used with gen3sis2::run_simulation().
#
# Key concepts encoded here (narrative + parameters):
# - species trait `albedo`: per-species reflectance value (0..1). Black daisies
#   have low albedo (absorb heat), white daisies have high albedo (reflect heat).
# - `bare_ground_albedo`: the albedo of habitat not covered by daisies.
# - `solar_luminosity`: a tunable external forcing controlling incoming energy.
# - Local temperature -> growth: species grow best when local temperature
#   matches their `opt_temp` trait; growth determines coverage which feeds back
#   into planetary albedo.

## GENERAL settings -------------------------------------------------------
# random seed for reproducibility
random_seed <- 42
# start and end times (these follow the vignette: -999..0 kyr)
start_time <- 999
end_time <- 0
# maximum numbers (safety limits)
max_number_of_species <- 20000
max_number_of_coexisting_species <- 5000

## Model parameters specific to Daisyworld ---------------------------------
# The following parameters control the simple planetary energy balance
# implemented outside this config (see notes below). You can tune them to
# replicate classic Daisyworld behaviour (e.g. increase solar_luminosity).
solar_luminosity <- 1.0   # relative incoming energy (1.0 is an arbitrary baseline)
bare_ground_albedo <- 0.4 # albedo of empty ground (between typical black/white values)

# trait names used by the simulation
# - `albedo`: per-species reflectance (0..1). Use ~0.2 for black daisies, ~0.75 for white.
# - `opt_temp`: preferred temperature for growth (°C)
# - `dispersal`: dispersal trait used by the kernel
trait_names <- c("dispersal", "opt_temp", "albedo")
# ## environmental ranges for scaling/visualisation (optional)
# environmental_ranges <- list("mean_temp" = c(-50, 50), "radiation" = c(0, 1))

## OBSERVER (optional) ----------------------------------------------------
end_of_timestep_observer <- function(data, vars, config) {
  # Minimal observer: do nothing heavy by default. You can implement a
  # per-timestep planetary albedo -> temperature update here if you want the
  # full Daisyworld feedback inside the simulation loop.
  invisible(NULL)
}

## INITIALIZATION ---------------------------------------------------------
# initial abundance in colonized cells
initial_abundance <- 10

# Daisyworld initialization: create two ancestor species (black & white daisies)
# Each ancestor is given an `albedo` trait (low for black, high for white) and
# an `opt_temp` close to the local mean temperature at the founding cell.
# This mirrors the classic setup where two daisy types compete spatially and
# affect planetary albedo via their coverage.
create_ancestor_species <- function(space, config) {
  new_species <- list()
  for(i in 1:2){
    initial_cells <- sample(rownames(space$coordinates),3)
    new_species[[i]] <- create_species(initial_cells, config)
    if (i == 1){# Black daisy: low albedo (absorbs heat)
      albedo <- 0.20
      t_opt <- 15
    } else{
      albedo <- 0.75
      t_opt <- 35
    }
    #set local adaptation to max optimal temp equals local temp
    new_species[[i]]$traits[ , "dispersal"] <- 1
    new_species[[i]]$traits[ , "opt_temp"] <- t_opt
    new_species[[i]]$traits[, "albedo"] <- albedo
  }

  return(new_species)
}

## DISPERSAL --------------------------------------------------------------
# basic dispersal kernel returning n dispersal distances
get_dispersal_values <- function(n, species, space, config) {
  # small exponential kernel for this toy example
  rexp(n, rate = 1/5)
}

## SPECIATION -------------------------------------------------------------
# divergence threshold (time units in the config; keep low for example)
divergence_threshold <- 10
get_divergence_factor <- function(species, cluster_indices, space, config) {
  # no special factor, return 1
  0
}

## MUTATION / EVOLUTION ---------------------------------------------------
apply_evolution <- function(species, cluster_indices, space, config) {
  traits <- species[["traits"]]
  return(traits)
}

## ECOLOGY ----------------------------------------------------------------
apply_ecology <- function(abundance, traits, space, config) {
  # THIS IS NOT FINISHED
  abundance_scale = 10
  abundance_threshold = 1
  #abundance threshold
  survive <- abundance>=abundance_threshold
  abundance[!survive] <- 0
  abundance <- (( 1-abs( traits[, "opt_temp"] - space[, "temp"]))*abundance_scale)*as.numeric(survive)
  #abundance threshold
  abundance[abundance<abundance_threshold] <- 0
  k <- ((space[,"area"]*(space[,"arid"]+0.1)*(space[,"temp"]+0.1))*abundance_scale^2)
  total_ab <- sum(abundance)
  subtract <- total_ab-k
  if (subtract > 0) {
    # print(paste("should:", k, "is:", total_ab, "DIFF:", round(subtract,0) ))
    while (total_ab>k){
      alive <- abundance>0
      loose <- sample(1:length(abundance[alive]),1)
      abundance[alive][loose] <- abundance[alive][loose]-1
      total_ab <- sum(abundance)
    }
    #set negative abundances to zero
    abundance[!alive] <- 0
  }

  return(abundance)
}

## NOTES -----------------------------------------------------------------
# To use:
# config_object <- gen3sis2::create_input_config(config_file = "path/to/config_daisyworld.R")
# dirs <- gen3sis2::prepare_directories(config_file = "path/to/config_daisyworld.R", input_directory = "path/to/inst/extdata/Daisyworld")
# gen3sis2::run_simulation(config = config_object, landscape = file.path(dirs$input, "spaces.rds"), output_directory = dirs$output)

# End of Daisyworld minimal config
