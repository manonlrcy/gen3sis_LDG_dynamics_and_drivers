# Load libraries
library(ncdf4)
library(gen3sis)
library(raster)

# Set directories
output_directory <- "set_your_land_var_repository"
data_directories <- c(
  PDI = "your_variable.nc",
  TEMP = "your_variable.nc",
  RAIN = "your_variable.nc",
  hydro = "your_variable.nc"
)

# Set variables
time_interval <- 1
geo_period <- 150
n_timesteps <- geo_period + time_interval

# Stack environmental variables
stack_variables <- function(directory, varname) {
  files <- list.files(directory, pattern = '*.nc', full.names = TRUE)
  sorted_files <- mixedsort(files)
  stack(sorted_files, varname = varname)
}

Precip <- stack_variables(data_directories$RAIN, "precip")
Phydiv <- stack_variables(data_directories$PDI, "phydiv")
Temp <- stack_variables(data_directories$TEMP, "temp")
Hydro <- stack_variables(data_directories$hydro, "hydrocat")

# Set names for time steps
time_steps <- seq(0, geo_period, by = time_interval)
names(Precip) <- names(Temp) <- names(Phydiv) <- names(Hydro) <- time_steps

# Set CRS
crs_value <- "+proj=longlat +datum=WGS84"
crs(Precip) <- crs(Temp) <- crs(Phydiv) <- crs(Hydro) <- crs_value

# Create area raster stack
area_ras_4D <- area(Hydro[[1]])
crs(area_ras_4D) <- crs_value
ts_area_stack_4D <- stack(area_ras_4D)

for (i in 1:n_timesteps) {
  ts_area_stack_4D[[i]] <- mask(area_ras_4D, Hydro[[i]])
}
names(ts_area_stack_4D) <- time_steps

# Make a list of environmental variables
landscape_list <- list(
  precip = unstack(Precip),
  temp = unstack(Temp),
  phydiv = unstack(Phydiv),
  hydro = unstack(Hydro),
  area = unstack(ts_area_stack_4D)
)

# Gen3sis Input
cost_function_water <- function(source, habitable_src, dest, habitable_dest) {
  if (!all(habitable_src, habitable_dest)) {
    return(2/1000)
  } else {
    return(1/1000)
  }
}

create_input_landscape(
  landscapes = landscape_list,
  timesteps = time_steps,
  cost_function = cost_function_water,
  directions = 8,
  output_directory = file.path(output_directory, "landscape1"),
  overwrite_output = TRUE,
  crs = crs_value,
  calculate_full_distance_matrices = FALSE,
  verbose = 1
)

# Dispersal including Phydiv - Model 1d
cost_function_M2 <- function(source, habitable_src, dest, habitable_dest) {
  if(!all(habitable_src, habitable_dest)) {
    return(2/1000)
  } else {
    return((max(source["phydiv"], dest["phydiv"])+1)/1000)
  }
}

create_input_landscape(
  landscapes = landscape_list,
  cost_function = cost_function_M2,
  directions = 8,
  output_directory = file.path(output_directory, "landscape2"),
  timesteps = time_steps,
  overwrite_output = TRUE,
  crs = crs_value,
  calculate_full_distance_matrices = FALSE,
  verbose = 1
)

# References
citation("gen3sis")

# Print session info
devtools::session_info()
