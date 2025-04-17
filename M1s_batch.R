########################
### General settings ###
########################

random_seed = params$seed
start_time = NA
end_time = NA
max_number_of_species = 1e5
max_number_of_coexisting_species = 1e5
initial_abundance =  30

# a list of traits to include with each species
trait_names = c("dispersal", "temp_opt", "precip_opt")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]

environmental_ranges = list("temp" = NA, "phydiv"=NA, "precip"=NA, "hydro"=NA, "area"=NA)

#########################
### Observer Function ###
#########################

end_of_timestep_observer = function(data, vars, config){
  #save stuff
  #browser()
  #plot_richness(data$all_species, data$landscape)
  save_species()
  save_richness()
  
  # make p/a matrices
  if(!file.exists(file.path(config$directories$output, "pa_matrices"))){dir.create(file.path(config$directories$output, "pa_matrices"))}
  
  # cell names
  all_cells <- rownames(data$landscape$coordinates)
  
  # get 0 for absence and 1 for presence in each grid cell
  all_species_presence <- do.call( cbind, lapply(data$all_species, FUN = function(x) {ifelse(all_cells %in% names(x$abundance), 1, 0)}))
  
  # colnames are species names
  colnames(all_species_presence ) <- unlist(lapply(data$all_species, function(x){x$id}))
  
  # column bind with x/y coordinates
  presence_absence_matrix <- cbind(data$landscape$coordinates, all_species_presence)
  
  saveRDS(presence_absence_matrix, file=file.path(config$directories$output,"pa_matrices",  paste0("pa_t_",vars$ti, ".rds")))
  
}

######################
### Initialization ###
######################

create_ancestor_species <- function(landscape, config) {
  range <- landscape$extent #extent of Africa an other way to put it : range <- landscape$extent
  #selec all coord inside that range (x and y into this range)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
    co[, "x"] <= range[2] &
    co[, "y"] >= range[3] &
    co[, "y"] <= range[4]
  
  new_species <- list()
  for(i in 1:30){
    initial_cells <- rownames(co)[selection]
    initial_cells <- sample(initial_cells, 1)
    new_species[[i]] <- create_species(initial_cells, config)
    #set local adaptation to max optimal temp equals local temp
    new_species[[i]]$traits[ , "temp_opt"] <- landscape$environment[initial_cells,"temp"]
    new_species[[i]]$traits[ , "precip_opt"] <- landscape$environment[initial_cells,"precip"]
    new_species[[i]]$traits[ , "dispersal"] <- 1 
    #plot_species_presence(landscape, species=new_species[[i]])
  }
  
  return(new_species)
}


#################
### Dispersal ###
#################

dispersal_scale <-  params$dispersal_scale
dispersal_shape <-  params$dispersal_shape

get_dispersal_values <- function(n, species, landscape, config) {
  values <- rweibull(n, shape = dispersal_shape, scale = dispersal_scale)
  return(values)
}

##################
### Speciation ###
##################

# threshold for genetic distance after which a speciation event takes place
divergence_threshold =  params$divergence_threshold

# factor by which the divergence is increased between geographically isolated population
# can also be a matrix between the different population clusters
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  divergences <- vector("numeric", length(unique(cluster_indices)))
  
  divergence_matrix <- matrix(0, ncol=length(unique(cluster_indices)),  nrow=length(unique(cluster_indices)))
  cluster_indices_order <- unique(cluster_indices)[order(unique(cluster_indices))]
  # order them, find pairwise phydiv dist
  for(i_x in 1:length(cluster_indices_order)){
    mean_phydiv_i <- mean(landscape$environment[names(species$abundance[which(cluster_indices %in% cluster_indices_order[i_x])]),"phydiv"], na.rm=T)
    for(j_x in 1:length(cluster_indices_order)){
      mean_phydiv_j <- mean(landscape$environment[names(species$abundance[which(cluster_indices %in% cluster_indices_order[j_x])]),"phydiv"], na.rm=T)
      divergence_matrix[i_x, j_x] <- (sum(mean_phydiv_i, mean_phydiv_j)/(2))
    }
  }
  rownames(divergence_matrix) <- cluster_indices_order
  colnames(divergence_matrix) <- cluster_indices_order
  
  if(length(divergence_matrix) == 1){ divergence_matrix <- as.numeric(divergence_matrix)}else{
    
    diag(divergence_matrix) <- 0
  }
  
  return(divergence_matrix)
}

#######################
### Trait Evolution ###
#######################

apply_evolution <- function(species, cluster_indices, landscape, config) {
  #browser()
  # cell names
  traits <- species[["traits"]]
  cells <- rownames( traits )
  
  #browser()
  
  # evolve traits for each cluster
  for(cluster_index in unique(cluster_indices)){
    
    cells_cluster <- cells[which(cluster_indices == cluster_index)]
    
    t_theta_cluster <- median(landscape$environment[cells_cluster,"temp"], na.rm=T) #clusters local temp values: Ts
    
    # evolve temperature and precipitation niche optima
    delta_temp_cluster <- abs(rnorm(1, mean = 0, sd = 0.005))
    
    # get the trait values
    t_opt_cluster <- traits[cells_cluster, "temp_opt"][1]   #clusters of temperature niche center: Ti
    
    # find the difference between environmenta and trait value
    # if env is less than niche it will be negative
    # if env is more than niche it will be positive
    diff_t <- t_theta_cluster - t_opt_cluster    # Ts - Ti 
    
    # new value will be in the direction of theta
    new_value_t <- t_opt_cluster + delta_temp_cluster * sign(diff_t)
    
    # however if the difference is small relative to the rate (sigma_e), the evolutionary change can overshoot the optima
    # if the new niche is more poorly adapted than the original then we need to change it (under the assumption that the niche is under positive selection)
    new_diff_t <-  t_theta_cluster - new_value_t 
    
    # so if the new difference is > than the original difference, draw the value back to match the original difference (this time on the other side of the optima)
    if(abs(new_diff_t) > abs(diff_t)){new_value_t = t_theta_cluster - (sign(new_diff_t) * abs(diff_t))}
    
    # add new values back to traits matrix
    traits[cells_cluster, "temp_opt"]  <- new_value_t
    
  }
  
  # set bounds between 0 and 1 so the species can;t evolve a niche beyond that present in the data (all temp is scaled between 0 and 1)
  if(any(traits[, "temp_opt"] > 1)){traits[which(traits[,"temp_opt"]>1), "temp_opt"] <- 1}
  if(any(traits[, "temp_opt"] < 0)){traits[which(traits[,"temp_opt"]<0), "temp_opt"] <- 0}
  
  return(traits)
  
}

#################################################
### Environmental and Ecological Interactions ###
#################################################
library(deSolve)

myfun <- function(t, N, parameters) {
  n <- length(N)
  dN <- rep(0, n)
  N <- sapply(N, FUN=function(x) max(0, x))
  with(as.list(c(N, parameters)),{
    for (i in 1:n){
      Ki <- exp(-(temp_opt[i]-temp_site)^2/(temp_omega^2)) * exp(-(precip_opt[i]-precip_site)^2/(precip_omega^2)) *(K - sum(N[-i]))
      dN[i] <- N[i] * (Ki - N[i])
    }
    return(list(dN))
  })
}

temp_omega <- params$ecology_omega # environmental filtering parameter [0, n]
precip_omega <- params$ecology_omega # environmental filtering parameter [0, n]

inflection <- 0.2 # extinction parameter: inflexion point of extinction probability curve
decay <- 100  # extinction parameter: decay rate of extinction probability curve

K_max <- 1
alpha =0.5
beta = 0.5

apply_ecology <- function(abundance, traits, landscape, config) {
  #browser()
  n <- length(abundance)
  # define the traits
  temp_opt <- traits[, "temp_opt"]
  precip_opt <- traits[, "precip_opt"]
  temp_site <- landscape[, "temp"]
  precip_site <- landscape[, "precip"]
  hydro_site <- landscape[, "hydro"]
  area_site <- landscape[, "area"]
  
  P = precip_site
  H = hydro_site
  area = area_site
  
  K <- area * pmin(K_max, (alpha * P/0.31 + beta * H/1.31) * K_max)
  
  # fill in parameters
  parameters <- list("temp_opt" = temp_opt,
                     "temp_site" = temp_site,
                     "precip_opt" = precip_opt,
                     "precip_site" = precip_site,
                     "hydro_site" = hydro_site,
                     "K" = K, 
                     "temp_omega" = temp_omega,
                     "precip_omega" = precip_omega,
                     "area_site"= area_site) 
  N0 <- rep(1, n)
  times <- c(0,100)
  out <- ode(y = N0, times = times, func = myfun, parms = parameters)
  
  abundance <- out[2,]
  
  # low abundance (<0.01) goes extinct (could also do the probabilistic extinction function here):
  # e.g., 
  #prob_extinction <- (1/(1+exp(-decay*(inflection - abundance))))
  #abundance[which(sapply(prob_extinction, FUN=function(x){sample(c(0,1),1, prob= c(x, 1-(x)))})==0)] <- 0
  
  abundance[which(abundance < 0.01)] <- 0
  
  return(abundance)
  
}
