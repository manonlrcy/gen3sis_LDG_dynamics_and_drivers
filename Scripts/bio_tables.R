library(tidyverse)
library(raster)
library(sp)
library(ape)
library(phytools)
library(dplyr)
library(network)

bio_tables <- function(path, out) {
  # Extract folder name
  fname = 1
  fname_keep <- tail(strsplit(path, "/")[[1]], n=fname+1)
  name_keep <- tail(fname_keep, 1)
  
  # Create output directory
  out_dir <- file.path(out, paste0(name_keep))
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  
  for (i in 150:1){
    # Reads presence-absence matrices for consecutive timesteps
    pa_t1 <- readRDS(file.path(paste0(path,"/pa_matrices","/pa_t_",i,".rds")))
    pa_t2 <- readRDS(file.path(paste0(path,"/pa_matrices","/pa_t_",i-1,".rds")))
    
    #####_________RICHNESS____________________________________________________________________________###
    # Converts matrix to dataframe
    pa_df <- as.data.frame(pa_t2)
    # Calculates species richness by summing rows (each location's total species)
    rich_t2 <- pa_df %>% mutate(richness = rowSums(.[c(3:ncol(pa_t2))])) 
    rich_2 <- rich_t2 %>% relocate(richness, .after = y)
    # Keeps only coordinates and richness columns
    richness <- rich_2[c(0:3)]  
    
    # Save richness data
    write_csv(richness, file=file.path(out_dir, paste0("/richness", i-1, ".csv")))
    
    #___________________________________________________________________________________________________#
    
    # Gets species names from both timesteps
    species_t1 <- colnames(pa_t1)[3:ncol(pa_t1)]
    species_t2 <- colnames(pa_t2)[3:ncol(pa_t2)]
    
    # Identifies new species that appear in t2 but not t1
    new_species <- species_t2[which(!species_t2 %in% species_t1)]
    
    # Identifies extinct species (those with zero occurrences)
    dead_species_t1 <- colnames(pa_t1)[colSums(pa_t1)==0]  
    dead_species_t2 <- colnames(pa_t2)[colSums(pa_t2)==0] 

    # Identifies species that exist in t1 and are dead in t2
    dead_t2 <- dead_species_t2[dead_species_t2 %in% colnames(pa_t1)]
    
    # New extinctions - species that were alive in t1 but extinct in t2
    new_dead_species <- dead_t2[which(!dead_t2 %in% dead_species_t1)]
    
    #####_________EXTINCTION____________________________________________________________________________###
    # Where species alive at t1 extinct at t2, extinction column sums these events. Coordinates of 
    # these events are found by looking at location of the extinct species at t1
    
    # If extinction occurs
    if(length(new_dead_species) > 0){
      # Gets coordinates from t1 of newly extinct species
      extinction_result <- tryCatch({
        dat <- pa_t1[, c("x", "y", new_dead_species)]
        data <- as.data.frame(dat)
        
        # Calculates extinction events at each location
        ext <- data %>% mutate(extinction = rowSums(.[c(3:ncol(dat))]))
        ext_f <- ext %>% relocate(extinction, .after = y)
        extinction <- ext_f[c(0:3)] 
        
        extinction  # Return value
      }, error = function(e) {
        # If error occurs, create empty dataset
        data <- as.data.frame(pa_t1)
        ext <- data[c(0:2)]  
        ext$extinction <- NA
        ext  # Return value
      })
      
      # Save extinction data
      write_csv(extinction_result, file=file.path(out_dir, paste0("/extinction", i-1, ".csv")))
    } else { 
      # If no extinction occurs, creates empty dataset with NA values
      data <- as.data.frame(pa_t1)
      ext <- data[c(0:2)] 
      ext$extinction <- NA
      
      # Save extinction data
      write_csv(ext, file=file.path(out_dir, paste0("/extinction", i-1, ".csv")))
    }
    
    #####_________SPECIATION____________________________________________________________________________###
    # Handle speciation events 
    
    if(length(new_species) > 0){
      # Gets coordinates of new species
      speciation_result <- tryCatch({
        dat2 <- pa_t2[, c("x", "y", new_species)]
        data2 <- as.data.frame(dat2)
        
        # Calculates speciation events at each location
        spe <- data2 %>% mutate(speciation = rowSums(.[c(3:ncol(dat2))]))
        spe_f <- spe %>% relocate(speciation, .after = y)
        speciation <- spe_f[c(0:3)]  
        
        speciation  # Return value
      }, error = function(e) {
        # If error occurs, create empty dataset
        data <- as.data.frame(pa_t2)
        spe <- data[c(0:2)]  
        spe$speciation <- NA
        spe  # Return value
      })
      
      # Save speciation data
      write_csv(speciation_result, file=file.path(out_dir, paste0("/speciation", i-1, ".csv")))
    } else {
      # If no speciation occurs, creates empty dataset
      data <- as.data.frame(pa_t2)
      spe <- data[c(0:2)] 
      spe$speciation <- NA
      
      # Save speciation data
      write_csv(spe, file=file.path(out_dir, paste0("/speciation", i-1, ".csv")))
    }
  }
}

main <- function() {
  dir <- "./outputs/"
  args <- commandArgs(trailingOnly = TRUE) 
  
  # Check if arguments were provided
  if(length(args) < 1) {
    stop("Error: Please provide a path to the input directory as an argument.")
  }
  
  bio_tables(path = args[1], out = file.path(dir))
}

main()