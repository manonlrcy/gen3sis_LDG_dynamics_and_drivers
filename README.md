# Data & Scripts for 'Deep time evolution of the Latitudinal Diversity Gradient: insights from mechanistic models'
Manon Lorcery, Laurent Husson, Tristan Salles, Sébastien Lavergne, Oskar Hagen, Alexander Skeels. Deep time evolution of the Latitudinal Diversity Gradient: insights from mechanistic models.

## Data
Parameter table

```
params_table.txt
```

## Scripts 
### Gen3sis Configuration:

```
M0_batch.R
M1s_batch.R
M1d_batch.R
Mie_batch.R
```

### Create Gen3sis Landscape:

```
landscape_input.R
```

### Create Speciation, Extinction & alpha Richness tables based on Presence & Abscence Matrices (PA matrices)
This R script processes sequences of presence-absence matrices of species across geographic locations over multiple time steps (from t=150 to t=1). For each consecutive pair of time steps, it calculates and exports three key biodiversity metrics: species richness (the total number of species at each location), extinction events (where species present in the earlier time step are absent in the later one), and speciation events (where new species appear in the later time step). The script identifies new species and extinct species between time points, computes these metrics for each geographic coordinate, and saves the results as CSV files in an organized directory structure. It's designed to track spatial patterns of biodiversity change through time, likely for ecological or evolutionary studies.

```
bio_tables.R
```

Feel free to ask any questions: lorcerymanon@gmail.com
