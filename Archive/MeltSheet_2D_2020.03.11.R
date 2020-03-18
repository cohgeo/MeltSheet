# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact.
# Updated 2020.03.11 CH

## SETUP -----------------------------------------------------------------------

# Clear all from workspace/environment.
rm(list = ls())

# Install required packages.
install.packages("lattice")
library(lattice)
# install.packages("RColorBrewer")
# library(RColorBrewer)
# install.packages("ggplot2")
# library(ggplot2)


## SET MODEL DIMENSIONS --------------------------------------------------------

# Define size of model domain. Code works for a square.
x.size <- 100000  # [m] 
z.size <- 50000  # [m] 

# Set resolution (number of nodes).  
x.num <- 100
z.num <- 50

# Define stepsize.
dx <- x.size / (x.num - 1)  # [m] 
dz <- z.size / (z.num - 1)  # [m]  

# Create vectors of node locations.
x <- seq(from = 0, to = x.size, by = dx)  # [m]
z <- seq(from = 0, to = z.size, by = dz)  # [m]


## SET CONSTANTS ---------------------------------------------------------------
# Define material properties, initialize arrays.

L <- 421000  # latent heat of fusion of diopside, [J/kg] 
# Suggested value: 421000 J/kg; Abramov and Kring (2007)

T.liquidus.C <- 1177  # liquidus temperature, [°C]
# Suggested value: 1177°C; Abramov and Kring (2007)
T.liquidus <- T.liquidus.C + 273.15  # liquidus temperature, [K]

T.solidus.C <- 997  # solidus temperature, [°C]
# Suggested value: 997°C; Abramov and Kring (2007)
T.solidus <- T.solidus.C + 273.15  # solidus temperature, [K]

Cp.nolatent <- 1000  # heat capacity with no latent heat, [J/(kg*K)]
# Suggested value: 1000 J/(kg/K), heat capacity of basement, melt, breccia; 
# Abramov and Kring (2007)

Cp.initial <- matrix(data = Cp.nolatent,   # initial heat capacity, [J/(kg*K)]
                     nrow = x.num,  
                     ncol = z.num)

Cp.prime <- Cp.initial + (L / (T.liquidus - T.solidus)) 
# heat capacity accounting for latent heat of fusion, [J/(kg*K)]
# Suggested value: 3339 J/(kg*K); Abramov and Kring (2007)
Cp.0 <- Cp.initial  # initialize matrix for initial model setup below

k <- matrix(data = 2.5,    # thermal conductivity, [W/(m*K)]
            nrow = x.num,  
            ncol = z.num)  
# Suggested value: 2.5 W/(m*K), thermal conductivity of basement, melt, 
# breccia; Abramov and Kring (2007)

A <- matrix(data = 2.2e-6,      # heat production, [W/(m^3)]
            nrow = x.num,  
            ncol = z.num)  
# Average A value from Jones (1988), table 6

rho <- matrix(data = 2700,   # density, [kg/(m^3)]
              nrow = x.num,  
              ncol = z.num)
# Suggested value: 2700 kg/(m^3), density of basement, melt, breccia; 
# Abramov and Kring (2007)

Kappa <- k / (rho * Cp.nolatent)  # thermal diffusivity, [(m^2)/s]
Kappa.min <- min(Kappa)      # [(m^2)/s] 

U <- matrix(data = 0,        # uplift, [m/s]
            nrow = x.num,  
            ncol = z.num)
# We assumed no tectonic uplift. 

# Input geothermal gradient in K/km.
a.Kkm <- 16.5  # [K/km]
# Average value between 12 and 21 K/km from Jones (1992)
a <- a.Kkm / 1000  # [K/m]

# Define time step.
dt.max <- (dx ^ 2) / (6 * Kappa.min)  # [s] 
dt.yr <- dt.max * 3.17098e-8 # [yr]
dt <- dt.max


## INITIAL T DISTRIBUTION ------------------------------------------------------
# Set up an initial condition of an instantaneously emplaced melt sheet at an 
# elevated temperature.

# Initialize a matrix for initial temperature conditions.
T.0 <- matrix(data = 0,
              nrow = x.num,  
              ncol = z.num)

# Set the top row to be air temperature.
# Choose a value for the air temperature.
T.surface <- 25  # [°C]
# Assign the first row of the initial matrix to be the air temperature.
T.0[ , 1] <- T.surface + 273.15

# Make a geotherm for the initialization matrix.
for (i in 2:dim(T.0)[2]) {
  T.0[ , i] <- T.0[1, 1] + a * x[i]
}

# Instantanteously emplace a melt sheet and set material properties of the melt
# sheet.
# Set dimensions of melt sheet.
meltsheetdiameter <- 30000  # [m]
meltsheetvolume <- 1000 # [km^3]
meltsheetvolume.m3 <- meltsheetvolume * (1000^3)
meltsheetthickness <- meltsheetvolume.m3 / (((meltsheetdiameter / 2)^2) * 3.14)  # [m]
# If the melt sheet is centered in the x direction in the model, find the 
# position in meters along the x axis.
xstart.melt <- mean(x) - (meltsheetdiameter / 2)  # x position, [m]
xend.melt <- mean(x) + (meltsheetdiameter / 2)  # x position, [m]
# Find the index in the vector that corresponds to the position above.
xstart.melt.index <- which(abs(x - xstart.melt) == min(abs(x - xstart.melt)))
xend.melt.index <- which(abs(x - xend.melt) == min(abs(x - xend.melt)))
# Make the melt sheet start on the second row of the matrix, below the air
# layer.
zstart.melt.index <- 2 
# If the melt sheet is meltsheetthickness and starts at the second row of the
# model, find the index that corresponds to the bottom of the melt sheet.
zend.melt <- z[2] + meltsheetthickness
zend.melt.index <- which(abs(z - zend.melt) == min(abs(z - zend.melt)))

# Set temperature of melt sheet.
T.meltsheet <- 1700  # [°C] 
  # Suggested value: 1700°C; Abramov and Kring (2007)

# Emplace half of central uplift temperature gradient. 
for (j in (xstart.melt.index - 0.1 * (meltsheetdiameter / 1000)):(median(x) / 1000)) {
  T.0[j, ] <- 900
  # T.0[j, ] <- T.0[1, 1] + a / 2 * z[j]
}

# Emplace other half of central uplift temperature gradient. Set material 
# properties of central uplift.
for (j in ((median(x) / 1000) + 1):(xend.melt.index + 0.1 * (meltsheetdiameter / 1000))) {
  T.0[j, ] <- 500
  # T.0[j, ] <- T.0[1, 1] + a / 2 * z[j]
}

# Emplace melt sheet. Set material properties of melt sheet.
for (i in xstart.melt.index:xend.melt.index) {
  for (j in zstart.melt.index:zend.melt.index) {
    T.0[i, j] <- T.meltsheet + 273.15
    # k[i, j]     <- 3.0    # thermal conductivity, [W/(m*K)]
    # A[i, j]     <- 1e-6   # heat production, [W/(m^3)]
    # rho[i, j]   <- 3300   # density, [kg/(m^3)]
    if (T.0[i, j] < T.liquidus & T.0[i, j] > T.solidus) {
      Cp.0[i, j] <- Cp.initial[i, j]   # heat capacity, [J/(kg*K)]
    } else {
      Cp.0[i, j] <- Cp.prime[i, j]  
    }
    # Kappa[i, j] <- k[i, j] / (rho[i, j] * Cp[i, j])  # [(m^2)/s]
  }
}

# Clean up extra variables.
rm(i, j)

## PLOT INITIAL CONDITIONS -----------------------------------------------------

# Plot the initial temperature matrix as a heatmap with the package lattice.
levelplot(T.0)


## CREATE FUNCTION OF THERMAL MODEL --------------------------------------------
# Make a function to solve the finite difference equation for n time steps.

RunModel <- function(n.timesteps) {
  
  # Initialize a list as a container for model returns at different times.
  T.n  <- list()
  
  # Place T.0 in the first position of the T.n list.
  T.n[[1]] <- T.0
  
  # Initialize a list as a container for model returns at different times.
  Cp  <- list()
  
  # Place T.0 in the first position of the T.n list.
  Cp[[1]] <- Cp.0
  
  for (t in 1:n.timesteps) {
    
    # Initialize T.n[[t + 1]]
    T.n[[t + 1]] <- matrix(data = NA,
                           nrow = x.num,
                           ncol = z.num)
    
    # Initialize Cp[[t + 1]]
    Cp[[t + 1]] <- matrix(data = NA,
                          nrow = x.num,
                          ncol = z.num)
    
    # Calculate T for interior of T.n[[t + 1]] matrix.
    for (i in 2:(x.num - 1)) {
      for (j in 2:(z.num - 1)) {
        
        # Set constant T border.
        T.n[[t + 1]][, 1] <- T.n[[t]][, 1]          # top
        T.n[[t + 1]][, z.num] <- T.n[[t]][, z.num]  # bottom
        T.n[[t + 1]][1, ] <- T.n[[t]][1, ]          # left
        T.n[[t + 1]][x.num, ] <- T.n[[t]][x.num, ]  # right
        
        # Calculate T in the interior of the model
        T.n[[t + 1]][i, j] <- T.n[[t]][i, j] +
          
          # conduction term
          ((Kappa[i, j] * dt) * (((T.n[[t]][(i + 1), j] - (2 * T.n[[t]][i, j]) + T.n[[t]][(i - 1), j]) / (dx ^ 2)) + ((T.n[[t]][i, (j + 1)] - (2 * T.n[[t]][i, j]) + T.n[[t]][i, (j - 1)]) / (dz ^ 2) ))) +
          
          # heat production term
          ((A[i, j] * dt) / (rho[i, j] * Cp[[t]][i, j])) +
          
          # advection/uplift term
          (U[i, j] * dt * ((T.n[[t]][i, (j + 1)] - T.n[[t]][i, j]) / dz))
      }
    }
    
    # Calculate Cp based on the most recent run of the model to be used in the
    # next iteration of the model
    for (i in 1:x.num) {
      for (j in 1:z.num) {
        
        Cp[[t + 1]][i, j] <- ifelse(T.n[[t + 1]][i, j] < T.liquidus && T.n[[t + 1]][i, j] > T.solidus, Cp.prime[i, j], Cp.initial[i, j])
        
      }
    }
  }
  objects.to.return <- list("T.n" = T.n, "Cp" = Cp)
  return(objects.to.return)
}


## RUN MODEL -------------------------------------------------------------------

# Solve the finite difference equation for n time steps.
n <- 150
model.results <- RunModel(n)
T.n <- model.results$T.n
Cp <- model.results$Cp


## PLOT MODEL RESULTS ----------------------------------------------------------

# Plot model results.
levelplot(T.0)
levelplot(T.n[[5]])
levelplot(T.n[[10]])
levelplot(T.n[[15]])
levelplot(T.n[[20]])
levelplot(T.n[[25]])
levelplot(T.n[[30]])
levelplot(T.n[[40]])
levelplot(T.n[[50]])
levelplot(T.n[[75]])
levelplot(T.n[[100]])
levelplot(T.n[[125]])
levelplot(T.n[[150]])
dt.yr*50


# Make plotting function.
levelplotCH <- function(InputMatrix) {
  OutputPlot <- levelplot(InputMatrix,
                          xlab = "distance",
                          ylab = "depth")
  return(OutputPlot)
}

# Testing plotting function.
levelplotCH(T.n[[30]])
