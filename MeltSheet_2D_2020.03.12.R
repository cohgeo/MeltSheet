# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact.
# Updated 2020.03.12 CH

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

# Define size of model domain. 
x.size <- 60000  # [m] 
z.size <- 10000  # [m] 

# Set resolution (number of nodes).  
x.num <- 300
z.num <- 100

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
# Suggested value: 2.2e-6 is average A value from Jones (1988), table 6

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
a.Kkm <- 12  # [K/km]
# Used low end of gradient from Jones (1992) (12 and 21 K/km)
a <- a.Kkm / 1000  # [K/m]

# Define time step.
dt.max.x <- (dx ^ 2) / (6 * Kappa.min)  # [s] 
dt.max.z <- (dz ^ 2) / (6 * Kappa.min)  # [s] 
dt.min <- min(dt.max.x, dt.max.z)
dt.yr <- dt.min * 3.17098e-8 # [yr]
dt <- dt.min


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
T.0[ , 1] <- T.surface + 273.15 # [°C + 273.15 to convert to K] 

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
meltsheetthickness <- meltsheetvolume.m3 / 
                      (((meltsheetdiameter / 2)^2) * 3.14)  # [m]
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
zend.melt.index <- round(zend.melt / dz, digits = 0)

# Set temperature of melt sheet.
T.meltsheet <- 1700 + 273.15 # [°C + 273.15 to convert to K] 
  # Suggested value: 1700°C; Abramov and Kring (2007)

# Elevate geotherm at central uplift.
# Set the temperature at the center and outer edge of the uplift.
T.uplift <- 900 + 273.15 # [°C + 273.15 to convert to K] 
T.edge <- 500 + 273.15 # [°C + 273.15 to convert to K] 
# Set the diameter of the central uplift to be 10% larger than the melt sheet 
# diameter on either side.
uplift.edge.low <- xstart.melt.index - 0.1 * (meltsheetdiameter / 1000)
uplift.edge.high <- xend.melt.index + 0.1 * (meltsheetdiameter / 1000)

# # Optional simple uplift temperature gradient.
# for (j in uplift.edge.low:uplift.edge.high) {
#   for (i in 2:z.num){
#     T.0[j, i] <- T.uplift
#   }
# }

# Optional simple uplift of the geotherm temperature gradient.
for (j in uplift.edge.low:uplift.edge.high) {
  for (i in 2:z.num){
    T.0[j, i] <- T.0[j, i] + 480
  }
}


# # For more complex uplift temperature gradient:
# # Set up quadratic.
# QA <- (T.edge + T.uplift) / (uplift.edge.low^2 + (median(x) / 1000)^2 + 
#       (((uplift.edge.low^2 - uplift.edge.high^2) / 
#       (uplift.edge.high - uplift.edge.low)) * (uplift.edge.low + 
#       (median(x) / 1000))))
# QB <- QA * (uplift.edge.low^2 - uplift.edge.high^2) / 
#       (uplift.edge.high - uplift.edge.low)
# QC <- (QA * (median(x) / 1000)^2) + ((QA * (uplift.edge.low^2 - 
#       uplift.edge.high^2) * (median(x) / 1000)) / (uplift.edge.high - 
#       uplift.edge.low)) - (T.uplift)
# 
# # Calculate temperatures for central uplift.
# for (j in uplift.edge.low:uplift.edge.high) {
#   for (i in 2:z.num) {
#     T.0[j, ] <- (QA * (x[j] / 1000)^2) + (QB * (x[j] / 1000)) + QC
#   }
# }

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
rm(i, j, QA, QB, QC)

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
        T.n[[t + 1]][, 1] <- T.n[[1]][1, 1]          # top
        T.n[[t + 1]][, z.num] <- T.n[[1]][x.num, z.num]  # bottom
        T.n[[t + 1]][x.num, ] <- T.n[[t]][x.num, ]  # right
        T.n[[t + 1]][1, ] <- T.n[[1]][1, ]          # left
        
        
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
    
    # Set constant Cp border.
    Cp[[t + 1]][1, ] <- Cp.initial[1, ] # left
    Cp[[t + 1]][x.num, ] <- Cp.initial[x.num, ] # right
    Cp[[t + 1]][, 1] <- Cp.initial[, 1] # top
    Cp[[t + 1]][, z.num] <- Cp.initial[, z.num] # bottom
    
    for (i in 2:(x.num - 1)) {
      for (j in 2:(z.num - 1)) {
        
        # Calculate Cp for interior
        Cp[[t + 1]][i, j] <- ifelse(T.n[[t + 1]][i, j] < T.liquidus && T.n[[t + 1]][i, j] > T.solidus, Cp.prime[i, j], Cp.initial[i, j])
      }
    }
  }
  objects.to.return <- list("T.n" = T.n, "Cp" = Cp)
  return(objects.to.return)
}


## RUN MODEL -------------------------------------------------------------------

# Solve the finite difference equation for n time steps.
n <- 10
model.results <- RunModel(n)
T.n <- model.results$T.n
Cp <- model.results$Cp

# 6 min to do 1000 model runs

## PLOT MODEL RESULTS ----------------------------------------------------------

# Plot model results.
levelplot(T.n[[1]])
levelplot(T.n[[25]])
levelplot(T.n[[50]])
levelplot(T.n[[75]])
levelplot(T.n[[100]])
levelplot(T.n[[200]])
levelplot(T.n[[300]])
levelplot(T.n[[400]])
levelplot(T.n[[500]])
levelplot(T.n[[600]])
levelplot(T.n[[700]])
levelplot(T.n[[800]])
levelplot(T.n[[900]])
levelplot(T.n[[1000]])


# Make plotting function.
levelplotCH <- function(InputMatrix, 
                        plot.title = "",
                        xnum = x.num,
                        znum = z.num,
                        xval = x.size,
                        zval = z.size) {
  
  OutputPlot <- levelplot(InputMatrix - 273.15,
                          
                          # Change axes.
                          # xlim = c(0, xsize),
                          ylim = c(znum, 0),
                          # row.values = zval,
                          # column.values = xval,
                          # 
                          # Plot labels.
                          xlab = "distance (number of nodes)", # x axis label
                          ylab = "depth (numbner of nodes)", # y axis label
                          # xlab = "distance", # x axis label
                          # ylab = "depth", # y axis label
                          main = plot.title, # plot title
                          
                          # Change the look of the plot.
                          # aspect = 1, # change aspect ratio of plot
                          contour = T, # add concour lines
                          
                          # Adjust colorbar.
                          # colorkey = c(400, 1200), # set colorbar limits
                          at = c(seq(0, 1800, length.out = 50)),
                          cuts = 10) # change the number of colors in colorbar
  return(OutputPlot)
}

# Testing plotting function.
levelplotCH(T.0)
levelplotCH(T.n[[1]])
levelplotCH(T.n[[5]])
levelplotCH(T.n[[25]])
levelplotCH(T.n[[50]])
levelplotCH(T.n[[75]])
levelplotCH(T.n[[100]])

# Possible function inputs.
# plot.title = paste((round(dt.yr * 5, digits = 0)), " yrs")

z ~ x * y