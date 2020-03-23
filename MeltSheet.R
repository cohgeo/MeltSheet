# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact.
# Updated 2020.03.18 CH

## SETUP -----------------------------------------------------------------------

# Clear all from workspace/environment.
rm(list = ls())

# Install and load required packages.
# Install and load lattice for plotting results.
if(!require(lattice)){install.packages("lattice")}
library(lattice)
# Install and load RColorBrewer for plot color schemes.
if(!require(RColorBrewer)){install.packages("RColorBrewer")}
library(RColorBrewer)


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
# All temperatures in °C should be converted to K.

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

# Set the top row of the model to be the air temperature.
# Choose a value for the air temperature.
T.surface.C <- 25  # surface temperature, [°C]
T.surface <- T.surface.C + 273.15  # surface temperature, [K]
# Assign the first row of the initial matrix to be the air temperature.
T.0[ , 1] <- T.surface  # [K] 

# Make a geotherm for the initialization matrix.
for (i in 2:dim(T.0)[2]) {
  T.0[ , i] <- T.0[1, 1] + a * x[i]
}

# Instantanteously emplace a melt sheet and set material properties of the melt
# sheet.
# Set dimensions of melt sheet.
meltsheetdiameter <- 30000  # [m]
meltsheetvolume <- 1000  # [km^3]
meltsheetvolume.m3 <- meltsheetvolume * (1000^3)
meltsheetthickness <- meltsheetvolume.m3 / 
                      (((meltsheetdiameter / 2)^2) * pi)  # [m]
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
T.meltsheet.C <- 1700  # temperature of melt sheet, [°C]
T.meltsheet <- T.meltsheet.C + 273.15  # temperature of melt sheet, [°C]
  # Suggested value: 1700°C; Abramov and Kring (2007)

# Elevate geotherm at central uplift.
# Set the temperature with which to raise the geotherm everywhere in the central
# uplift. 
# The height of uplift is dependent on the size of the crater. 40 km of uplift
# is an estimate based on the size of the crater associated with a 30 km in 
# diameter melt sheet.
height.of.uplift <- 40  # amount of uplift, [km]
T.uplift <- height.of.uplift * a.Kkm

# Optional simple uplift of the geotherm temperature gradient.
for (j in xstart.melt.index:xend.melt.index) {
  for (i in 2:z.num){
    T.0[j, i] <- T.0[j, i] + T.uplift
  }
}

# Emplace melt sheet. Set material properties of melt sheet.
for (i in xstart.melt.index:xend.melt.index) {
  for (j in zstart.melt.index:zend.melt.index) {
    T.0[i, j] <- T.meltsheet
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
n <- 5  # Set number of time steps.
model.results <- RunModel(n)  # Run model for n timesteps.
T.n <- model.results$T.n  # Save the model results of temperature, [K]


## PLOT 2D MODEL RESULTS -------------------------------------------------------

# Make plotting function to plot model results in °C using the lattice package.
levelplotCH <- function(IterationNumber,
                        x.axis.in.m = x,
                        z.axis.in.m = z,
                        dtyr = dt.yr,
                        xnum = x.num,
                        znum = z.num,
                        xval = x.size,
                        zval = z.size) {
  
  InputMatrix <- T.n[[IterationNumber]]
                                                  
  OutputPlot <- levelplot(InputMatrix - 273.15,
                         
                          # Change x and y axis scales and labeling.
                          # Change axis labels.
                          scales = list(x = list(at = seq(from = 0, to = xval, by = 50),
                                                 labels = seq(from = 0, to = (xval / 1000), by = 10)),
                                        y = list(at = seq(from = 0, to = zval, by = 20),
                                                 labels = seq(from = 0, to = (zval / 1000), by = 2))),
                          # Flip y axis so depth increases down.
                          # Cut off the bottom part of the model because of 
                          # edge effects.
                          ylim = c(80, 0),
                          
                          # Plot labels.
                          xlab = "distance (km)", # x axis label
                          ylab = "depth (km)", # y axis label
                          # Make plot title number of years of model run.
                          main = paste("t =", 
                                       round(((IterationNumber - 1) * dtyr), 
                                             digits = 0), 
                                       "years", 
                                       sep = " "),
                          
                          # Change the look of the plot.
                          aspect = zval/xval, # change aspect ratio of plot
                          # Change the color scheme.
                          col.regions = colorRampPalette(brewer.pal(9, 
                                                                    'Purples')),

                          # Adjust colorbar and contour lines.
                          # contour = T,  # Add contour lines.
                          # colorkey = list(labels = "T °C"),
                          # Set colorbar limits.
                          at = c(seq(0, 1800, 
                                     # Set density of contour interval on plot.
                                     length.out = 20))) 
  return(OutputPlot)
}

# Test plotting function.
levelplotCH(1)
levelplotCH(4)



# Plot results using base R plotting.

# Function for making a colorbar.
# Modified from Aurélien Madouasse, https://aurelienmadouasse.wordpress.com/2012/01/13/legend-for-a-continuous-color-scale-in-r/

legend.col <- function(col, lev){
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}


MeltSheetPlot <- function(IterationNumber) {
  # Make plot.
        # Input data into plot.
  image(x = x / 1000,  # Set x axis scale.
        y = z / 1000,  # Set y axis scale.
        z = (T.n[[IterationNumber]]) - 273.15,  # Set values to be plotted.
        
        # Set x and y axes and labels.
        xlim = c(0, (x.size / 1000)),
        ylim = c((z.size / 1000), 0),  # Make sure depth axis is reversed.
        xlab = "Distance (km)",
        ylab = "Depth (km)",
        
        # Set title.
        main = paste("t =", round(((IterationNumber - 1) * dt.yr), digits = 0), "years", sep = " "),
        
        # Set limits of color range.
        zlim = c(0, 1800),
        
        # Plot in grayscale.
        col = gray.colors(10, start = 0.9, end = 0.2))
  
  # Add colorbar.
  legend.col(col = gray.colors(10, start = 0.9, end = 0.2), 
             lev = c(0, 1800))
}


MeltSheetPlot(4)







## PLOT 1D MODEL RESULTS -------------------------------------------------------

# Extract 1D results through center of 2D model results.


# Plot 1D results.