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

# Make plotting function to plot model results in °C.
levelplotCH <- function(InputMatrix, 
                        plot.title = "",
                        x.axis.in.m = x,
                        z.axis.in.m = z,
                        xnum = x.num,
                        znum = z.num,
                        xval = x.size,
                        zval = z.size) {
  
  # Rename columns and rows so scales of x and y axes make sense.
  # rownames(InputMatrix) <- x.axis.in.km
  # colnames(InputMatrix) <- z.axis.in.km 
  
  # Make plot.
  # OutputPlot <- levelplot(InputMatrix ~ (x.axis.in.m / 1000) * (z.axis.in.m / 1000),
                 
                          OutputPlot <- levelplot((InputMatrix) ~ (x.axis.in.m / 1000) * (z.axis.in.m / 1000),
                                                  
                          # OutputPlot <- levelplot(InputMatrix - 273.15,
                                                  
                          
                          
                          # For the functions documented here, the formula is generally of the form y ~ x | g1 * g2 * ... (or equivalently, y ~ x | g1 + g2 + ...), indicating that plots of y (on the y-axis) versus x (on the x-axis) should be produced conditional on the variables g1, g2, .... Here x and y are the primary variables, and g1, g2, ... are the conditioning variables. The conditioning variables may be omitted to give a formula of the form y ~ x, in which case the plot will consist of a single panel with the full dataset. The formula can also involve expressions, e.g., sqrt(), log(), etc. See the data argument below for rules regarding evaluation of the terms in the formula.
                          
                          # or the formula method, a formula of the form z ~ x * y | g1 * g2 * ..., where z is a numeric response, and x, y are numeric values evaluated on a rectangular grid. g1, g2, ... are optional conditional variables, and must be either factors or shingles if present.
                          
                          # Change axes.
                          # xlim = c(0, xsize),
                          # ylim = c((max(z)), 0),
                          # row.values = zval,
                          # column.values = xval,
                          # 
                          # Plot labels.
                          xlab = "distance (km)", # x axis label
                          ylab = "depth (km)", # y axis label
                          # xlab = "distance", # x axis label
                          # ylab = "depth", # y axis label
                          main = plot.title, # plot title
                          # key=list(title="Three Cylinder Options"),
                          
                          # Change the look of the plot.
                          # aspect = 1, # change aspect ratio of plot
                          contour = T, # add contour lines
                          # Change the color scheme.
                          col.regions = colorRampPalette(brewer.pal(9, 'Purples')),
                          # colorkey = list(labels = "T °C"),
                          # scales = list(log = "e"),
                          
                          # Adjust colorbar.
                          # Set colorbar limits
                          at = c(seq(0, 1800, length.out = 50))) 
  return(OutputPlot)
}

# Test plotting function.
levelplotCH(T.n[[1]])


levelplotCH(T.n[[5]])
levelplotCH(T.n[[25]])
levelplotCH(T.n[[50]])
levelplotCH(T.n[[75]])
levelplotCH(T.n[[100]])

# Possible function inputs.
# plot.title = paste((round(dt.yr * 5, digits = 0)), " yrs")



# Using base R plotting

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
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), 
         c("0",".5","1"), fill = gray.colors(10, start = 0.9, end = 0.2), xpd = NA)
}


MeltSheetPlot(1)




## PLOT 1D MODEL RESULTS -------------------------------------------------------

# Extract 1D results through center of 2D model results.


# Plot 1D results.