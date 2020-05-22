# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact.
# Updated 2020.05.15 CH

## SETUP -----------------------------------------------------------------------

# Clear all from workspace/environment.
rm(list = ls())

# Install and load required packages.
# Install and load lattice for plotting results.
if(!require(lattice)){install.packages("lattice")}
library(lattice)
# Install and load grDevices for plot color schemes.
if(!require(grDevices)){install.packages("grDevices")}
library(grDevices)

## SET MODEL DIMENSIONS --------------------------------------------------------

# Define size of model domain. 
x.size <- 60000  # [m] 
z.size <- 20000  # [m] 

# Set resolution (number of nodes).  
x.num <- 300
z.num <- 200

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

T.liquidus.C <- 1195  # liquidus temperature, [°C]
# Suggested value: 1195°C from MELTS modeling
T.liquidus <- T.liquidus.C + 273.15  # liquidus temperature, [K]

T.solidus.C <- 796  # solidus temperature, [°C]
# Suggested value: 796°C from MELTS modeling
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

# Choose our own dt (has to be bigger than dt.min).
dt.choice.yr <- 100  # [yr]
dt <- dt.choice.yr / 3.17098e-8  # [s]

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
meltsheetvolume <- 800  # [km^3]
# Option 1: 800 km^3
# Option 2: 1000 km^3
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
T.meltsheet.C <- 2370  # temperature of melt sheet, [°C]
T.meltsheet <- T.meltsheet.C + 273.15  # temperature of melt sheet, [°C]
# Option 1: 1700°C; Abramov and Kring (2007)
# Option 2: 2000°C
# Option 3: 2370°C; Timms et al. (2017)


# Elevate the temperature at central uplift.

# # Option 1: raise geotherm everywhere beneath the melt sheet.
# # Set the temperature with which to raise the geotherm everywhere in the central
# # uplift. 
# # The height of uplift is dependent on the size of the crater. 40 km of uplift
# # is an estimate based on the size of the crater associated with a 30 km in 
# # diameter melt sheet.
# height.of.uplift <- 20  # amount of uplift, [km]
# # Option 1: 20 km
# # Option 2: 40 km
# T.uplift <- height.of.uplift * a.Kkm
# # Simple uplift of the geotherm temperature gradient.
# for (j in xstart.melt.index:xend.melt.index) {
#   for (i in 2:z.num){
#     T.0[j, i] <- T.0[j, i] + T.uplift
#   }
# }


# Option 2: raise geotherm in half ellipse below melt sheet, 
# Set melt sheet base depth in m as a variable to place uplift below melt sheet.
base.meltsheet.m <- meltsheetthickness + dz  # [m], add dz for air layer at top
# Set maximum depth of uplift region for ellipse calculations.
uplift.depth.m <- meltsheetthickness * 3

# Scale the temperature by position from center of ellipse.
# Calculate the position of half an ellipse of uplift under the melt sheet.
for (i in 1:x.num) {
  for (j in 1:z.num) {
    # Calculate the value at the node to determine if it inside or outside the 
    # ellipse.
    # Option 1: ellipse with the same diameter as the melt sheet:
    # node.value <- (((x[[i]] - mean(x))^2) / ((meltsheetdiameter / 2)^2)) + 
    #   (((z[[j]] - base.meltsheet.m)^2) / (uplift.depth.m^2))
    # Option 2: circle with the same diameter as the melt sheet:
    # node.value <- (((x[[i]] - mean(x))^2) / ((meltsheetdiameter / 2)^2)) + 
    #   (((z[[j]] - base.meltsheet.m)^2) / ((meltsheetdiameter / 2) ^ 2))
    # Option 3: Large circle with diameter 2 times larger than melt sheet 
    # diameter, centered above melt sheet:
    node.value <- (((x[[i]] - mean(x))^2) / ((meltsheetdiameter)^2)) + 
      (((z[[j]] - (-meltsheetdiameter / 2))^2) / ((meltsheetdiameter) ^ 2))
    # Assign a value to the uplift matrix that depends on whether the node is
    # inside out outside the uplift ellipse.
    if (node.value <= 1) {
      # Inside the ellipse, scale the temperature position to decrease from 
      # center of ellipse. 
      T.0[i, j] <- 1200 + -node.value * 600 
    } else {
      # Outside the ellipse, keep the geotherm the same temperature.
      T.0[i, j] <- T.0[i, j] 
    }
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

# Plot initial conditions to verify that the model is set up correctly.
levelplot(T.0 - 273.15,
          # Change x and y axis scales and labeling.
          scales = list(x = list(labels = seq(from = 0, 
                                              to = (x.size / 1000), 
                                              by = 10)),
                        y = list(labels = seq(from = 0, 
                                              to = (z.size / 1000), 
                                              by = 2))),
          # Flip y axis so depth increases down.
          # Cut off the bottom part of the model because of edge effects.
          ylim = c(80, 0),
          # Plot labels.
          xlab = "distance (km)", # x axis label
          ylab = "depth (km)", # y axis label
          # Make plot title number of years of model run.
          main = "t = 0",
          # Change aspect ratio of plot
          aspect = (z.size / 2) / x.size,
          # Change the color scheme.
          col.regions = rev(hcl.colors(n = 25, palette = 'Reds 2')),
          # Set colorbar limits and density of contour interval on plot.
          at = c(seq(0, 2500, length.out = 25))) 

# Plot initial conditions to verify that the model is set up correctly with a 
# more restricted T range.
levelplot(T.0 - 273.15,
          # Change x and y axis scales and labeling.
          scales = list(x = list(labels = seq(from = 0, 
                                              to = (x.size / 1000), 
                                              by = 10)),
                        y = list(labels = seq(from = 0, 
                                              to = (z.size / 1000), 
                                              by = 2))),
          # Flip y axis so depth increases down.
          # Cut off the bottom part of the model because of edge effects.
          ylim = c(80, 0),
          # Plot labels.
          xlab = "distance (km)", # x axis label
          ylab = "depth (km)", # y axis label
          # Make plot title number of years of model run.
          main = "t = 0",
          # Change aspect ratio of plot
          aspect = (z.size / 2) / x.size,
          # Change the color scheme.
          col.regions = rev(hcl.colors(n = 25, palette = 'Reds 2')),
          # Set colorbar limits and density of contour interval on plot.
          at = c(seq(0, 1000, length.out = 10))) 



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
n <- 1000  # Set number of time steps.
model.results <- RunModel(n)  # Run model for n timesteps.
T.n <- model.results$T.n  # Save the model results of temperature, [K]

# # Save model results.
# saveRDS(T.n, file = "/Users/claireharrigan/Dropbox/IGL + Research/Other projects/MeltSheet/Results/MeltSheet Results_2020.05.19/T.n_2370C_800km3_40km.RData")
# # # Load in previous model run results.
# # T.n.1700 <- list()
# T.n <- readRDS("/Users/claireharrigan/Dropbox/IGL + Research/Other projects/MeltSheet/Results/MeltSheet Results_2020.05.19/T.n_1700C_800km3_40km.RData")

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
  
  # Specify dataset by the model iteration number.
  InputMatrix <- T.n[[IterationNumber]]
  # Make plot.                                                
  OutputPlot <- levelplot(InputMatrix - 273.15,
    # Change x and y axis scales and labeling.
    scales = list(x = list(labels = seq(from = 0, to = (xval / 1000), by = 10)),
                  y = list(labels = seq(from = 0, to = (zval / 1000), by = 2))),
    # Flip y axis so depth increases down.
    # Cut off the bottom part of the model because of edge effects.
    ylim = c(80, 0),
    # Plot labels.
    xlab = "distance (km)", # x axis label
    ylab = "depth (km)", # y axis label
    # Make plot title number of years of model run.
    main = paste("t =", 
                 (round(((IterationNumber - 1) * dt.choice.yr), 
                       digits = 0)) / 1000, 
                 "ka",
                 sep = " "),
    # Change aspect ratio of plot
    aspect = (zval/2)/xval,
    # Change the color scheme.
    col.regions = rev(hcl.colors(n = 25, palette = 'Reds 2')),
    # col.regions = rev(hcl.colors(n = 25, palette = 'Grays')),  # Plot in grayscale.
    # contour = T,  # Add contour lines.
    # Set colorbar limits and density of contour interval on plot.
    at = c(seq(0, 2500, length.out = 25))) 
  return(OutputPlot)
}

# Plot results.
levelplotCH(1)
levelplotCH(5)
levelplotCH(26)
levelplotCH(61)
levelplotCH(101)
levelplotCH(501)
levelplotCH(1001)

# Output plots as a PDF.
# Page 1 (0, 50, 100, 200)
pdf("/Users/claireharrigan/Dropbox/IGL + Research/Other projects/MeltSheet/Results/MeltSheet Results_2020.05.21/MeltSheet_2Dplot_Morokweng_2020.05.21_01.pdf",
  width = 8.5, height = 11,
  useDingbats = FALSE)
print(levelplotCH(1), split = c(1,1,1,4), more = T) 
print(levelplotCH(51), split = c(1,2,1,4), more = T)
print(levelplotCH(101), split = c(1,3,1,4), more = T)
print(levelplotCH(201), split = c(1,4,1,4))
dev.off()

# Page 2 (300, 400, 500, 600)
pdf("/Users/claireharrigan/Dropbox/IGL + Research/Other projects/MeltSheet/Results/MeltSheet Results_2020.05.21/MeltSheet_2Dplot_Morokweng_2020.05.21_02.pdf",
    width = 8.5, height = 11,
    useDingbats = FALSE)
print(levelplotCH(301), split = c(1,1,1,4), more = T) 
print(levelplotCH(401), split = c(1,2,1,4), more = T)
print(levelplotCH(501), split = c(1,3,1,4), more = T)
print(levelplotCH(601), split = c(1,4,1,4))
dev.off()

# Page 3 (700, 800, 900, 1000)
pdf("/Users/claireharrigan/Dropbox/IGL + Research/Other projects/MeltSheet/Results/MeltSheet Results_2020.05.21/MeltSheet_2Dplot_Morokweng_2020.05.21_03.pdf",
    width = 8.5, height = 11,
    useDingbats = FALSE)
print(levelplotCH(701), split = c(1,1,1,4), more = T) 
print(levelplotCH(801), split = c(1,2,1,4), more = T)
print(levelplotCH(901), split = c(1,3,1,4), more = T)
print(levelplotCH(1001), split = c(1,4,1,4))
dev.off()

## PLOT 1D MODEL RESULTS -------------------------------------------------------

# Extract 1D results through center of 2D model results.

# Initialize list to hold center values from T.n matricies.
T.1D <- list()
# Iterate through T.n list, extract T values through center of melt sheet, save
# in list T.1D.
for (i in 1:length(T.n)) {
  T.1D[[i]] <- T.n[[i]][(x.num / 2), ]
}
# Iterate through T.n list, extract T values through 70.3% away from the center 
# of melt sheet (location of borehole M3), save in list T.1D.
# Find the node location of M3.
# M3 <- round((mean(x) - ((meltsheetdiameter / 2) * 0.703)) / dx, digits = 0)
# # Iterate through all matrices in list
# for (i in 1:length(T.n)) {
#   T.1D[[i]] <- T.n[[i]][(M3), ]
# }




# Extract limited results to plot.
T.1D.lim <- list()
T.1D.lim[[1]] <- T.1D[[1]]
T.1D.lim[[2]] <- T.1D[[26]]
T.1D.lim[[3]] <- T.1D[[51]]
T.1D.lim[[4]] <- T.1D[[76]]
T.1D.lim[[5]] <- T.1D[[101]]
T.1D.lim[[6]] <- T.1D[[201]]
T.1D.lim[[7]] <- T.1D[[301]]
T.1D.lim[[8]] <- T.1D[[401]]
T.1D.lim[[9]] <- T.1D[[501]]
T.1D.lim[[10]] <- T.1D[[601]]
T.1D.lim[[11]] <- T.1D[[701]]


# Plot 1D results using base R. 

# Set melt sheet base depth in km as a variable.
base.meltsheet <- (meltsheetthickness + dz) / 1000 # add dz for air layer at top

# Make plot.
plot.new()  # Start new plot.
# Make plot.
plot(range(0, 2500), range(0, 4),
     type = "n",
     xlab = "Temperature (°C)",
     ylab = "Depth (km)",
     ylim = c(4, 0),
     yaxs = "i",
     xlim = c(0, 2500),
     xaxs = "i")
# Set number of colors based on number of model results being plotted.
colors <- (hcl.colors(n = length(T.1D.lim), palette = "Spectral"))
# Add lines.
for (i in 1:length(T.1D.lim)) {
  tempprofile <- T.1D.lim[[i]] - 273.15
  lines(tempprofile, z / 1000,
        type = "l",
        lwd = 1.75,
        col = colors[i],
        lty = "solid")
}
# Add title.
title("1D model, center of melt sheet\nT.meltsheet = 2370°C, V.meltsheet = 800 km^3, \nuplift geometry: large circle centered above model domain",
      cex.main = 0.8)
# Add vertical lunes for liquidus and solidus.
lines(c(T.liquidus.C, T.liquidus.C), c(0, (z.size / 1000)),
      lwd = 2,
      col = "gray")
lines(c(T.solidus.C, T.solidus.C), c(0, (z.size / 1000)),
      lwd = 2,
      col = "gray")
text(x = T.liquidus.C + 30, y = 3.5, labels = "liquidus", srt = 270)
text(x = T.solidus.C + 30, y = 3.5, labels = "solidus", srt = 270)
# Add horizontal lines for sample horizons.
# M3-1
lines(c(0, 2500), c(base.meltsheet - 0.7566, base.meltsheet - 0.7566),
      lwd = 1,
      col = "black")
# M3-2
lines(c(0, 2500), c(base.meltsheet - 0.4711, base.meltsheet - 0.4711),
      lwd = 1,
      col = "black")
# M3-3
lines(c(0, 2500), c(base.meltsheet - 0.2599, base.meltsheet - 0.2599),
      lwd = 1,
      col = "black")
# M3-4
lines(c(0, 2500), c(base.meltsheet - 0.1716, base.meltsheet - 0.1716),
      lwd = 1,
      col = "black")
# M3-6
lines(c(0, 2500), c(base.meltsheet - 0.1081, base.meltsheet - 0.1081),
      lwd = 1,
      col = "black")
# base of melt sheet
lines(c(0, 2500), c(base.meltsheet, base.meltsheet),
      lwd = 1,
      col = "blue")
# Create a vector of years associated with lines plotted.
yrsplotted <- vector()
modeliteration <- vector()
# yrsplotted <- c(1, 26, 51, 76, 101, 151, 201, 251, 301, 351, 401, 451, 501, 551, 601, 651, 701, 751, 801, 851, 901, 951, 1001)
# yrsplotted <- c(1, 26, 51, 76, 101, 151, 201, 251, 301, 351, 401, 451, 501, 751, 1001, 1251, 1501, 1751, 2001)
yrsplotted <- c(1, 26, 51, 76, 101, 201, 301, 401, 501, 601, 701)
modeliteration <- yrsplotted - 1
yrsplotted <- round(((yrsplotted - 1) * dt.choice.yr), digits = 0) / 1000
# Add a legend. 
legend("bottomright", 
       # paste("t = ", yrsplotted, " years, model run #", modeliteration, sep = ""), 
       paste(yrsplotted), 
       cex = 0.8, 
       col = colors, 
       lty = "solid",
       lwd = 1.75,
       title = "time elapsed (ka)")


## PLOT MODEL RESULTS THROUGH TIME FOR EACH SAMPLE -----------------------------

# Extract a temperature value for each time step (each matrix in the T.n list)
# at a given depth and distance from the center of the melt sheet (sample 
# location).

# Initialize matrix to hold temperature values from T.n matricies for all 
# samples.
T.M3 <- matrix(data = NA,
               nrow = 5, # a row for each sample  
               ncol = n + 1)  # a column for each time step
# Rename T.M3 rows by sample.
rownames(T.M3) <- c("M3-1", "M3-2", "M3-3", "M3-4", "M3-6")
# Rename T.M3 columns by time step (in ka).
colnames(T.M3) <- paste((round(((seq(from = 0, to = n)) * dt.choice.yr), 
                               digits = 0)) / 1000, "ka", sep = " ")

# Find the node location (x axis) of M3 borehole (70.3% away from the center of
#   the meltsheet).
M3 <- round((mean(x) - ((meltsheetdiameter / 2) * 0.703)) / dx, digits = 0)

# Determine which row of the model (z axis) results matrices corresponds to each
# sample depth.
# Set melt sheet base depth in km as a variable.
base.meltsheet <- (meltsheetthickness + dz) / 1000 # add dz for air layer at top
# Initialize vector to hold sample depths.
z.samples <- vector()
# Find nodal depth of M3-1.
z.samples[[1]] <- round(((base.meltsheet - 0.7566) * 1000) / dz, digits = 0)
# M3-2
z.samples[[2]] <- round(((base.meltsheet - 0.4711) * 1000) / dz, digits = 0)
# M3-3
z.samples[[3]] <- round(((base.meltsheet - 0.2599) * 1000) / dz, digits = 0)
# M3-4
z.samples[[4]] <- round(((base.meltsheet - 0.1716) * 1000) / dz, digits = 0)
# M3-6
z.samples[[5]] <- round(((base.meltsheet - 0.1081) * 1000) / dz, digits = 0)


# Iterate through T.n list, extract T values at the M3 borehole for each sample,
# save in matrix T.M3.
for (n in 1:length(z.samples)) {
  for (i in 1:length(T.n)) {
    T.M3[n, i] <- T.n[[i]][z.samples[[n]], M3]
  }
}

# Plot T-t path for each sample and indicate when each passes through the
solidus.




## OUTPUT RESULTS --------------------------------------------------------------

# Make table of model parameters to output as .csv file.

# Create and populate a data frame of model parameters (MP).
MP <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
colnames(MP) <- c("Variable", "Units", "Value")

# Model domain and resolution.
MP[1, ]  <- c("Model domain and resolution", "", "")
MP[2, ]  <- c("Modeled domain, distance", "km", x.size / 1000)
MP[3, ]  <- c("Step size, distance", "m", round(dx, digits = 1))
MP[4, ]  <- c("Modeled domain, depth", "km", z.size / 1000)
MP[5, ]  <- c("Step size, depth", "m", round(dz, digits = 1))
MP[6, ]  <- c("Time step", "years", round(dt.choice.yr, digits = 1))

# Melt sheet size.
MP[7, ]  <- c("Melt sheet dimensions", "", "")
MP[8, ]  <- c("Melt sheet diameter", "km", meltsheetdiameter / 1000)
MP[9, ]  <- c("Melt sheet volume", "km^3", meltsheetvolume)
MP[10, ] <- c("Melt sheet thickness", "km", round(meltsheetthickness / 1000, digits = 2))
MP[11, ] <- c("Depth excavated during impact (central uplift)", "km", height.of.uplift)

# Temperatures.
MP[12, ] <- c("Temperatures", "", "")
MP[13, ] <- c("Melt sheet temperature at emplacement", "°C", T.meltsheet.C)
MP[14, ] <- c("Surface temperature", "°C", T.surface.C)
MP[15, ] <- c("Temperature, liquidus", "°C", T.liquidus.C)
MP[16, ] <- c("Temperature, solidus", "°C", T.solidus.C)

# Melt sheet and basement rock properties.
MP[17, ] <- c("Physical properties of basement and melt sheet", "", "")
MP[18, ] <- c("Geothermal gradient", "°C/km", a.Kkm)
MP[19, ] <- c("Latent heat of fusion of diopside", "J/kg", L)
MP[20, ] <- c("Heat capacity, no latent heat effects", "J/(kg*K)", Cp.nolatent)
MP[21, ] <- c("Heat capacity, accounting for latent heat of fusion", "J/(kg*K)", round(Cp.prime[1, 1], digits = 1))
MP[22, ] <- c("Thermal conductivity", "W/(m*K)", k[1, 1])
MP[23, ] <- c("Heat production, basement and melt sheet", "W/(m^3)", A[1, 1])
MP[24, ] <- c("Density, basement and melt sheet", "kg/(m^3)", rho[1,1])
MP[25, ] <- c("Thermal diffusivity, basement and melt", "(m^2)/s", round(Kappa[1, 1], digits = 7))
MP[26, ] <- c("Regional tectonic uplift", "m/s", U[1, 1])

# Write table to csv.
write.csv(MP, file = "/Users/claireharrigan/Dropbox/IGL + Research/Other projects/MeltSheet/Results/ModelParameters_2020.03.24.csv")








