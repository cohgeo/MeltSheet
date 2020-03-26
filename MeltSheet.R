# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact.
# Updated 2020.03.24 CH

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
n <- 500  # Set number of time steps.
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
  
  # Specify dataset by the model iteration number.
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
                                       "years, model iteration =",
                                       IterationNumber - 1,
                                       sep = " "),
                          
                          # Change the look of the plot.
                          aspect = zval/xval, # change aspect ratio of plot
                          # Change the color scheme.
                          col.regions = colorRampPalette(brewer.pal(9,
                                        'Purples')),
                                        # 'Greys')),  # Plot in grayscale.

                          # Adjust colorbar and contour lines.
                          # contour = T,  # Add contour lines.
                          # colorkey = list(labels = "T °C"),
                          # Set colorbar limits.
                          at = c(seq(0, 1800, 
                                     # Set density of contour interval on plot.
                                     length.out = 20))) 
  return(OutputPlot)
}

# Plot results.
levelplotCH(1)
levelplotCH(26)
levelplotCH(51)
levelplotCH(101)
levelplotCH(126)
levelplotCH(151)
levelplotCH(176)
levelplotCH(201)
levelplotCH(500)

# Output plots as a PDF.
# Page 1 (0, 25, 50, 100)
pdf("/Users/claireharrigan/Desktop/MeltSheetPlots/MeltSheet_2Dplots_0-25-50-100.pdf",
  width = 8.5, height = 11,
  useDingbats = FALSE)
print(levelplotCH(1), split = c(1,1,1,4), more = T) 
print(levelplotCH(26), split = c(1,2,1,4), more = T)
print(levelplotCH(51), split = c(1,3,1,4), more = T)
print(levelplotCH(101), split = c(1,4,1,4))
dev.off()

# Page 2 (125, 150, 175, 200)
pdf("/Users/claireharrigan/Desktop/MeltSheetPlots/MeltSheet_2Dplots_125-150-175-200.pdf",
    width = 8.5, height = 11,
    useDingbats = FALSE)
print(levelplotCH(126), split = c(1,1,1,4), more = T) 
print(levelplotCH(151), split = c(1,2,1,4), more = T)
print(levelplotCH(176), split = c(1,3,1,4), more = T)
print(levelplotCH(201), split = c(1,4,1,4))
dev.off()

# Page 3 (225, 250, 275, 300)
pdf("/Users/claireharrigan/Desktop/MeltSheetPlots/MeltSheet_2Dplots_225-250-275-300.pdf",
    width = 8.5, height = 11,
    useDingbats = FALSE)
print(levelplotCH(226), split = c(1,1,1,4), more = T) 
print(levelplotCH(251), split = c(1,2,1,4), more = T)
print(levelplotCH(276), split = c(1,3,1,4), more = T)
print(levelplotCH(301), split = c(1,4,1,4))
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

# Extract limited results to plot.
T.1D.lim <- list()
T.1D.lim[[1]] <- T.1D[[1]]
T.1D.lim[[2]] <- T.1D[[26]]
T.1D.lim[[3]] <- T.1D[[51]]
T.1D.lim[[4]] <- T.1D[[101]]
T.1D.lim[[5]] <- T.1D[[126]]
T.1D.lim[[6]] <- T.1D[[151]]
T.1D.lim[[7]] <- T.1D[[176]]
T.1D.lim[[8]] <- T.1D[[201]]
T.1D.lim[[9]] <- T.1D[[226]]
T.1D.lim[[10]] <- T.1D[[251]]
T.1D.lim[[11]] <- T.1D[[276]]
T.1D.lim[[12]] <- T.1D[[301]]
T.1D.lim[[13]] <- T.1D[[401]]
T.1D.lim[[14]] <- T.1D[[501]]


# Plot all 1D results using base R. 

# Set axis limits for plots.
xrange <- range(0, 1800)  # T range
yrange <- range(0, 10)  # depth in km

# Make plot.
plot.new()  # Start new plot.
plot(xrange, yrange,
     type = "n",
     xlab = "Temperature (°C)",
     ylab = "Depth (km)",
     ylim = c(10, 0),)
# Set number of colors based on number of model results being plotted.
colors <- rainbow(length(T.1D))  # Set number of colors.
# Add lines.
for (i in 1:length(T.1D)) {
  tempprofile <- T.1D[[i]] - 273.15
  lines(tempprofile, z / 1000,
        type = "l",
        lwd = 1.5,
        col = colors[i],
        lty = "solid")
}
# Add title.
title("1D model of temperature through center of the melt sheet")
# Add a legend. 
# legend(1600, 8, 1:length(T.1D), cex=0.8, col=colors, lty = "solid", title = "Model run")
# Add vertical lunes for liquidus and solidus.
lines(c(T.liquidus.C, T.liquidus.C), c(0, (z.size / 1000)))
lines(c(T.solidus.C, T.solidus.C), c(0, (z.size / 1000)))


# Plot key 1D results using base R. 

# Make plot.
plot.new()  # Start new plot.
# Make plot.
plot(xrange, yrange,
     type = "n",
     xlab = "Temperature (°C)",
     ylab = "Depth (km)",
     ylim = c(10, 0),
     xlim = c(0, 1800))
# Set number of colors based on number of model results being plotted.
colors <- rainbow(length(T.1D.lim))  # Set number of colors.
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
title("1D model of temperature through center of the melt sheet")
# Add vertical lunes for liquidus and solidus.
lines(c(T.liquidus.C, T.liquidus.C), c(0, (z.size / 1000)),
      lwd = 2,
      col = "gray")
lines(c(T.solidus.C, T.solidus.C), c(0, (z.size / 1000)),
      lwd = 2,
      col = "gray")
text(x = T.liquidus.C + 30, y = 3.5, labels = "liquidus", srt = 270)
text(x = T.solidus.C + 30, y = 3.5, labels = "solidus", srt = 270)
# Create a vector of years associated with lines plotted.
yrsplotted <- vector()
modeliteration <- vector()
yrsplotted <- c(1, 26, 51, 101, 126, 151, 176, 201, 226, 251, 276, 301, 401, 501)
modeliteration <- yrsplotted - 1
yrsplotted <- round(((yrsplotted - 1) * dt.yr), digits = 0)
# Add a legend. 
legend(900, 5.75, 
       paste("t = ", yrsplotted, " years, model run #", modeliteration, sep = ""), 
       cex = 0.8, 
       col = colors, 
       lty = "solid",
       lwd = 1.75,
       title = "Years elapsed")


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
MP[6, ]  <- c("Time step", "years", round(dt.yr, digits = 1))

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








