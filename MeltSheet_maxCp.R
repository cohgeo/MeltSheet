# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact. 
# In this case we model the cooling of the Morokweng melt sheet with a maximum
# Cp.prime value.
# This code is designed to run in R.
# Updated 2020.12.07 CH

## SETUP -----------------------------------------------------------------------

# Clear all from workspace/environment.
  # rm(list = ls())

# Install and load required packages.
  # Install and load lattice for plotting results.
  if(!require(lattice)){install.packages("lattice")}
  library(lattice)
  # Install and load grDevices for plot color schemes.
  if(!require(grDevices)){install.packages("grDevices")}
  library(grDevices)


## SET MODEL DIMENSIONS --------------------------------------------------------

# Define grid size and node spacing for the model.

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


## SET MODEL PARAMETERS --------------------------------------------------------

# Define material properties and initialize arrays for the model area and for
# the melt sheet. 
# All temperatures in °C should be converted to K before running the model.

# Set and calculate constants.
# Set latent heat of fusion of diopside.
  # Suggested value: 421000 J/kg; Abramov and Kring (2007)
  L <- 421000  # [J/kg] 
# Set liquidus temperature.
  # Suggested value for regular model: 1195°C from MELTS modeling.
  # Value to maximize Cp.prime: 1108°C (felsic, 2 wt% H20)
  T.liquidus.C <- 1108  # [°C]
  T.liquidus <- T.liquidus.C + 273.15  # [K]
# Set solidus temperature.
  # Suggested value for regular model: 796°C from MELTS modeling.
  # Value to maximize Cp.prime: 825°C (felsic, 0.1 wt% H20)
  T.solidus.C <- 825  # [°C]
  T.solidus <- T.solidus.C + 273.15  # [K]
# Set Zr saturation temperature from MELTS modeling (not adjusted for changing 
# solidus and liquidus to maximize Cp.prime)
  T.Zrsat.C <- 880  # [°C]
# Set heat capacity with no latent heat.
  # Suggested value: 1000 J/(kg*K), heat capacity of basement, melt, breccia; 
  # Abramov and Kring (2007)
  Cp.nolatent <- 1000  # [J/(kg*K)]
# Set initial heat capacity.
  Cp.initial <- matrix(data = Cp.nolatent, nrow = x.num, 
                       ncol = z.num)  # [J/(kg*K)]
# Calculate heat capacity accounting for latent heat of fusion.
  # Suggested value: 3339 J/(kg*K); Abramov and Kring (2007)
  Cp.prime <- Cp.initial + (L / (T.liquidus - T.solidus))  # [J/(kg*K)]
# Initialize matrix for initial model setup below.
  Cp.0 <- Cp.initial  
# Set thermal conductivity.
  # Suggested value: 2.5 W/(m*K), thermal conductivity of basement, melt, 
  # breccia; Abramov and Kring (2007)
  k <- matrix(data = 2.5, nrow = x.num, ncol = z.num)  #[W/(m*K)]
# Set heat productiion.
  # Suggested value: 2.2e-6 is average A value from Jones (1988), table 6
  A <- matrix(data = 2.2e-6, nrow = x.num, ncol = z.num)  # [W/(m^3)]
# Set density.
  # Suggested value: 2700 kg/(m^3), density of basement, melt, breccia; 
  # Abramov and Kring (2007)
  rho <- matrix(data = 2700, nrow = x.num, ncol = z.num)  # [kg/(m^3)]
# Calculate thermal diffusivity.
  Kappa <- k / (rho * Cp.nolatent)   # [(m^2)/s]
  Kappa.min <- min(Kappa)  # [(m^2)/s] 
# Set tectonic uplift rate.
  # We assumed no tectonic uplift. 
  U <- matrix(data = 0, nrow = x.num, ncol = z.num)  # [m/s]
# Set geothermal gradient.
  # Suggested value: 12 K/km, low end of gradient from Jones (1992) 
  a.Kkm <- 12  # [K/km]
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

# Set material properties of the melt sheet.
# Set the temperature of melt sheet.
  # Suggested value: 2370°C; Timms et al. (2017)  
  # Another option: 1700°C; Abramov and Kring (2007)
  T.meltsheet.C <- 2370  # [°C]
  T.meltsheet <- T.meltsheet.C + 273.15  # [K]
# Set dimensions of melt sheet.
  meltsheetdiameter <- 30000  # [m]
  meltsheetvolume <- 800  # [km^3]
# Calculate melt sheet thickness.
  meltsheetvolume.m3 <- meltsheetvolume * (1000 ^ 3)
  meltsheetthickness <- meltsheetvolume.m3 / 
    (((meltsheetdiameter / 2) ^ 2) * pi)  # [m]
# Find the position of the melt sheet in the modeling domain,
  # If the melt sheet is centered in the x direction in the model, find the 
  # position in meters along the x axis of the left (start) and right (end) 
  # sides of the meltsheet.
  xstart.melt <- mean(x) - (meltsheetdiameter / 2)  # [m]
  xend.melt <- mean(x) + (meltsheetdiameter / 2)  # [m]
  # Find the index in the vector that corresponds to the positions above.
  xstart.melt.index <- which(abs(x - xstart.melt) == min(abs(x - xstart.melt)))
  xend.melt.index <- which(abs(x - xend.melt) == min(abs(x - xend.melt)))
  # Find the z positions of the melt sheet. Make the melt sheet start on the 
  # second row of the matrix, below the air layer.
  zstart.melt.index <- 2 
  # If the melt sheet is meltsheetthickness and starts at the second row of the
  # model, find the index that corresponds to the bottom of the melt sheet.
  zend.melt <- z[2] + meltsheetthickness
  zend.melt.index <- round(zend.melt / dz, digits = 0)
  # Rename zend.melt for clarity later on.
  base.MS <- zend.melt / 1000  # [km]
# Input the sample height above base of melt sheet.
  M3.1 <- 0.7532
  M3.2 <- 0.4696
  M3.3 <- 0.2596
  M3.4 <- 0.1701
  M3.6 <- 0.1050
  
# Set a value for the air (surface) temperature.
  T.surface.C <- 25  # [°C]
  T.surface <- T.surface.C + 273.15  # [K]
  
# # Save table of parameters as a .csv file. Uncomment if needed.
# # Create and populate a data frame of model parameters (MP).
#   MP <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3))
#   colnames(MP) <- c("Variable", "Units", "Value")
#   # Model domain and resolution.
#   MP[1, ]  <- c("Model domain and resolution", "", "")
#   MP[2, ]  <- c("Modeled domain, distance", "km", x.size / 1000)
#   MP[3, ]  <- c("Model step size, distance", "m", round(dx, digits = 1))
#   MP[4, ]  <- c("Modeled domain, depth", "km", z.size / 1000)
#   MP[5, ]  <- c("Model step size, depth", "m", round(dz, digits = 1))
#   MP[6, ]  <- c("Model time step", "years", round(dt.choice.yr, digits = 1))
#   # Melt sheet size.
#   MP[7, ]  <- c("Melt sheet dimensions", "", "")
#   MP[8, ]  <- c("Melt sheet diameter", "km", meltsheetdiameter / 1000)
#   MP[9, ]  <- c("Melt sheet volume", "km^3", meltsheetvolume)
#   MP[10, ] <- c("Melt sheet thickness", "km", round(meltsheetthickness / 1000, digits = 2))
#   # Melt sheet and basement rock properties.
#   MP[11, ] <- c("Physical properties of basement and melt sheet", "", "")
#   MP[12, ] <- c("Geothermal gradient", "°C/km", a.Kkm)
#   MP[13, ] <- c("Latent heat of fusion of diopside", "J/kg", L)
#   MP[14, ] <- c("Heat capacity, no latent heat effects", "J/(kg*K)", Cp.nolatent)
#   MP[15, ] <- c("Heat capacity, accounting for latent heat of fusion", "J/(kg*K)", round(Cp.prime[1, 1], digits = 1))
#   MP[16, ] <- c("Thermal conductivity", "W/(m*K)", k[1, 1])
#   MP[17, ] <- c("Heat production, basement and melt sheet", "W/(m^3)", A[1, 1])
#   MP[18, ] <- c("Density, basement and melt sheet", "kg/(m^3)", rho[1,1])
#   MP[19, ] <- c("Thermal diffusivity, basement and melt", "(m^2)/s", round(Kappa[1, 1], digits = 7))
#   MP[20, ] <- c("Tectonic uplift", "m/s", U[1, 1])
#   # Temperatures.
#   MP[21, ] <- c("Temperatures", "", "")
#   MP[22, ] <- c("Temperature, melt sheet at time = 0", "°C", T.meltsheet.C)
#   MP[23, ] <- c("Temperature, surface", "°C", T.surface.C)
#   MP[24, ] <- c("Temperature, liquidus", "°C", T.liquidus.C)
#   MP[25, ] <- c("Temperature, solidus", "°C", T.solidus.C)
#   MP[26, ] <- c("Temperature, Zr saturation", "°C", T.Zrsat.C)
# # Write table to csv. Change path for your working directory.
#   write.csv(MP, file = "Morokweng Results/ModelParameters_2020.05.28.csv")
  

## INITIAL T DISTRIBUTION ------------------------------------------------------

# Set up an initial condition of an instantaneously emplaced melt sheet at an 
# elevated temperature and a central uplift region.

# Initialize a matrix for initial temperature conditions.
  T.0 <- matrix(data = 0, nrow = x.num, ncol = z.num)  # [K]
# Assign the first row of the initial matrix to be the air temperature.
  T.0[ , 1] <- T.surface  # [K] 
# Make a geotherm for the initialization matrix.
  for (i in 2:z.num) {
    T.0[ , i] <- T.0[1, 1] + a * z[i]
  }
  
# Elevate the temperature at central uplift.
  # Elevate the temperature at the central uplift by creating an ellipse with a
  # temperature gradient that is centered above the model domain (in the sky).
  # Here we used a large circle with a diameter two times larger than the melt 
  # sheet diameter, centered above the melt sheet. We scaled the temperature by 
  # position from center of ellipse and replaced these temperatures on relevant
  # nodes of the T.0 initialization matrix.
  for (i in 1:x.num) {
    for (j in 1:z.num) {
      # Calculate the value at the node to determine if it inside or outside the 
      # ellipse.
      node.value <- (((x[[i]] - mean(x))^2) / ((meltsheetdiameter)^2)) + 
        (((z[[j]] - (-meltsheetdiameter / 2))^2) / ((meltsheetdiameter) ^ 2))
      # Assign a value to the uplift matrix that depends on whether the node is
      # inside out outside the uplift ellipse.
      if (node.value <= 1) {
        # Inside the ellipse, scale the temperature position to decrease from 
        # center of ellipse. 
        T.0[i, j] <- 1600 + -node.value * 1150
      } else {
        # Outside the ellipse, keep the geotherm the same temperature.
        T.0[i, j] <- T.0[i, j] 
      }
    }
  }
# Check to see how the uplift geometry and scaling looks (uncomment next line).
  # levelplot(T.0 - 273.15, col.regions = rev(hcl.colors(n = 25, palette = 'Spectral')))

# Emplace the melt sheet on the initial conditions matrix. 
  # Set material properties of the melt sheet.
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
# Reapply air layer to top row of the model to cover the overwriting by the code
# that inserts the central uplift region.
  T.0[, 1] <- T.surface  # [K]

# Plot initial conditions to verify that the model is set up correctly.
levelplot(T.0 - 273.15,
  # Change x and y axis scales and labeling.
  scales = list(x = list(labels = seq(from = 0, to = (x.size / 1000), by = 10)),
              y = list(labels = seq(from = 0, to = (z.size / 1000), by = 5))),
  # Flip y axis so depth increases down.
  ylim = c(200, 0),
  # Add plot labels and title.
  xlab = "distance (km)", 
  ylab = "depth (km)", 
  main = "t = 0",
  # Set the aspect ratio of plot.
  aspect = (z.size / 2) / x.size,
  # Set the color scheme.
  col.regions = hcl.colors(n = 12.5, palette = "Blue-Red"),
  # Turn on contour lines.
  contour = T,
  # Set colorbar limits and density of contour interval on plot.
  # Uncomment the next line to color all temperatures.
  # at = c(seq(0, 2400, length.out = 12.5)))
  # Or, using the next line: color a restricted range of temperatures to better 
  # see the gradient of the central uplift region.
  at = c(seq(0, 1000, length.out = 11)))


## CREATE FUNCTION OF THERMAL MODEL --------------------------------------------

# Make a function to solve the finite difference equation for n time steps.
  RunModel <- function(n.timesteps) {
    # Initialize a list as a container for model returns at different times.
    T.n <- list()
    # Place T.0 in the first position of the T.n list.
    T.n[[1]] <- T.0
    # Initialize a list as a container for model returns at different times.
    Cp <- list()
    # Place T.0 in the first position of the T.n list.
    Cp[[1]] <- Cp.0
    for (t in 1:n.timesteps) {
      # Initialize T.n[[t + 1]]
      T.n[[t + 1]] <- matrix(data = NA, nrow = x.num, ncol = z.num)
      # Initialize Cp[[t + 1]]
      Cp[[t + 1]] <- matrix(data = NA, nrow = x.num, ncol = z.num)
      # Calculate T for interior of T.n[[t + 1]] matrix.
      for (i in 2:(x.num - 1)) {
        for (j in 2:(z.num - 1)) {
          # Set constant T border.
          T.n[[t + 1]][, 1] <- T.n[[1]][1, 1]  # top
          T.n[[t + 1]][, z.num] <- T.n[[1]][x.num, z.num]  # bottom
          T.n[[t + 1]][x.num, ] <- T.n[[t]][x.num, ]  # right
          T.n[[t + 1]][1, ] <- T.n[[1]][1, ]  # left
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
      # next iteration of the model.
      # Set constant Cp border.
      Cp[[t + 1]][1, ] <- Cp.initial[1, ]  # left
      Cp[[t + 1]][x.num, ] <- Cp.initial[x.num, ]  # right
      Cp[[t + 1]][, 1] <- Cp.initial[, 1]  # top
      Cp[[t + 1]][, z.num] <- Cp.initial[, z.num]  # bottom
      # Calculate Cp for interior of the model.
      for (i in 2:(x.num - 1)) {
        for (j in 2:(z.num - 1)) {
          Cp[[t + 1]][i, j] <- ifelse(T.n[[t + 1]][i, j] < T.liquidus && T.n[[t + 1]][i, j] > T.solidus, Cp.prime[i, j], Cp.initial[i, j])
        }
      }
    }
    objects.to.return <- list("T.n" = T.n, "Cp" = Cp)
    return(objects.to.return)
  }


## RUN MODEL -------------------------------------------------------------------

# Solve the finite difference equation for n time steps.
  # Set number of time steps. 
  n <- 800 
  # Run model for n timesteps.
  model.results <- RunModel(n)  
  # Save the model results of temperature to R environment.
  T.n <- model.results$T.n  # [K]

# Save model results to computer. Change path for your working directory.
  saveRDS(T.n, file = "Morokweng Results/MeltSheet Results_2020.12.07/T.n_Cp.prime.max.RData")
# Load in previous model run results. Change path for your working directory.
  # T.n <- readRDS("Morokweng Results/MeltSheet Results_2020.06.23/T.n_1700.RData")
  
  
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
    ylim = c(100, 0),
    # Add plot labels and title.
    xlab = "distance (km)", # x axis label
    ylab = "depth (km)", # y axis label
    main = paste("t =", (round(((IterationNumber - 1) * dt.choice.yr),
                 digits = 0)) / 1000, "ka", sep = " "),
    # Set the aspect ratio of plot.
    aspect = (zval / 2) / xval,
    # Set the color scheme.
    col.regions = hcl.colors(n = 12.5, palette = "Blue-Red"),
    # Comment out the line above and uncomment the lind below to plot in 
    # grayscale.
    # col.regions = rev(hcl.colors(n = 25, palette = 'Grays')),  
    contour = T,  # Add contour lines.
    # Set colorbar limits and density of contour interval on plot.
    at = c(seq(0, 2400, length.out = 12.5))) 
  return(OutputPlot)
  }

# Make plotting function to plot model results in °C zoomed in on the top four 
# kilometers of the model. See function above for code comments.
  levelplotCHzoom <- function(IterationNumber,
                              x.axis.in.m = x,
                              z.axis.in.m = z,
                              dtyr = dt.yr,
                              xnum = x.num,
                              znum = z.num,
                              xval = x.size,
                              zval = z.size) {
    InputMatrix <- T.n[[IterationNumber]]          
    OutputPlot <- levelplot(InputMatrix - 273.15,
      scales = list(x = list(labels = seq(from = 0, to = (xval / 1000), 
                                          by = 10)),
                    y = list(labels = seq(from = 0, to = (zval / 1000), 
                                          by = 1))),
      ylim = c(40, 0),
      xlab = "distance (km)",
      ylab = "depth (km)",
      main = paste("t =", (round(((IterationNumber - 1) * dt.choice.yr), 
                   digits = 0)) / 1000, "ka", sep = " "),
      aspect = (zval / 2) / xval,
      col.regions = hcl.colors(n = 12.5, palette = "Blue-Red"),
      contour = T,
      at = c(seq(0, 2400, length.out = 12.5)))
    return(OutputPlot)
  }

# Make plotting function to plot model results in °C zoomed in on the top four 
# kilometers of the model with only contour lines, no color. 
# See function above for code comments.
  levelplotCHzoom.nocolor <- function(IterationNumber,
                                      x.axis.in.m = x,
                                      z.axis.in.m = z,
                                      dtyr = dt.yr,
                                      xnum = x.num,
                                      znum = z.num,
                                      xval = x.size,
                                      zval = z.size) {
    InputMatrix <- T.n[[IterationNumber]]                               
    OutputPlot <- levelplot(InputMatrix - 273.15,
      scales = list(x = list(labels = seq(from = 0, to = (xval / 1000), 
                                          by = 10)),
                    y = list(labels = seq(from = 0, to = (zval / 1000), 
                                         by = 1))),
      ylim = c(40, 0),
      xlab = "distance (km)", 
      ylab = "depth (km)",
      main = paste("t =", (round(((IterationNumber - 1) * dt.choice.yr), 
                  digits = 0)) / 1000, "ka", sep = " "),
      aspect = (zval / 2) / xval,
      region = FALSE,
      # col.regions = hcl.colors(n = 12.5, palette = "Blue-Red"),
      contour = T, 
      at = c(seq(0, 2400, length.out = 12.5))) 
    return(OutputPlot)
  }

# Check to see that the model ran correctly and that the plotting functions
# work (uncomment the next lines).
  # levelplotCH(1)
  # levelplotCH(51)
  # levelplotCH(801)
  # levelplotCHzoom(1)
  # levelplotCHzoom(51)
  # levelplotCHzoom(801)
  # levelplotCHzoom.nocolor(1)
  # levelplotCHzoom.nocolor(51)
  # levelplotCHzoom.nocolor(801)

# Output plots as a PDF.
  # pdf("Morokweng Results/MeltSheet Results_2020.06.03/MeltSheet_Morokweng_2D results_2020.06.03.pdf",
  # width = 8.5, height = 11)
  # keyplots <- c(1, 51, 101, 151, 201, 251, 301, 401, 501)
  # levelplotCH(1)
  # for (n in 1:length(keyplots)) {
  #   print(levelplotCHzoom(keyplots[[n]]))
  # }
  # for (n in 1:length(keyplots)) {
  #   print(levelplotCHzoom.nocolor(keyplots[[n]]))
  # }
  # dev.off()


## EXTRACT AND PLOT 1D PROFILE RESULTS FOR M3 ----------------------------------

# Extract 1D results at the location of the M3 borehole.

# Iterate through T.n list, extract temperature values through 70.3% away from 
# the center of  the melt sheet (location of borehole M3), save in list T.1D.  
  # Initialize list to hold center values from T.n matricies.
  T.1D <- list()
  # Find the node location of M3.
  M3 <- round((mean(x) - ((meltsheetdiameter / 2) * 0.703)) / dx, digits = 0)
  # Iterate through all matrices in list T.n to extract temperature values at M3.
  for (i in 1:length(T.n)) {
    T.1D[[i]] <- T.n[[i]][(M3), ]
  }
# Extract key results to plot.
  # Initialize list to hold key results.
  T.1D.lim <- list()
  # Set key results.
  keyresults <- c(1, 51, 101, 151, 201, 251, 301, 401, 501, 601, 701, 801)
  # Iterate through keyresults and save a list of those temperature profiles.
  for (n in 1:length(keyresults)) {
    T.1D.lim[[n]] <- T.1D[[keyresults[[n]]]]
  }

# Plot 1D profile through M3 borehole. 
  plot.new()  
  plot(range(0, 2500), range(0, 4),
       type = "n",
       xlab = "temperature (°C)",
       ylab = "depth (km)",
       ylim = c(4, 0),
       yaxs = "i",
       xlim = c(0, 2500),
       xaxs = "i")
  # Set number of colors based on number of model results being plotted.
  colors <- (hcl.colors(n = length(T.1D.lim), palette = "Spectral"))
  # Add lines.
  for (i in 1:length(T.1D.lim)) {
    tempprofile <- T.1D.lim[[i]] - 273.15
    lines(tempprofile, z / 1000, lwd = 1.75,  col = colors[i])
  }
  # Add title.
  title("1D profile of borehole M3", cex.main = 0.8)
  # Add vertical lines for liquidus and solidus.
  lines(c(T.liquidus.C, T.liquidus.C), c(0, 4), lwd = 2, col = "gray")
  lines(c(T.solidus.C, T.solidus.C), c(0, 4), lwd = 2, col = "gray")
  text(x = T.liquidus.C + 30, y = 3.5, labels = "liquidus", srt = 270)
  text(x = T.solidus.C + 30, y = 3.5, labels = "solidus", srt = 270)
  # Add horizontal lines for sample horizons.
  lines(c(0, 2500), c(base.MS - M3.1, base.MS - M3.1))  # M3-1
  lines(c(0, 2500), c(base.MS - M3.2, base.MS - M3.2))  # M3-2
  lines(c(0, 2500), c(base.MS - M3.3, base.MS - M3.3))  # M3-3
  lines(c(0, 2500), c(base.MS - M3.4, base.MS - M3.4))  # M3-4
  lines(c(0, 2500), c(base.MS - M3.6, base.MS - M3.6))  # M3-6
  lines(c(0, 2500), c(base.MS, base.MS), col = "blue")  # base of melt sheet
  # Add a legend. 
  legend("bottomright", 
         paste(round(((keyresults - 1) * dt.choice.yr), digits = 0) / 1000), 
         cex = 0.8, col = colors, lty = "solid", lwd = 1.75,
         title = "time elapsed (ka)")

# Zoom in on samples to see when each sample passes through the solidus.
  # Extract limited results to plot.
  T.1D.lim2 <- list()
  lim2.iterations <- c(105, 106, 107,  # M3-1
                       171, 172, 173,  # M3-2
                       217, 218, 219,  # M3-3
                       235, 236, 237,  # M3-4
                       249, 250, 251)  # M3-6
  for (n in 1:length(lim2.iterations)){
    T.1D.lim2[[n]] <- T.1D[[lim2.iterations[[n]]]]
  }
  # Make plot.
  plot.new()  
  plot(range(800, 850), range(0.4, 1.2),
       type = "n",
       xlab = "Temperature (°C)",
       ylab = "Depth (km)",
       ylim = c(1.2, 0.4),
       yaxs = "i",
       xlim = c(800, 850),
       xaxs = "i")
  colors <- (hcl.colors(n = length(T.1D.lim2), palette = "Spectral"))
  for (i in 1:length(T.1D.lim2)) {
    tempprofile <- T.1D.lim2[[i]] - 273.15
    lines(tempprofile, z / 1000, type = "l", lwd = 1.75, col = colors[i],
          lty = "solid")
  }
  title("1D profile of borehole M3", cex.main = 0.8)
  lines(c(T.solidus.C, T.solidus.C), c(0, 1.5), lwd = 2, col = "gray")
  text(x = T.solidus.C + 1, y = 0.6, labels = "solidus", srt = 270)  
  lines(c(0, 2500), c(base.MS - M3.1, base.MS - M3.1))  # M3-1
  lines(c(0, 2500), c(base.MS - M3.2, base.MS - M3.2))  # M3-2
  lines(c(0, 2500), c(base.MS - M3.3, base.MS - M3.3))  # M3-3
  lines(c(0, 2500), c(base.MS - M3.4, base.MS - M3.4))  # M3-4
  lines(c(0, 2500), c(base.MS - M3.6, base.MS - M3.6))  # M3-6
  legend("bottomright", paste((lim2.iterations - 1) / 10), cex = 0.8, 
         col = colors, lty = "solid", lwd = 1.75, title = "time elapsed (ka)")
 

## EXTRACT AND PLOT 1D PROFILE RESULTS FOR CENTER OF MELT SHEET ----------------
  
# Extract 1D results through center of 2D model results.
 
# Iterate through T.n list, extract temperature values through center of melt 
# sheet, save in list T.1D.center.
  # Initialize list to hold center values from T.n matricies.
  T.1D.center <- list()
  # Iterate through all matrices in list T.n to extract temperature values at 
  # the center of the melt sheet.
  for (i in 1:length(T.n)) {
    T.1D.center[[i]] <- T.n[[i]][(x.num / 2), ]
  }
# Extract key results to plot.
  # Initialize list to hold key results.
  T.1D.center.lim <- list()
  # Set key results.
  keyresults.center <- c(1, 51, 101, 151, 201, 251, 301, 401, 501, 601, 701, 801)
  # Iterate through keyresults and save a list of those temperature profiles.
  for (n in 1:length(keyresults.center)) {
    T.1D.center.lim[[n]] <- T.1D.center[[keyresults.center[[n]]]]
  }
  
# Plot 1D profile through the center of the melt sheet. 
  plot.new()  
  plot(range(0, 2500), range(0, 4),
       type = "n",
       xlab = "temperature (°C)",
       ylab = "depth (km)",
       ylim = c(4, 0),
       yaxs = "i",
       xlim = c(0, 2500),
       xaxs = "i")
  colors <- (hcl.colors(n = length(T.1D.center.lim), palette = "Spectral"))
  for (i in 1:length(T.1D.center.lim)) {
    tempprofile <- T.1D.center.lim[[i]] - 273.15
    lines(tempprofile, z / 1000, lwd = 1.75,  col = colors[i])
  }
  title("1D profile of center of melt sheet", cex.main = 0.8)
  lines(c(T.liquidus.C, T.liquidus.C), c(0, 4), lwd = 2, col = "gray")
  lines(c(T.solidus.C, T.solidus.C), c(0, 4), lwd = 2, col = "gray")
  text(x = T.liquidus.C + 30, y = 3.5, labels = "liquidus", srt = 270)
  text(x = T.solidus.C + 30, y = 3.5, labels = "solidus", srt = 270)
  lines(c(0, 2500), c(base.MS - 0.7566, base.MS - 0.7566))  # M3-1
  lines(c(0, 2500), c(base.MS - 0.4711, base.MS - 0.4711))  # M3-2
  lines(c(0, 2500), c(base.MS - 0.2599, base.MS - 0.2599))  # M3-3
  lines(c(0, 2500), c(base.MS - 0.1716, base.MS - 0.1716))  # M3-4
  lines(c(0, 2500), c(base.MS - 0.1081, base.MS - 0.1081))  # M3-6
  lines(c(0, 2500), c(base.MS, base.MS), col = "blue")  # base of melt sheet
  lines(c(0, 2500), 
        c(2 * meltsheetthickness / 1000, 2 * meltsheetthickness / 1000), 
        col = "red")  # cutoff for melt availability for the melt sheet system
  legend("bottomright", 
         paste(round(((keyresults.center - 1) * dt.choice.yr), digits = 0) / 1000), 
         cex = 0.8, col = colors, lty = "solid", lwd = 1.75,
         title = "time elapsed (ka)")
  
# Find when 2 * meltsheetthickness passed through the solidus at the center of
# the model.
  T.1D.lim.meltavail <- list()
  lim.meltavail.iterations <- c(695, 696, 697, 698, 699)
  # lim.meltavail.iterations <- c(541, 542, 543, 544, 545, 546)
  for (n in 1:length(lim.meltavail.iterations)){
    T.1D.lim.meltavail[[n]] <- T.1D.center[[lim.meltavail.iterations[[n]]]]
  }
  # Make plot.
  plot.new()  
  plot(range(750, 800), range(2.30, 2.24),
       type = "n",
       xlab = "Temperature (°C)",
       ylab = "Depth (km)",
       ylim = c(2.30, 2.24),
       yaxs = "i",
       xlim = c(750, 800),
       xaxs = "i")
  colors <- (hcl.colors(n = length(T.1D.lim.meltavail), palette = "Spectral"))
  for (i in 1:length(T.1D.lim.meltavail)) {
    tempprofile <- T.1D.lim.meltavail[[i]] - 273.15
    lines(tempprofile, z / 1000, type = "l", lwd = 1.75, col = colors[i],
          lty = "solid")
  }
  title("1D profile of center of melt sheet", cex.main = 0.8)
  lines(c(T.solidus.C, T.solidus.C), c(2.4, 2.1), lwd = 2, col = "gray")
  text(x = T.solidus.C - 1, y = 2.295, labels = "solidus", srt = 270)  
  lines(c(0, 2500), 
        c(2 * meltsheetthickness / 1000, 2 * meltsheetthickness / 1000), 
        col = "red")  # cutoff for melt availability for the melt sheet system
  legend("bottomright", paste((lim.meltavail.iterations - 1) / 10), cex = 0.8, 
         col = colors, lty = "solid", lwd = 1.75, title = "time elapsed (ka)")  
  
  
## PLOT MODEL RESULTS THROUGH TIME FOR EACH SAMPLE -----------------------------
  
# Extract a temperature value for each time step (each matrix in the T.n list)
# at a given depth and distance from the center of the melt sheet (sample 
# location).
  # Initialize matrix to hold temperature values from T.n matricies for all 
  # samples.
  T.M3 <- matrix(data = NA,
                 nrow = 5, # a row for each sample  
                 ncol = length(T.n))  # a column for each time step
  # Rename T.M3 rows by sample.
  rownames(T.M3) <- c("M3-1", "M3-2", "M3-3", "M3-4", "M3-6")
# Determine which row of the model (z axis) results matrices corresponds to each
# sample depth.
  # Initialize vector to hold sample depths.
  z.samples <- vector()
  # Find nodal depth of M3-1.
  z.samples[[1]] <- round(((base.MS - M3.1) * 1000) / dz, digits = 0)  # M3-1
  z.samples[[2]] <- round(((base.MS - M3.2) * 1000) / dz, digits = 0)  # M3-2
  z.samples[[3]] <- round(((base.MS - M3.3) * 1000) / dz, digits = 0)  # M3-3
  z.samples[[4]] <- round(((base.MS - M3.4) * 1000) / dz, digits = 0)  # M3-4
  z.samples[[5]] <- round(((base.MS - M3.6) * 1000) / dz, digits = 0)  # M3-6
# Iterate through T.n list, extract T values at the M3 borehole for each sample,
# save in matrix T.M3.
  for (n in 1:length(z.samples)) {
    for (i in 1:length(T.n)) {
      T.M3[n, i] <- T.n[[i]][M3, z.samples[[n]]]
    }
  }
# Make a vector of values for x axis.
  xax.M3 <- vector()
  xax.M3 <- 0:(length(T.n) - 1)
  xax.M3 <- xax.M3 * dt.choice.yr / 1000 # ka
# Add a row to T.M3 of time in ka.
  T.M3 <- rbind(T.M3, xax.M3)
  rownames(T.M3)[6] <- "elapsed time (ka)"

# Plot T-t path for each sample with a zoomed in inset plot.
  plot.new()
  par(fig = c(0, 1, 0, 1))
  plot(range(0, 80), range(0,2500),
       type = "n",
       xlab = "elapsed time (ka)",
       ylab = "temperature (°C)",
       xaxs = "i",
       yaxs = "i")
  colors <- (hcl.colors(n = (dim(T.M3)[1] - 1), palette = "Temps"))
  for (i in 1:((dim(T.M3)[1]) - 1)) {
    lines(x = T.M3[6, ], y = T.M3[i, ] - 273.15,
          type = "l", lwd = 1.75, col = colors[i], lty = "solid")
  }
  title("Time-temperature path of M3 borehole samples", cex.main = 0.8)
  lines(c(0, 80), c(T.liquidus.C, T.liquidus.C), lwd = 2, col = "gray")
  lines(c(0, 80), c(T.solidus.C, T.solidus.C), lwd = 2, col = "gray")
  text(x = 60, y = T.liquidus.C + 30, labels = "liquidus")
  text(x = 60, y = T.solidus.C + 30, labels = "solidus")
  legend("bottomleft", rownames(T.M3)[1:((dim(T.M3)[1]) - 1)], 
         cex = 0.8, col = colors, lty = "solid", lwd = 1.75, title = "sample")
  # Make inset figure.
  par(fig = c(0.25, 1, 0.5, 1), new = TRUE)
  plot(range(5, 30), range(700,900),
       type = "n",
       xlab = "elapsed time (ka)",
       ylab = "temperature (°C)",
       xaxs = "i",
       yaxs = "i")
  colors <- (hcl.colors(n = (dim(T.M3)[1] - 1), palette = "Temps"))
  for (i in 1:((dim(T.M3)[1]) - 1)) {
    lines(x = T.M3[6, ], y = T.M3[i, ] - 273.15,
          type = "l", lwd = 1.75,col = colors[i], lty = "solid")
  }
  lines(c(5, 30), c(T.solidus.C, T.solidus.C), lwd = 2, col = "gray")
  text(x = 27, y = T.solidus.C + 5, labels = "solidus")

  
## PLOT U-Pb VS. THERMAL MODEL AGES WITH DEPTH ---------------------------------
  
# Use relative.ages data frame to plot time-depth plots for U-Pb weighted mean 
# ages and thermal model ages
  
# Find thermal model anchoring point (mean variance of upper four samples).
  anchor.model <- mean(relative.ages$solidus.ka[1:4])
# Find U-Pb anchoring point of upper four samples.
  # anchor.UPb.simple <- mean(relative.ages$UPb.age[1:4])
  # Find max and min ages of y axis.
  max.age <- max(relative.ages$UPb.age[1:4]) + 
    5 * median(relative.ages$UPb.uncertainty[1:4])
  min.age <- min(relative.ages$UPb.age[1:4]) - 
    5 * median(relative.ages$UPb.uncertainty[1:4])
  # Use compound.prob function to create a summed pdf from min.age to max.age.
  # INPUTS:      ages = a vector of ages
  #              sigs = a vector of 2 sigma uncertainties
  # OUTPUTS:     P    = probability values
  compound.prob <- function(ages, sigs){
    # Divide 2 sigma uncertainties to get 1 sigma uncertainties.
    sigs1 <- sigs / 2  
    x <- seq(min(min.age), max(max.age), length = ((max.age - min.age) * 100))  
    interval <- matrix(0, nrow = length(x), ncol = length(ages))
    for (i in 1:length(ages)) {
      interval[, i] <- dnorm(x, ages[i], 
                             sigs1[i]) / length(ages) * mean(diff(x))
    }
    P <- data.frame(x = x,probability = apply(interval, 1, sum))
    return(P)
  }
  # Calculate summed PDFs of weighted mean ages the compound.prob function.
  P <- compound.prob(relative.ages$UPb.age[1:4], 
                     relative.ages$UPb.uncertainty[1:4]) 
  # Scale the column by the max value in the column so the y axis will have a 
  # max value of 1.
  P$probability.scaled <- P$probability / max(P$probability) 
# Make plot of summed pdf to visualize.
  plot.new()  
  par(fig = c(0, 1, 0, 1))
  plot(range(145.98, 146.1), range(0, 2),
       type = "n",
       xlab = "age (Ma)",
       ylab = "relative probability",
       ylim = c(0, 1.1),
       yaxs = "i",
       xlim = c(146.1, 145.98),
       xaxs = "i")
  lines(x = P$x,
        y = P$probability.scaled,
        lwd = 2,
        col = "blue")
  title("Summed PDF of weighted mean U-Pb ages",
        cex.main = 0.8)
# Find mean of variance of weighted mean U-Pb ages and uncertainties.
  anchor.UPb <- P[(match(max(P$probability.scaled), P$probability.scaled)), 1]
  
# Make plot of U-Pb ages by stratigraphic height.
  plot.new()
  plot(range(145.98, 146.1), range(0.2, 1.4),
       type = "n",
       xlab = "time (Ma)",
       ylab = "stratigraphic position of samples",
       xaxs = "i",
       yaxs = "i",
       ylim = c(1.4, 0.2),
       xlim = c(146.1, 145.98))
  # Set number of colors based on number of model results being plotted.
  colors <- (hcl.colors(n = 5, palette = "Temps"))
  # Add error bars for U-Pb ages.
  arrows(relative.ages$UPb.age - 
           relative.ages$UPb.uncertainty, relative.ages$strat.pos,
         relative.ages$UPb.age + 
           relative.ages$UPb.uncertainty, relative.ages$strat.pos,
         length = 0.05, angle = 90, code = 3)
  # Add line with points for U-Pb time.
  lines(x = relative.ages$UPb.age, y = relative.ages$strat.pos)
  points(x = relative.ages$UPb.age, y = relative.ages$strat.pos,
         pch = 22, lwd = 0.5, col = "black", cex = 1, bg = colors)
  # Add anchor line.
  lines(x = c(anchor.UPb, anchor.UPb), y = c(0.2, 1.4))
  
# Make plot of solidus thermal model ages by stratigraphic height.
  # Need to cover range from 146.10 to 145.98 Ma (0.12 Ma)
  plot.new()
  plot(range(0, 0.12), range(0.2, 1.4),
       type = "n",
       xlab = "time (Ma)",
       ylab = "stratigraphic position of samples",
       xaxs = "i",
       yaxs = "i",
       ylim = c(1.4, 0.2),
       xlim = c(0, 0.12))
  # Set number of colors based on number of model results being plotted.
  colors <- (hcl.colors(n = 5, palette = "Temps"))
  # Add line with points for U-Pb time.
  lines(x = relative.ages$solidus.ka / 1000, y = relative.ages$strat.pos)
  points(x = relative.ages$solidus.ka / 1000, y = relative.ages$strat.pos,
         pch = 21, lwd = 0.5, col = "black", cex = 1, bg = colors)
  # Add anchor line.
  lines(x = c(anchor.model / 1000, anchor.model / 1000), y = c(0.2, 1.4))
  
# Make plot of Zr saturation thermal model ages by stratigraphic height.
  # Need to cover range from 146.10 to 145.98 Ma (0.12 Ma)
  plot.new()
  plot(range(0, 0.12), range(0.2, 1.4),
       type = "n",
       xlab = "time (Ma)",
       ylab = "stratigraphic position of samples",
       xaxs = "i",
       yaxs = "i",
       ylim = c(1.4, 0.2),
       xlim = c(0, 0.12))
  # Set number of colors based on number of model results being plotted.
  colors <- (hcl.colors(n = 5, palette = "Temps"))
  # Add line with points for U-Pb time.
  lines(x = relative.ages$Zrsat.ka / 1000, y = relative.ages$strat.pos)
  points(x = relative.ages$Zrsat.ka / 1000, y = relative.ages$strat.pos,
         pch = 21, lwd = 0.5, col = "black", cex = 1, bg = colors)
  # Add anchor line.
  lines(x = c(anchor.model / 1000, anchor.model / 1000), y = c(0.2, 1.4))

  