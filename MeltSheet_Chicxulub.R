# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact.
# In this case we model the cooling of the Chicxulub. melt sheet.
# This code is designed to run in R.
# Updated 2020.10.30 CH


## SETUP -----------------------------------------------------------------------

# Clear all from workspace/environment.
  # rm(list = ls())  # Uncomment if needed

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
  x.size <- 120000  # [m] 
  z.size <- 20000  # [m]
# Set resolution (number of nodes).  
  x.num <- 600
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
  # Suggested value: 421000 J/kg; Birch et al. (1942)
  L <- 421000  # [J/kg] 
# Set liquidus temperature.
  # Suggested value: 1177°C; Abramov and Kring (2007)
  T.liquidus.C <- 1177  # [°C]
  T.liquidus <- T.liquidus.C + 273.15  # [K]
# Set solidus temperature.
  # Suggested value: 997°C; Abramov and Kring (2007)
  T.solidus.C <- 997  # [°C]
  T.solidus <- T.solidus.C + 273.15  # [K]
# Set Zr saturation temperature unknown, set to 0 to prevent issues in code.
  T.Zrsat.C <- 0  # [°C]
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
  # Suggested value: 13 K/km, Abramov and Kring (2007)
  a.Kkm <- 13  # [K/km]
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
  meltsheetdiameter <- 60000  # [m]
  meltsheetthickness <- 3000  # [m]
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
  
# Set a value for the air (surface) temperature.
  T.surface.C <- 25  # [°C]
  T.surface <- T.surface.C + 273.15  # [K]
  
# # Save table of parameters as a .csv file.
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
# write.csv(MP, file = "Chicxulub Results/MeltSheetChicxulub_2020.09.22/ModelParameters_2370_2020.09.22.csv")
  

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
        # T.0[i, j] <- 1600 + -node.value * 1150
        T.0[i, j] <- 1700 + -node.value * 1150 
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
  scales = list(x = list(labels = seq(from = 0, to = (x.size / 1000), by = 20)),
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
  # Next line: color all temperatures.
  # at = c(seq(0, 2400, length.out = 12.5)))
  # Next line: color a restricted range of temperatures to better see the
  # gradient of the central uplift region.
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
  n <- 10000 
  # Run model for n timesteps.
  model.results <- RunModel(n)  
  # Save the model results of temperature to R environment.
  T.n <- model.results$T.n  # [K]

# Save model results to computer. Change path for your working directory.
  # saveRDS(T.n, file = "Chicxulub results/MeltSheetChicxulub_2020.09.22/T.n_2370_0to10000.RData")
# Load in previous model run results. Change path for your working directory.
  # T.n <- readRDS("Chicxulub Results/MeltSheetChicxulub Results_2020.09.22/T.n_1700_0to10000.RData")

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
    scales = list(x = list(labels = seq(from = 0, to = (xval / 1000), by = 20)),
                  y = list(labels = seq(from = 0, to = (zval / 1000), by = 5))),
    # Flip y axis so depth increases down.
    # Cut off the bottom part of the model because of edge effects.
    ylim = c(200, 0),
    # Add plot labels and title.
    xlab = "distance (km)", # x axis label
    ylab = "depth (km)", # y axis label
    main = paste("t =", (round(((IterationNumber - 1) * dt.choice.yr),
                 digits = 0)) / 1000, "ka", sep = " "),
    # Set the aspect ratio of plot.
    aspect = (zval / 2) / xval,
    # Set the color scheme.
    col.regions = hcl.colors(n = 12.5, palette = "Blue-Red"),
    # Commount out the line above and uncomment the lind below to plot in 
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
                                          by = 20)),
                    y = list(labels = seq(from = 0, to = (zval / 1000), 
                                          by = 2))),
      ylim = c(100, 0),
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
                                          by = 20)),
                    y = list(labels = seq(from = 0, to = (zval / 1000), 
                                         by = 2))),
      ylim = c(100, 0),
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
  # levelplotCH(2)
  # levelplotCH(5)
  # levelplotCH(10)
  # levelplotCH(51)
  # levelplotCH(101)
  # levelplotCHzoom(1)
  # levelplotCHzoom(51)
  # levelplotCHzoom(101)
  # levelplotCHzoom(1001)
  # levelplotCHzoom.nocolor(1)
  # levelplotCHzoom.nocolor(51)
  # levelplotCHzoom.nocolor(101)

# Output plots as a PDF. Change path for your working directory.
  # pdf("Chicxulub results/MeltSheetChicxulub_2020.09.22/MeltSheet_Chicxulub_2D results_2020.09.22.pdf",
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


## EXTRACT 1D PROFILE RESULTS FOR CENTER OF MELT SHEET -------------------------
  
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

# Save model results to computer. Change path for your working directory.
  # saveRDS(T.1D.center, file = "Chicxulub results/MeltSheetChicxulub_2020.09.22/T.1D.center_1700_0to10000b.RData")

  
## RUN MODEL FOR SECOND ITERATION ----------------------------------------------

# Run model.
  # Replace the starting conditions with the last conditions from the previous
  # model run.
  T.0 <- model.results$T.n[[length(T.n)]]
  Cp.0 <- model.results$Cp[[length(T.n)]]
  # Remove T.n from the environment.
  rm(T.n)
  # Set the number of model iterations: 10002 to 15003
  n <- 5000 
  # Run model for n timesteps.
  model.results <- RunModel(n)  
  # Save the model results of temperature to R environment.
  T.n <- model.results$T.n  # [K]
  
# Save model results to computer. Change path for your working directory.
  # saveRDS(T.n, file = "Chicxulub results/MeltSheetChicxulub_2020.09.22/T.1D.center_1700_10002to15003.RData")

# Extract 1D results and append to previous 1D results.
  # Iterate through T.n list, extract temperature values through center of melt 
  # sheet, save in list T.1D.center.2.
  # Initialize list to hold center values from T.n matricies.
  T.1D.center.2 <- list()
  # Iterate through all matrices in list T.n to extract temperature values at 
  # the center of the melt sheet.
  for (i in 1:length(T.n)) {
    T.1D.center.2[[i]] <- T.n[[i]][(x.num / 2), ]
  }

# Combine T.1D.center and T.1D.center.2 by appending T.1D.center.2 to the bottom of the T.1D.center list.
  T.1D.center <- c(T.1D.center, T.1D.center.2)
  
# Remove repeats (model runs used as initial conditions for the subsequent 
  # model).
  # n = 10000 gives T.n with length  10001, so 10002 is a repeat
  # 10002 + 5000 = 15002, so 15003 is a repeat
  # should end up with total n + 1 (for T.0)
  T.1D.center <- T.1D.center[-c(10002, 15003)]
  
# Save model results to computer. Change path for your working directory.
  # saveRDS(T.1D.center, file = "Chicxulub results/MeltSheetChicxulub_2020.09.22/T.1D.center_1700_0to15001.RData")
  
  
## PLOT 1D PROFILE RESULTS FOR CENTER OF MELT SHEET ----------------------------
  
# Extract key results to plot.
  # Initialize list to hold key results.
  T.1D.center.lim <- list()
  # Set key results.
  # keyresults.center <- c(1, 1001, 2001, 3001, 4001, 5001, 6001, 7001, 8001, 
  #                        9001, 10001, 11001, 12001, 13001, 14001, 15001)
  keyresults.center <- c(1, 1001, 2001, 3001, 4001, 5001, 6001, 7001, 8001, 
                         9001, 10001)
  # Iterate through keyresults and save a list of those temperature profiles.
  for (n in 1:length(keyresults.center)) {
    T.1D.center.lim[[n]] <- T.1D.center[[keyresults.center[[n]]]]
  }
  
# Plot 1D profile through the center of the melt sheet. 
  plot.new()  
  plot(range(0, 2500), range(0, 10),
       type = "n",
       xlab = "temperature (°C)",
       ylab = "depth (km)",
       ylim = c(10, 0),
       yaxs = "i",
       xlim = c(0, 2500),
       xaxs = "i")
  colors <- (hcl.colors(n = length(T.1D.center.lim), palette = "Spectral"))
  for (i in 1:length(T.1D.center.lim)) {
    tempprofile <- T.1D.center.lim[[i]] - 273.15
    lines(tempprofile, z / 1000, lwd = 1.75,  col = colors[i])
  }
  title("1D profile of center of melt sheet (Chicxulub)", cex.main = 0.8)
  lines(c(T.liquidus.C, T.liquidus.C), c(0, 10), lwd = 2, col = "gray")
  lines(c(T.solidus.C, T.solidus.C), c(0, 10), lwd = 2, col = "gray")
  lines(c(800, 800), c(0, 10), lwd = 2, col = "gray")
  lines(c(0, 2500), c((2 * meltsheetthickness / 1000), 
                      (2 * meltsheetthickness / 1000)), 
        lwd = 2, col = "black", lty = "dashed")
  text(x = T.liquidus.C + 30, y = 8, labels = "liquidus = 1177°C", srt = 270)
  text(x = T.solidus.C + 30, y = 8, labels = "solidus = 997°C", srt = 270)
  text(x = 800 + 30, y = 8, labels = "800°C", srt = 270)
  legend("bottomright", 
    paste(round(((keyresults.center - 1) * dt.choice.yr), digits = 0) / 1000), 
    cex = 0.8, col = colors, lty = "solid", lwd = 1.75,
    title = "time elapsed (ka)")

  