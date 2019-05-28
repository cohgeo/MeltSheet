# This script calculates a steady state geotherm based on a finite difference 
# model with constant heat flow with equations derived in class on 2/27/19 and 
# plots the geotherm.
# Updated 2019.03.05 CH.


## SETUP -----------------------------------------------------------------------

# Clear environment.
  rm(list = ls())

# Load required packages 

  # Load ggplot2.
  # install.packages("ggplot2")
  library(ggplot2)

  # Load reshape2.
  # install.packages("reshape2")
  library(reshape2)

  # Load RColorBrewer.
  # install.packages("RColorBrewer")
  library(RColorBrewer)


## INPUT AND CALCULATE PARAMETERS ----------------------------------------------

# Input constants and values. Convert values to the correct units for 
# calculations.

  # Input thermal conductivity in J/(s*K*m) or W/(m*K).
  k <- 3  # [W/(m*K)]

  # Input density in kg/(m^3).
  rho <- 2800  # [kg/(m^3)]

  # Input heat capacity in J/(kg*K).
  C <- 1e+03  # [J/(kg*K)]

  # Input heat production in W/(m^3).
  A <- 1e-6  # [W/(m^3)]

  # Input vertical velocity / uplift in mm/yr.
  U.mmyr <- 0  # [mm/yr]
  U <- U.mmyr / (31500000 * 1000)  # [m/s]
  
  # Input geothermal gradient in K/km.
  a.Kkm <- 20  # [K/km]
  a <- a.Kkm / 1000  # [K/m]

# Calculate parameters based on constants input above.

  # Calculate thermal diffusivity in (m^2)/s
  K <- k / (C * rho)  # [(m^2)/s]
  
  
## SET BOUNDARY CONDTIONS ------------------------------------------------------

# Input boundary condtions.
  
  # Input the surface temperature.
  Ts <- 20  # [°C]
  
  # Input the basal temperature.
  Tb <- 1020  # [°C]
  

## SET MODEL STEP SIZES --------------------------------------------------------
  
# Input the distance step for the model. Calculate the stable time step for the 
# model. Convert values to the correct units for calculations.
  
  # Input the change in distance.
  delta.x <- 1000  # [m]
  # delta.x2 <- delta.x ^ 2  # [m^2]
  
  # Calculate the change in time based on the stability criterion for the model.
  delta.t.s <- (delta.x ^ 2) / (2.9574 * K)  # [s]
  delta.t.Ma <- delta.t.s / 31560000000000  # [Ma]
  # K.delta.t <- K * delta.t.s   # [m^2]
  
  # Create a vector of distances over which to evaluate the model.
  max.depth.km <- 50  # Set maximum depth of model in km.
  x <- vector()  # Initialize vector.
  x <- seq(from = 0, to = (max.depth.km * 1000), by = delta.x)
  
  # Create a vector of time over which to evaluate the model.
  max.time.Ma <- 21  # Set maximum time of model in Ma.
  t.Ma <- vector()  # Initialize vector.
  t.Ma <- seq(from = 0, to = max.time.Ma, by = delta.t.Ma)
  t.s <- vector()  # Initialize vector.
  t.s <- seq(from = 0, to = (max.time.Ma * 31560000000000), by = delta.t.s)
  

## SETUP MATRIX TO HOLD RESULTS OF MODEL ---------------------------------------
  
  # Initialize matrix to hold results of the model.
  TC <- matrix(data = NA, nrow = length(t.s), ncol = length(x))  
  
  # Populate first column of the matrix.
  TC[, 1] <- Ts
  
  
## PERTURB THE STARTING CONDITIONS ---------------------------------------------
  
# Create the starting conditions of the model.
  
  # Populate first row of the matrix (unperturbed condtions).
  for (i in 2:dim(TC)[2]) {
    TC[1, i] <- TC[1, 1] + a * x[i]
  }
  
  # Overwrite part of first row of the matrix for the perturbation. 
  # Uncomment and motify relevant lines below.
  # TC[1, 5] <- 500
  # TC[1, 5:7] <- 500
  TC[1, 5:7] <- c(500, 600, 500)
  
  
## SOLVE FINITE DIFERENCE MODEL ------------------------------------------------

  # Populate the last column of the matix.
  TC[, dim(TC)[2]] <- TC[1, dim(TC)[2]]

  # Loop though the cells in the matrix and calculate model values.
  # Initialize vector used in for loop.
  TCvector <- vector()
  # Run for loop.
  for (j in 2:dim(TC)[1]) {  
    for (i in 2:(dim(TC)[2] - 1)) {
      TCvector[[1]] <- Ts
      TCvector[[i]] <- TC[j - 1, i] + 
        (K * delta.t.s) * (((TC[j - 1, i + 1] - 
                               (2 * TC[j - 1, i]) + 
                               TC[j - 1, i - 1]) / (delta.x ^ 2)) 
                           + (A / k) 
                           + ((U * (TC[j - 1, i + 1] - 
                                      TC[j - 1, i])) / (K * delta.x)))
    }
    TC[j, c(2:dim(TC)[2] - 1)] <- TCvector
  }
  

## SET UP DATA FRAME FOR PLOTTING ----------------------------------------------  
  
# Pull times of interest from results matrix.

  # Round the Ma time vector.
  t.Ma.round <- round(t.Ma, digits = 2)

  # Initialize data frame.
  TC.DF <- data.frame()
  
  # Build data frame, pulling rows of times of interest from TC.
  TC.DF <- as.data.frame(t(data.frame("t.0"  = TC[(match(0, t.Ma.round)), ], 
                                    "t.0.01" = TC[(match(0.01, t.Ma.round)), ], 
                                    "t.0.1"  = TC[(match(0.1, t.Ma.round)), ], 
                                    "t.1"    = TC[(match(1, t.Ma.round)), ], 
                                    "t.5"    = TC[(match(5, t.Ma.round)), ], 
                                    "t.10"   = TC[(match(10, t.Ma.round)), ], 
                                    "t.20"   = TC[(match(20, t.Ma.round)), ])))
  
  # Rename columns based on depth in km associated with each starting value.
  colnames(TC.DF) <- x / 1000 # TC.DF[1, ]  # could also use starting values
  
  # Add time step column for plotting.
  TC.DF$time.Ma <- c(0, 0.01, 0.1, 1, 5, 10, 20)
  
  # Melt data frame to prepare for plotting.
  TC.melt <- melt(TC.DF, id = "time.Ma")
  
  # Rename columns.
  colnames(TC.melt) <- c("time.Ma", "depth.km", "Temperature.C")
  
  # Sort data frame to make plotting easier.
  TC.melt <-  TC.melt[order(TC.melt$time.Ma),]
  
  # Make depth.km column numeric.
  TC.melt$depth.km <- as.numeric(TC.melt$depth.km)  
  
  # Make time.Ma column into factors.
  TC.melt$time.Ma <- as.factor(TC.melt$time.Ma)  
  

## PLOT GEOTHERMS --------------------------------------------------------------
  
  # Make plot.
  ggplot() +
    
    # Plot results of model.
    geom_path(TC.melt,
              mapping = aes(x = Temperature.C,
                            y = depth.km,
                            color = time.Ma),
              size = 0.75) +
    
    # Make y axis plot in reverse order (0 depth at top of plot).
    scale_y_reverse() +
    
    # Put x axist at top of plot.
    scale_x_continuous(position = "top") +
    
    # Change to simpler ggplot2 theme.
    theme_bw() +
    
    # Change color scale.
    scale_color_brewer(palette = "YlOrRd") +
    
    # Change legend title.
    guides(color = guide_legend(title = "Time (Ma)")) +
    
    # Title plot and label axes.
    labs(title = "Finite difference model geotherms",  # Title plot.
         x = "Temperature (°C)",  # Label x axis.
         y = "Depth (km)")  # Label y axis.
  
