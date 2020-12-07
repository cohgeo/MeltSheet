# This script explores how changing the liquidus and solidus changes the heat
# capacity.
# This code is designed to run in R.
# Updated 2020.12.07 CH

## SETUP -----------------------------------------------------------------------

# Clear all from workspace/environment.
# rm(list = ls())


## SET MODEL PARAMETERS --------------------------------------------------------

# All temperatures in °C should be converted to K before running the model.

# Set latent heat of fusion of diopside.
  # Suggested value: 421000 J/kg; Abramov and Kring (2007)
  L <- 421000  # [J/kg] 
  
# Import table of solidus, Zr sat., and liquids temperatures from MELTS.
  MELT.T <- read.csv("Morokweng Results/MeltSheet Results_2020.12.07/Melts_All_Runs.csv")
  
# Possible liquidus temperatures.
  # Suggested value: 1195°C from MELTS modeling (used in model).
  T.liq.min <- min(MELT.T$T.liquidus.C)   # [°C]
  T.liq.max <- max(MELT.T$T.liquidus.C)   # [°C]
  
# Possible solidus temperatures.
  # Suggested value: 796°C from MELTS modeling (used in model).
  T.sol.min <- min(MELT.T$T.solidus.C)   # [°C]
  T.sol.max <- max(MELT.T$T.solidus.C)   # [°C]
  
# Possible pairing for maximum and minimum ranges between liquidus and solidus.
  T.range.max <- T.liq.max - T.sol.min   # [°C]
  T.range.min <- T.liq.min - T.sol.max   # [°C]
  
# Set heat capacity with no latent heat.
  # Suggested value: 1000 J/(kg/K), heat capacity of basement, melt, breccia; 
  # Abramov and Kring (2007)
  Cp.nolatent <- 1000  # [J/(kg*K)]
  # Set initial heat capacity.
  Cp.initial <- Cp.nolatent
  
# Calculate heat capacity accounting for latent heat of fusion.
  # Suggested value: 3339 J/(kg*K); Abramov and Kring (2007)
  # Suggested value: 2055 J/(kg*K); based on T.liquidus = 1195°C and 
  # T.solidus = 796°C from MELTS modeling.
  Cp.prime.max <- Cp.initial + (L / T.range.min)  # [J/(kg*K)]
  Cp.prime.min <- Cp.initial + (L / T.range.max)  # [J/(kg*K)]
  Cp.prime.model <- Cp.initial + (L / (1195 - 796))  # [J/(kg*K)]
  
# Make plot of range 
  # Initialize vector to hold range between T.liquidus and T.solidus.
  Range <- 0
  # Create sequence from min to max range.
  Range <- seq(from = T.range.min, to = T.range.max, by = 1)
  # Initialize Cp.prime vector.
  Cp.prime <- 0
  # Calculate Cp.prime for each value in Range vector.
  for (i in 1:length(Range)){
    Cp.prime[[i]] <- Cp.initial + (L / Range[[i]]) 
  }

# Plot possible Cp.prime values.
  plot(x = Range, y = Cp.prime,
       type = "line")
  points(x = (1195 - 796), y = Cp.prime.model)
  