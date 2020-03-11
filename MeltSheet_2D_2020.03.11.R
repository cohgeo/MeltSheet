# This script creates a finite difference model of the cooling of a melt sheet
# caused by an impact.
# Updated 2020.03.11 CH

## SETUP -----------------------------------------------------------------------

# Clear all from workspace/environment.
rm(list = ls())

# Install required packages.


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

```{r Constants}
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
                     nrow = z.num,  
                     ncol = x.num)

Cp.prime <- Cp.initial + (L / (T.liquidus - T.solidus)) 
# heat capacity accounting for latent heat of fusion, [J/(kg*K)]
# Suggested value: 3339 J/(kg*K); Abramov and Kring (2007)
Cp.0 <- Cp.initial  # initialize matrix for initial model setup below

k <- matrix(data = 2.5,    # thermal conductivity, [W/(m*K)]
            nrow = z.num,  
            ncol = x.num)  
# Suggested value: 2.5 W/(m*K), thermal conductivity of basement, melt, 
# breccia; Abramov and Kring (2007)

A <- matrix(data = 2.2e-6,      # heat production, [W/(m^3)]
            nrow = z.num,  
            ncol = x.num)  
# Average A value from Jones (1988), table 6

rho <- matrix(data = 2700,   # density, [kg/(m^3)]
              nrow = z.num,  
              ncol = x.num)
# Suggested value: 2700 kg/(m^3), density of basement, melt, breccia; 
# Abramov and Kring (2007)

Kappa <- k / (rho * Cp.nolatent)  # thermal diffusivity, [(m^2)/s]
Kappa.min <- min(Kappa)      # [(m^2)/s] 

U <- matrix(data = 0,        # uplift, [m/s]
            nrow = z.num,
            ncol = x.num)
# We assumed no tectonic uplift. 

# Input geothermal gradient in K/km.
a.Kkm <- 16.5  # [K/km]
# Average value between 12 and 21 K/km from Jones (1992)
a <- a.Kkm / 1000  # [K/m]

# Define time step.
dt.max <- (dx ^ 2) / (6 * Kappa.min)  # [s] 
dt.yr <- dt.max * 3.17098e-8 # [yr]
dt <- dt.max
```

```{r InitializeTDistribution_plume}
# Set up an initial condition of a melt sheet at an elevated temperature.

# Initialize a matrix for initial temperature conditions.
T.0 <- matrix(data = 0,
              nrow = z.num,
              ncol = x.num)

# Set the top row to be air temperature.
# Choose a value for the air temperature.
T.surface <- 25  # [°C]
# Assign the first row of the initial matrix to be the air temperature.
T.0[1, ] <- T.surface + 273.15

# Make a geotherm for the initialization matrix.
for (j in 2:dim(T.0)[1]) {
  T.0[j, ] <- T.0[1, 1] + a * x[j]
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

# Emplace melt sheet. Set material properties of melt sheet.
for (i in zstart.melt.index:zend.melt.index) {
  for (j in xstart.melt.index:xend.melt.index) {
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

# Clean up.
rm(i, j)
```

### Initial Condtions  

**Model Parameters:**  
  dt = `r round((dt * 3.17098e-8 / 1000000), digits = 3)` Myr  
dx = `r round(dx, digits = 1)` m  
dz = `r round(dz, digits = 1)` m  

**Initial perterbation:**  
  rectangular intrusion with `r (meltsheetdiameter / 2 / 1000)` km radius and `r (meltsheetthickness / 1000)` km thickness at an elevated temperature of `r T.meltsheet`°C intrudes into rock with a geothermal gradient of `r a.Kkm` °C/km

```{r PlotInitialConditions_heatmap}
# Plot initial conditions.

# Load packages plotly and Rcolorbrewer for plotting. (GGK added next four lines on 02July2019 as issue with library(plotly) wiht my paths)
# myPaths <- .libPaths()
# myPaths <- c(myPaths, "C:/R stuff")
# .libPaths(myPaths)

library(ggplot2)  
library(plotly)
library(RColorBrewer)

# Make function for plotting.
PlotlyMeltSheet <- function(temperatures) {
  # Make matrix for hovertext labels.
  text.matrix <- round(matrix(unlist(temperatures), 
                              ncol = dim(temperatures)[2], 
                              byrow = FALSE) - 273.15, 
                       digits = 0)
  
  # Make plot.
  plot_ly(x = ~(x / 1000),
          y = ~(z / 1000), 
          z = ~(temperatures - 273.15),
          type = "heatmap",
          colors = rev(brewer.pal(9,"RdYlBu")),
          hoverinfo = "text",
          text = text.matrix) %>%
    add_trace(y = ((z[[2]] - z[[1]])/ 2 / 1000),
              type = "scatter",
              mode = "lines",
              line = list(color = "black")) %>%
    colorbar(title = "Temperature (°C)",
             limits = c(0, 2000))  %>%
    layout(xaxis = list(title = "Distance (km)"),
           yaxis = list(title = "Depth (km)",
                        autorange = "reversed"))
}

# Plot initial condtions as a heatmap.
PlotlyMeltSheet(T.0)
```

```{r StaticPlot}
# Plot a matrix as a heatmap with ggplot.
# Load required packages.
library(reshape2)
library(ggplot2)

# Make plot
ggplot(melt(T.0), 
       aes(Var1,
           Var2,
           fill = value)) + 
  geom_raster()

# Plot a matrix as a heatmap with base R.
# heatmap(T.0) # Doesn't work

# Plot a matrix as a heatmap with the package lattice.
install.packages("lattice")
library(lattice)
levelplot(T.0)
```



```{r Model_SetIterations}
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
                           nrow = z.num,
                           ncol = x.num)
    
    # Initialize Cp[[t + 1]]
    Cp[[t + 1]] <- matrix(data = NA,
                          nrow = z.num,
                          ncol = x.num)
    
    # Calculate T for interior of T.n[[t + 1]] matrix.
    for (i in 2:(x.num - 1)) {
      for (j in 2:(z.num - 1)) {
        
        # Set constant T border.
        T.n[[t + 1]][1, ] <- T.n[[t]][1, ]          # top
        T.n[[t + 1]][z.num, ] <- T.n[[t]][z.num, ]  # bottom
        T.n[[t + 1]][, 1] <- T.n[[t]][, 1]          # left
        T.n[[t + 1]][, x.num] <- T.n[[t]][, x.num]  # right
        
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

# Solve the finite difference equation for n time steps.
n <- 2
model.results <- RunModel(n)
T.n <- model.results$T.n
Cp <- model.results$Cp
```

### Results after `r eval(n)` iterations

Maximum temperature change between model run n = `r n` (time = `r round((n * dt * 3.17098e-8 / 1000000), digits = 2)` Myr) and initial condition (time = 0 Myr):  
  `r  round(((max(T.n[[length(T.n)]] - T.n[[1]])) - 273.15), digits = 2)` °C   

```{r PlotModel_heatmap}
# Plot model results after n iterations.
PlotlyMeltSheet(T.n[[length(T.n)]])
```

```{r Setup_m_iterations}
# Solve the finite difference equation for m time steps.
m <- 5
T.m <- RunModel(m)
```

### Results after `r eval(m)` iterations

Maximum temperature change between model run n = `r m` (time = `r round((m * dt * 3.17098e-8 / 1000000), digits = 2)` Myr) and initial condition (time = 0 Myr):  
  `r  round(((max(T.m[[length(T.m)]] - T.m[[1]])) - 273.15), digits = 2)` °C   

```{r PlotModel_heatmap_m}
# Plot model results after m iterations.
PlotlyMeltSheet(T.m[[length(T.m)]])
```


```{r Setup_o_iterations}
# Solve the finite difference equation for m time steps.
o <- 9
T.o <- RunModel(o)
```

### Results after `r eval(o)` iterations

Maximum temperature change between model run n = `r o` (time = `r round((o * dt * 3.17098e-8 / 1000000), digits = 2)` Myr) and initial condition (time = 0 Myr):  
  `r  round(((max(T.o[[length(T.o)]] - T.o[[1]])) - 273.15), digits = 2)` °C   

```{r PlotModel_heatmap_o}
# Plot model results after m iterations.
PlotlyMeltSheet(T.o[[length(T.o)]])
```