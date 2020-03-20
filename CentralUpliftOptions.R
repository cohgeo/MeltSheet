# This script provides code for other ways to elevate the temperature at the
# central uplift under the melt sheet.
# Updated 2020.03.18 CH.

## SET T EVERYWHERE ------------------------------------------------------------
# Set the temperature of the central uplift to be the same everywhere.

# Set the temperature for the uplift.
T.uplift <- 900 + 273.15 # [°C + 273.15 to convert to K] 

for (j in uplift.edge.low:uplift.edge.high) {
  for (i in 2:z.num){
    T.0[j, i] <- T.uplift
  }
}


## QUADRATIC -------------------------------------------------------------------
# Increase the temperature with a quadratic shape (highest T in center, rapid
# decrease in T to edge of central uplift).

# Set the temperature at the center and outer edge of the uplift.
T.uplift <- 900 + 273.15 # [°C + 273.15 to convert to K] 
T.edge <- 500 + 273.15 # [°C + 273.15 to convert to K] 
# Set the diameter of the central uplift to be 10% larger than the melt sheet 
# diameter on either side.
uplift.edge.low <- xstart.melt.index - 0.1 * (meltsheetdiameter / 1000)
uplift.edge.high <- xend.melt.index + 0.1 * (meltsheetdiameter / 1000)

# Set up quadratic.
QA <- (T.edge + T.uplift) / (uplift.edge.low^2 + (median(x) / 1000)^2 +
      (((uplift.edge.low^2 - uplift.edge.high^2) /
      (uplift.edge.high - uplift.edge.low)) * (uplift.edge.low +
      (median(x) / 1000))))
QB <- QA * (uplift.edge.low^2 - uplift.edge.high^2) /
      (uplift.edge.high - uplift.edge.low)
QC <- (QA * (median(x) / 1000)^2) + ((QA * (uplift.edge.low^2 -
      uplift.edge.high^2) * (median(x) / 1000)) / (uplift.edge.high -
      uplift.edge.low)) - (T.uplift)

# Calculate temperatures for central uplift.
for (j in uplift.edge.low:uplift.edge.high) {
  for (i in 2:z.num) {
    T.0[j, ] <- (QA * (x[j] / 1000)^2) + (QB * (x[j] / 1000)) + QC
  }
}
