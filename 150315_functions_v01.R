# originally created 2015-03-15: this file contains functions used in the analysis script
# 2015-03-15: excluded points based on nuclear area, total cell area, and cell.score
# values were chosed based on histograms of the raw.data

# need to make well assignments more intuitive (eg. Excel worksheet notation)

# define minimum and maximum sizes (manually determined for now)
cell.min <- 150
cell.max <- 5000

clean.raw <- function(x) { # where x is the data frame for the raw data
  names(x)[7] <- "area"
  names(x)[6] <- "cell.score"
  names(x)[8] <- "nucleus"
  names(x)[9] <- "cytoplasm"
  #   temp <- subset(x, select = c(Well.Name, Site.ID, Cell.ID, area))
  temp <- x
  temp$area <- as.numeric(as.character(temp$area))
  temp$nucleus <- as.numeric(as.character(temp$nucleus))
  temp$cytoplasm <- as.numeric(as.character(temp$cytoplasm))
  temp$diameter <- 2*sqrt(temp$area/pi)
  temp <- temp[(temp$cell.score == 12), ]
  temp <- temp[!((temp$nucleus > 300) | (temp$area > cell.max)), ]
  temp <- na.omit(temp)
#   temp <- temp[!((temp$area < cell.min) | (temp$area > cell.max)), ]  
  return(temp)
}

assign.treatment.all <- function(plate, well) { # (plate a) where x is the vector of $Well.Name
  temp.a <- rep(0, length(plate)) # dose of drug a
  temp.b <- rep(NA, length(plate)) # name of drug a
  temp.c <- rep("none", length(plate)) # name of drug b
  
  # assign AZD conc. based on columns on the plate
  temp.a[grepl("02", well)] <- 0
  temp.a[grepl("03", well)] <- 10
  temp.a[grepl("04", well)] <- 30
  temp.a[grepl("05", well)] <- 100
  temp.a[grepl("06", well)] <- 300
  temp.a[grepl("07", well)] <- 1000
  temp.a[(grepl("09", well) & grepl("2792", plate))] <- 0
  temp.a[(grepl("10", well) & grepl("2792", plate))] <- 10
  temp.a[(grepl("11", well) & grepl("2792", plate))] <- 30
  temp.a[(grepl("09", well) & grepl("2791", plate))] <- 100
  temp.a[(grepl("10", well) & grepl("2791", plate))] <- 300
  temp.a[(grepl("11", well) & grepl("2791", plate))] <- 1000
  
  # assign AZD  
  temp.b[grepl("2791|2792", plate)] <- "AZD" 
  temp.b[(grepl("02|03|04|05|06|07|", well) & grepl("2793|2795", plate))] <- "AZD"
  
  # assign CSA
  temp.b[(grepl("09", well) & grepl("2795", plate))] <- "CSA"
  temp.b[(grepl("08", well) & grepl("2793", plate))] <- "CSA"
  
  # assign ET-1
  temp.c[(grepl("09|10|11", well) & grepl("2792|2791", plate))] <- "ET-1"
  
  # assign FGF-23
  temp.c[(grepl("02|03|04|05|06|07|08", well) & grepl("2793", plate))] <- "FGF-23"
  
  # assign PE
  temp.c[(grepl("02|03|04|05|06|07", well) & grepl("2792", plate))] <- "PE"
  
  # assign IGF-1
  temp.c[(grepl("02|03|04|05|06|07", well) & grepl("2791", plate))] <- "IGF-1"  
  
  return(as.data.frame(cbind(temp.b, temp.a, temp.c)))
}

sum.area <- function(x) { # where x is the data frame to be analyzed
  temp <- aggregate(area ~ drug.a * dose.a * drug.b, data = x, mean)
  temp$sd <- aggregate(area ~ drug.a * dose.a * drug.b, data = x, sd)$area
  temp$n <- aggregate(area ~ drug.a * dose.a * drug.b, data = x, length)$area
  temp$sem <- temp$sd / sqrt(temp$n)
  return(temp)
}

sum.diameter <- function(x) { # where x is the data frame to be analyzed
  temp <- aggregate(diameter ~ drug.a * dose.a * drug.b, data = x, mean)
  temp$sd <- aggregate(diameter ~ drug.a * dose.a * drug.b, data = x, sd)$diameter
  temp$n <- aggregate(diameter ~ drug.a * dose.a * drug.b, data = x, length)$diameter
  temp$sem <- temp$sd / sqrt(temp$n)
  return(temp)
}

sum.intensity <- function(x) { # where x is the data frame to be analyzed
  temp <- aggregate(average.intensity ~ drug.a * dose.a * drug.b, data = x, mean)
  temp$sd <- aggregate(average.intensity ~ drug.a * dose.a * drug.b, data = x, sd)$average.intensity
  temp$n <- aggregate(average.intensity ~ drug.a * dose.a * drug.b, data = x, length)$average.intensity
  temp$sem <- temp$sd / sqrt(temp$n)
  return(temp)
}