# NRVMs plated on 2015-03-??; stained on 2015-03-13
# cells imaged on 2015-03-14
# cells segmented & data exported on 2015-03-15
# script copied/created on 2015-03-15

# loads raw data files and makes plots - particle sizes to exclude were
# determined manually by looking at several images from several groups

# desired outputs - cell area & diameter vs. PE dose vs. plating density, mean
# cells per well, number of cells per experimenta condition

# TODO: ANOVA, LM, Box plots, histograms of different experimental groups
#

library(ggplot2)

rm(list = ls())

# define minimum and maximum sizes (manually determined for now)
cell.min <- 150
cell.max <- 1000

clean.raw <- function(x) { # where x is the data frame for the raw data
  names(x)[6] <- "area"
#   temp <- subset(x, select = c(Well.Name, Site.ID, Cell.ID, area))
  temp <- x
  temp$area <- as.numeric(temp$area)
  temp$diameter <- sqrt(temp$area/pi)
  temp <- temp[!((temp$area < cell.min) | (temp$area > cell.max)), ]
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
  temp.a[(grepl("09", well) & grepl("2768", plate))] <- 0
  temp.a[(grepl("10", well) & grepl("2768", plate))] <- 10
  temp.a[(grepl("11", well) & grepl("2768", plate))] <- 30
  temp.a[(grepl("09", well) & grepl("2767", plate))] <- 100
  temp.a[(grepl("10", well) & grepl("2767", plate))] <- 300
  temp.a[(grepl("11", well) & grepl("2767", plate))] <- 1000
  
  # assign AZD  
  temp.b[grepl("2767|2768", plate)] <- "AZD" 
  temp.b[(grepl("02|03|04|05|06|07|", well) & grepl("2769|2770", plate))] <- "AZD"
  
  # assign CSA
  temp.b[(grepl("09", well) & grepl("2770", plate))] <- "CSA"
  temp.b[(grepl("08", well) & grepl("2769", plate))] <- "CSA"
  
  # assign ET-1
  temp.c[(grepl("09|10|11", well) & grepl("2768|2767", plate))] <- "ET-1"
  
  # assign FGF-23
  temp.c[(grepl("02|03|04|05|06|07|08", well) & grepl("2769", plate))] <- "FGF-23"
  
  # assign PE
  temp.c[(grepl("02|03|04|05|06|07", well) & grepl("2768", plate))] <- "PE"
  
  # assign IGF-1
  temp.c[(grepl("02|03|04|05|06|07", well) & grepl("2767", plate))] <- "IGF-1"  
  
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

# load tab delimited files
raw.all <- read.delim("150306_all_plates.txt") # all plates

# clean up raw data, remove outliers, and calculate "diameter"
temp.all <- clean.raw(raw.all)

# assign treatments based on Well.name
temp.all[, c("drug.a", "dose.a", "drug.b")] <- assign.treatment.all(temp.all$Plate.ID, temp.all$Well.Name)

temp.all$dose.a <- as.numeric(as.character(temp.all$dose.a))

# calculate summary stats
area <- sum.area(temp.all)
diameter <- sum.diameter(temp.all)

# plots

ggplot(data = temp.all, aes(x = dose.a, y = area)) +
  geom_bar(stat = "identity", position = position_dodge()) + theme(legend.position = "right") +
  xlab("AZD dose (nM) ") + ylab("Mean area (um^2)") + 
  ggtitle("Area vs. AZD dose") + theme_bw()
#   geom_errorbar(aes(ymax = area + AI_norm.sem, ymin = AI_norm - AI_norm.sem,), position = "dodge")

# write processed raw data to file to check processing
write.table(temp.all, "all_plates_processed.txt")

# ------ old code ------
# 
# assign.treatment.b <- function(x) { # (plate b) where x is the vector of $Well.Name
#   temp.a <- rep(NA, length(x)) # dose of drug a
#   temp.b <- rep(NA, length(x)) # name of drug a
#   temp.c <- rep(NA, length(x)) # name of drug b
#   temp.a[grepl("02|09", x)] <- 0 # assign AZD conc. based on columns on the plate
#   temp.a[grepl("03|10", x)] <- 10
#   temp.a[grepl("04|11", x)] <- 30
#   temp.a[grepl("05", x)] <- 100
#   temp.a[grepl("06", x)] <- 300
#   temp.a[grepl("07", x)] <- 1000
#   temp.b <- rep("AZD", length(x)) # AZD used in all wells
#   temp.c <- rep("PE", length(x)) # initialize with PE
#   temp.c[grepl("09|10|11", x)] <- "ET1" # assign columns 9-11 to ET-1
#   return(as.data.frame(cbind(temp.b, temp.a, temp.c)))
# }
# 
# assign.treatment.c <- function(x) { # (plate c) where x is the vector of $Well.Name
#   temp.a <- rep(NA, length(x)) # dose of drug a
#   temp.b <- rep(NA, length(x)) # name of drug a
#   temp.c <- rep(NA, length(x)) # name of drug b
#   temp.a[grepl("02|08", x)] <- 0 # assign AZD conc. based on columns on the plate
#   temp.a[grepl("03", x)] <- 10
#   temp.a[grepl("04", x)] <- 30
#   temp.a[grepl("05", x)] <- 100
#   temp.a[grepl("06", x)] <- 300
#   temp.a[grepl("07", x)] <- 1000
#   temp.b <- rep("AZD", length(x)) # initialize with AZD
#   temp.b[grepl("08", x)] <- "CSA"
#   temp.c <- rep("FGF", length(x)) # initialize with IGF-1
#   return(as.data.frame(cbind(temp.b, temp.a, temp.c)))
# }
# 
# assign.treatment.d <- function(x) { # (plate d) where x is the vector of $Well.Name
#   temp.a <- rep(NA, length(x)) # dose of drug a
#   temp.b <- rep(NA, length(x)) # name of drug a
#   temp.c <- rep(NA, length(x)) # name of drug b
#   temp.a[grepl("02|09", x)] <- 0 # assign AZD conc. based on columns on the plate
#   temp.a[grepl("03", x)] <- 10
#   temp.a[grepl("04", x)] <- 30
#   temp.a[grepl("05", x)] <- 100
#   temp.a[grepl("06", x)] <- 300
#   temp.a[grepl("07", x)] <- 1000
#   temp.b <- rep("AZD", length(x)) # initialize with AZD
#   temp.b[grepl("09", x)] <- "CSA"
#   temp.c <- rep("None", length(x)) # initialize with NA
#   return(as.data.frame(cbind(temp.b, temp.a, temp.c)))
# }


