# NRVMs plated on 2015-03-??; stained on 2015-03-13
# cells imaged on 2015-03-14
# cells segmented & data exported on 2015-03-15
# script copied/created on 2015-03-15

# 2015-03-15: moved functions to separate script file
# 2015-03-15: removed points based on histograms

# loads raw data files and makes plots - particle sizes to exclude were
# determined manually by looking at several images from several groups

# desired outputs - cell area & diameter vs. PE dose vs. plating density, mean
# cells per well, number of cells per experimenta condition

# TODO: ANOVA, LM, Box plots, histograms of different experimental groups
# plots to make - 0nm case only; diameter, area, cell # vs. log(dose response)


library(ggplot2)
library(lattice)

rm(list = ls())
graphics.off()

source("150315_functions_v01.R")

# load tab delimited files
raw.all <- read.delim("150315_all_AZ_plates.txt") # all plates

# clean up raw data, remove outliers, and calculate "diameter"
temp.all <- clean.raw(raw.all)

# assign treatments based on Well.name
temp.all[, c("drug.a", "dose.a", "drug.b")] <- assign.treatment.all(temp.all$Plate.ID, temp.all$Well.Name)

temp.all$dose.a <- as.numeric(as.character(temp.all$dose.a))

# calculate summary stats
area <- sum.area(temp.all)
diameter <- sum.diameter(temp.all)

# create data frame for 0nM groups
diameter.zero <- diameter[(diameter$dose.a == 0), ]
area.zero <- area[(area$dose.a == 0), ]

# plot histograms
hist(as.numeric(raw.all$W1.Stained.Area..MultiWaveScoring.), 100)
hist(as.numeric(raw.all$Cell..W2.Stained.Area..MultiWaveScoring.), 100)

hist(temp.all$nucleus, 100)
hist(temp.all$cytoplasm, 100)
hist(temp.all$area, 100)

# plots

# ggplot(data = temp.all, aes(x = dose.a, y = area)) +
#   geom_bar(stat = "identity", position = position_dodge()) + theme(legend.position = "right") +
#   xlab("AZD dose (nM) ") + ylab("Mean area (um^2)") + 
#   ggtitle("Area vs. AZD dose") + theme_bw()
#   geom_errorbar(aes(ymax = area + AI_norm.sem, ymin = AI_norm - AI_norm.sem,), position = "dodge")

# write processed raw data to file to check processing
write.table(temp.all, "all_plates_processed.txt")
write.csv(area, "150315_area.csv")
write.csv(diameter, "150315_diameter.csv")

# ------ old code ------
# 


