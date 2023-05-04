#Code to plot pedigrees using kinship2
#Juliana Acosta-Uribe 2022

#1. Install and load packages
install.packages("kinship2")
library(kinship2)

#2. Set your working directory (the folder where you have your documents)
setwd("~/Documents/Kosik Lab -UCSB/ReD-Lat/Training/Week_8")

#3. Load your data
Family_data <- read.delim("Pedigree.txt")

#4. Create the pedigree object
Family_ped <- pedigree(Family_data$id, 
                         Family_data$dadid,
                         Family_data$momid,
                         Family_data$sex,
                         Family_data$redhair)

#5. Plot your pedigree to check its ok
print.pedigree(Family_ped)
  
plot(Family_ped, cex = 0.7)

#6. Save it as an eps file
setEPS()
postscript("Family_ped.eps")
plot(Family_ped)
dev.off()

## Now you can open the "Family_ped.eps" 
