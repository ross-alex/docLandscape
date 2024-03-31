## Comparison between spring and summer secchi depths from the same sampling year
# AJR: 2024-03-20

library(tidyverse)
library(here)

# Read data
secDat <- read_csv(here("data","summerSecchis_bsm.csv")) %>% 
  select(WbyLID,`Water Chem Sample Date`, `Water Chem Secchi Depth (Spring)`,`Limno Sample Date`, `Limno Secchi Depth (Summer)`) %>% 
  rename(waterbodyID = WbyLID, 
         waterChem_dateSample = `Water Chem Sample Date`, 
         limno_dateSample = `Limno Sample Date`,
         springSecchi = `Water Chem Secchi Depth (Spring)`, 
         summerSecchi = `Limno Secchi Depth (Summer)`) %>% 
  relocate(limno_dateSample, .after = waterChem_dateSample)

# Remove rows with NAs and where waterChem_dateSample == limno_dateSample
secDat_complete <- secDat %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  filter(waterChem_dateSample != limno_dateSample)

# T-test between spring and summer secchi's
secTest <- t.test(secDat_complete$springSecchi,secDat_complete$summerSecchi)
secTest
dim(secDat_complete)

# Visualize data, boxplot
boxplot(secDat_complete$springSecchi,secDat_complete$summerSecchi, names = c("Spring","Summer"),
        ylab="Secchi Depth", las = 1)

# Summer secchi's are slightly higher (mean spring = 3.37, mean summer = 3.7), so lakes appear clearer during the summer
# We could add a correction factor to ARU Secchi's (add 0.33cm to BsM data to translate spring model data for use
# in summer predictions?). 
#** Should actually correct ARU "summer" data to "spring" data if a correction is to be made (remove 33cm)

# It should be noted, these are for ALL BsM lakes, and have not been filtered to the model data

## Only lakes in our contemporary dataset
statDat_wide <- read_csv(here("data","BsM_WideSampleModelData_20240311.csv")) %>% #This dataset is slightly reduced
  filter(log(TKN) >4) #outliers <4                            #to remove all rows that had at 
                                                              #at least one NA. See 'modelData_assembly_20240311.R'
                                                              # script for details (bottom chunk)

#Pull out waterbodyIDs from our study
bsmUnique <- tibble(statDat_wide$waterbodyID) 
colnames(bsmUnique) <- "waterbodyID"

# Reduce data set and re-test
# Remove rows with NAs and where waterChem_dateSample == limno_dateSample
secDat_reduced <- secDat_complete %>% 
  semi_join(bsmUnique)

# T-test between spring and summer secchi's
secTest_red <- t.test(secDat_reduced$springSecchi,secDat_reduced$summerSecchi)
secTest_red
dim(secDat_reduced)

# Visualize data, boxplot
boxplot(secDat_reduced$springSecchi,secDat_reduced$summerSecchi, names = c("Spring","Summer"),
        ylab="Secchi Depth", las = 1)

## Same results; significant test and same mean difference (0.32 compared to 0.33)
# in summer predictions?)