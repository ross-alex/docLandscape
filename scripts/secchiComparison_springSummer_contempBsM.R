## Comparison between spring and summer secchi depths from the same sampling year
# AJR: 2024-03-20

library(tidyverse)
library(here)

# Read data
secDat <- read_csv(here("data","summerSecchis_bsm.csv")) %>% 
  select(WbyName,`Water Chem Sample Date`, `Water Chem Secchi Depth (Spring)`,`Limno Sample Date`, `Limno Secchi Depth (Summer)`) %>% 
  rename(waterbodyID = WbyName, 
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
# in summer predictions?)