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
  filter(waterChem_dateSample != limno_dateSample) %>% 
  mutate(secDatDiff = springSecchi-summerSecchi)

# T-test between spring and summer secchi's
secTest <- t.test(secDat_complete$springSecchi,secDat_complete$summerSecchi, paired = TRUE)
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

## Only lakes in our contemporary dataset                                                           # script for details (bottom chunk)
mainDf <- read_csv("data/mainAnalyticalDataset_20240716.csv")

#Pull out waterbodyIDs from our study
bsmUnique <- tibble(mainDf$waterbodyID) %>% 
  distinct()
colnames(bsmUnique) <- "waterbodyID"

# Reduce data set and re-test
# Remove rows with NAs and where waterChem_dateSample == limno_dateSample
secDat_reduced <- secDat_complete %>% 
  semi_join(bsmUnique)

# T-test between spring and summer secchi's
secTest_red <- t.test(secDat_reduced$springSecchi,secDat_reduced$summerSecchi, paired=TRUE)
secTest_red
dim(secDat_reduced)

# Visualize data, boxplot
boxplot(secDat_reduced$springSecchi,secDat_reduced$summerSecchi, names = c("Spring","Summer"),
        ylab="Secchi Depth", las = 1)

## Significant test and mean difference of 0.37 (summer is clearer)
range(secDat_reduced$secDatDiff)
0.37/17.2 # The overall mean difference between spring and summer secchi depths are 2.2% of the range of overall secchi differences
# between paired spring/summer secchi observations, so despite the significant test, the degree of difference between
# secchi readings is small relative to the among-lake variations. Thus, we chose to use Spring secchi's because there was
# (1) more data available, (2) a shorter time-window of when sampling was taken making between-lake comparisons more 
# relevant to one another, (3) taken at the same time in which water chemsitry (e.g., DOC) was sampled, and (4) taken
# by a consistent group of personel which significantly reduces possible variation caused by differences in Secchi observers


## Try a chi square test - are higher spring/summer more frequent?
# Step 1: Create a new column to compare spring and summer Secchi values
secDat_comp <- secDat_reduced %>%
  mutate(Comparison = case_when(
    springSecchi > summerSecchi ~ "Spring Higher",
    springSecchi < summerSecchi ~ "Summer Higher",
    TRUE ~ "Equal"
  ))

# Step 2: Create a contingency table
contingency_table <- table(secDat_comp$Comparison)

# Step 3: Perform the chi-square test
chi_square_result <- chisq.test(contingency_table)

# Print the results of the chi-square test
print(chi_square_result)

## What does a 0.3 m difference in secchi do for predicted DOC? ####
# The model used to estimate DOC was:
bestmodel2 <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
                      data = mainDF_BsM, family = Gamma(link = "log"))

predictedDOC_preds_dfSum <- mainDf %>% 
  summarise(meanSecchi = mean(secchiDepth),
            meanMaxDepth = mean(maxDepth), 
            meanLat = mean(lat))

# Now get estimates from other script for model coefficients 
docLinearPred <- -2.8086626 + (-0.0795200*-10.3) # Since scale = F all estimates are at 0. Essentially 
# calculating the intercept +/- mean secchiDepth (which is 0) 
exp(docLinearPred) # a 0.3m change in secchi at mean variables is only worth 0.0014 mg/L DOC

# Coefficients from your model
intercept <- -2.8086626
beta_secchi <- -0.0795200  

# Change in Secchi depth (centered value)
delta_secchi <- 0.3

# Calculate the change in log-transformed DOC
delta_log_DOC <- beta_secchi * delta_secchi

# Calculate the multiplicative change in DOC
multiplicative_change_DOC <- exp(delta_log_DOC)

# Calculate the predicted DOC concentration after the change in Secchi depth
# Assuming the original DOC is based on the intercept (when all predictors are at their mean)
original_DOC <- exp(intercept)  # The original DOC concentration when Secchi, MaxDepth, and Latitude are at their means

# New DOC concentration after 0.3m change in Secchi depth
new_DOC <- multiplicative_change_DOC * original_DOC

# Output the results
cat("Original DOC concentration:", original_DOC, "mg/L\n")
cat("New DOC concentration after 0.3m Secchi depth change:", new_DOC, "mg/L\n")
cat("Difference that 0.3m Secchi affects predicted DOC at predictor means:", 
    original_DOC - new_DOC, "mg/L\n")

exp(-2.81 + (-0.080*0.3 - 3.55) + (-0.0059*0 - 29.56) + (0.098*0) + (0.0003*(0.3 - 3.55)*(0 - 29.56)))


