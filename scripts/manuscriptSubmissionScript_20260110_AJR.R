## Last script update. This reflects results as they're presented in the MS
# AJR: 2024-09-01

#**ALEX, note from March 10,2025 - maybe add elevation into predictive model. Easy to get info
#* otherwise, maybe add watershed slope (less easy, but useful predictor)

# Load libraries
library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(MuMIn)
library(visreg)
library(patchwork)
library(svglite)
library(sf)
library(ghibli)
source("scripts/theme_DOC.R") #load ggplot theme


## Read Data ####
combIDs <- read_csv(here("data/commonDataIDs_aruToBsm_20240307.csv")) #has common IDs between BsM and ARU
lakeElevationDat <- read_csv("data/BsM_ARU_elevation_20240320.csv")  ## Extract elevation data
statDat_wide <- read_csv("data/statDatWide_updatedCyc3_20240520.csv") %>% 
  filter(log(TKN) >4) %>%  # Clear outliers < 4 (only two observations)                       
  left_join(select(combIDs, bsmWaterbodyID, commonID), by = c("waterbodyID" = "bsmWaterbodyID"))
# Ontario shapefile: 
ontario_shapefile <- st_read("data/Province/Province.shp")
ontario_shapefile <-st_transform(ontario_shapefile, crs = 4326)
ontario_ohnWaterbody_latlon <- read_csv("data/ohnWaterbodies_csvWithCoords.csv") %>% 
  mutate(geoRegion = if_else(Latitude > 47, "North (> 47° N)", "South (< 47° N)"))

## All data (AHI and BsM)
mainDF <- read_csv("data/dataForDryadRepo_landscapeDOC_AR_20241122.csv") %>% 
  left_join(select(combIDs, bsmWaterbodyID, commonID), by = "commonID") %>% 
  left_join(select(lakeElevationDat,bsmWaterbodyID,elevation), by="bsmWaterbodyID")

# Create BsM only DF
mainDF_onlyBsM <- mainDF %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE))) %>% 
  arrange(year)

mainDF_onlyBsM$lakeTrophicStatus <- factor(mainDF_onlyBsM$lakeTrophicStatus, 
                                           levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting


# Correlation between DOC and Secchi, and TDP ####
cor.test(log(mainDF_onlyBsM$secchiDepth),log(mainDF_onlyBsM$DOC))

pCor_secchiDOC <- ggplot(mainDF_onlyBsM) +
  geom_point(aes(x=secchiDepth,y=DOC), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_smooth(aes(x=(secchiDepth),y=(DOC)), method = "lm") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  scale_x_continuous(trans='log', labels = function(x) signif(x, 2)) +
  scale_y_continuous(trans='log', labels = function(x) signif(x, 2)) +
  ylab("DOC (mg/L)") +
  xlab("Secchi depth (m)") +
  # ylim(c(0,21))+
  # xlim(c(0,21)) +
  theme_DOC()
pCor_secchiDOC

# ggsave(filename = ("plotsTables/plots/correlationSecchiDOC_20260121.svg"),device = "svg", plot = pCor_secchiDOC,
#        width = 8, height = 6,bg="white")

# Correlation between DOC and TDP 
cor.test(log(mainDF_onlyBsM$TDP),log(mainDF_onlyBsM$DOC))

pCor_TDPDOC <- ggplot(mainDF_onlyBsM) +
  geom_point(aes(x=TDP,y=DOC), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_smooth(aes(x=(TDP),y=(DOC)), method = "lm") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  scale_x_continuous(trans='log', labels = function(x) signif(x, 3)) +
  scale_y_continuous(trans='log', labels = function(x) signif(x, 3)) +
  ylab("DOC (mg/L)") +
  xlab("Total dissolved phosphorus (ug/L)") +
  # ylim(c(0,21))+
  # xlim(c(0,21)) +
  theme_DOC() 
pCor_TDPDOC

environCovariateCorP <- pCor_secchiDOC + pCor_TDPDOC + plot_layout(ncol = 2)

# ggsave(filename = ("plotsTables/plots/environCovariateCor_20260121.svg"),device = "svg", plot = environCovariateCorP,
#        width = 8, height = 6,bg="white")


## Change in DOC and Secchi; ONLY BsM ####
## DOC
# Using IC approach, assemble different models to look for changes in DOC between 2008-2023
M1a <- glmmTMB(DOC ~ scaled_yearSample + (1|commonID),
               data = mainDF_onlyBsM,
               family = Gamma(link = "log"))
M2a <- glmmTMB(DOC ~ scaled_yearSample + (scaled_yearSample|commonID),
               data = mainDF_onlyBsM,
               family = Gamma(link = "log"))          
M3a = update(M1a, family = lognormal(link="log"))
M4a = update(M2a, family = lognormal(link="log"))
AIC(M1a,M2a,M3a,M4a)
# Gamma models are best; random intercept marginally better than random slopes

# Parameterize best model for DOC over time
bestmodel3a = M1a

# Diagnostics
DHARMa::simulateResiduals(bestmodel3a, plot=T)
hist(residuals(bestmodel3a), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumDOCa <- summary(bestmodel3a)) 
exp(sumDOCa$coefficients$cond[, "Estimate"]) #DOC decrease by 0.0032% per year; NS
r.squaredGLMM(bestmodel3a)          # pretty much all lake effect

mainDF_onlyBsM$DOC_fit = predict( bestmodel3a, type = "response", re.form = NA)
mainDF_onlyBsM$DOC_fitWRandom = predict( bestmodel3a, type = "response", re.form = NULL)  # To plot just predicted lines

# Plot it:
pDOCTime_bsm <- ggplot(data=mainDF_onlyBsM,aes(x=yearSample,y=DOC)) +
  geom_jitter(aes(x=yearSample, y = DOC, fill=DOC), alpha=0.7, size=2.5, shape=21, width = 0.25) +
  #geom_smooth(method = "lm") +      # Don't plot line as relationship not significant
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  labs(y="DOC (mg/L)",x="Year") +
  theme_DOC()+
  theme(legend.position = "none")
pDOCTime_bsm

# ggsave(filename = ("plotsTables/plots/DOCTime_bsm_20260121.svg"),device = "svg", plot = pDOCTime_bsm,
#        width = 8, height = 6)

## Secchi - do we see consistent trends with DOC?
secchiTimeMod <- glmmTMB(secchiDepth ~ scaled_yearSample+ (1|commonID), 
                         data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
secchiTimeModb <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
                          data = mainDF_onlyBsM,
                          family = Gamma(link = "log"))          
secchiTimeModc = update(secchiTimeMod, family = lognormal(link="log"))
secchiTimeModd = update(secchiTimeModb, family = lognormal(link="log"))
AIC(secchiTimeMod,secchiTimeModb,secchiTimeModc,secchiTimeModd)

secchiBsmMod_year <- secchiTimeModb

# Diagnostics
DHARMa::simulateResiduals(secchiBsmMod_year, plot=T)
hist(residuals(secchiBsmMod_year), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumSecchiYear <- summary(secchiBsmMod_year)) 
exp(sumSecchiYear$coefficients$cond[, "Estimate"]) #Secchi decrease by 1.27% per year; significant
r.squaredGLMM(secchiBsmMod_year)          # pretty much all lake effect

mainDF_onlyBsM$Secchi_fit = predict(secchiBsmMod_year, type = "response", re.form = NA)
mainDF_onlyBsM$Secchi_fitWRandom = predict(secchiBsmMod_year, type = "response", re.form = NULL, se.fit = TRUE)  # To plot just predicted lines

secchiBsmYearPlot <- ggplot(mainDF_onlyBsM)+
  geom_jitter(aes(x=yearSample, y = secchiDepth, fill=secchiDepth), alpha=0.7, size=2.5, shape=21, width = 0.25)+
  geom_smooth(aes(x=yearSample, y = secchiDepth),    # Use raw data, create a model - no partial effects as only one predictor (all marginal anyways)
              colour = "blue",
              method="glm",
              formula = y~x,
              se = TRUE,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  theme(legend.position = "none")
secchiBsmYearPlot

# ggsave(filename = ("plotsTables/plots/secchiTime_bsm_20260121.svg"),device = "svg", plot = secchiBsmYearPlot,
#        width = 8, height = 6)


## Plot random slopes associated w. Secchi year model
#Add fitted values from the mixed model (conditional = include random effects)
mainDF_onlyBsM$fit_cond_secchiYear <- predict(
  secchiTimeModb,
  type = "response",   # get on original gamma scale
  re.form = NULL       # include random effects
)

#Remove instances where lakes were only sampled once to emphasize lakes w. repeated samples
plotDF_secYearBsM <- mainDF_onlyBsM %>%
  group_by(commonID) %>%
  filter(n() > 1) %>%
  ungroup()

#Plot lake-specific fitted lines
pSecchiBsM_randomSlopes <- ggplot(plotDF_secYearBsM, aes(x = yearSample, y = Secchi_fitWRandom)) +
  geom_point(aes(fill = Secchi_fitWRandom),
            alpha = 0.7, size = 2.5, shape = 21) +
  geom_line(aes(y = fit_cond_secchiYear, colour = commonID), alpha=0.7) +   # <- one slope per lake
  scale_x_continuous("Year") +
  scale_y_continuous("Secchi depth (m)") +
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  scale_colour_viridis_d() +
  theme_DOC() +
  theme(legend.position = "none")
pSecchiBsM_randomSlopes

# ggsave(filename = ("plotsTables/plots/secchiTime_bsm_randomSlopes_20260121.svg"),device = "svg", plot = pSecchiBsM_randomSlopes,
#        width = 8, height = 6)


# Summarize modelled differences in Secchi over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_secNoEnv <- ggplot_build(secchiBsmYearPlot)
# Extract the data related to geom_smooth
lineData_secNoEnv <- plotData_secNoEnv$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumsecNoEnv_timeBsM <- lineData_secNoEnv %>%
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)


## Change in DOC and Secchi with only BsM
docBsMTimeDat <- mainDF_onlyBsM %>% 
  group_by(commonID) %>% 
  summarise(initDOC = DOC[which.min(yearSample)],
            finalDOC = DOC[which.max(yearSample)],
            diffDOC = finalDOC-initDOC,
            yearDiffDOC = max(yearSample) - min(yearSample),
            diffDOCperYear = diffDOC/yearDiffDOC,
            initSecchi = secchiDepth[which.min(yearSample)],
            finalSecchi = secchiDepth[which.max(yearSample)],
            diffSecchi = finalSecchi-initSecchi,
            diffSecchiperYear = diffSecchi/yearDiffDOC,
            sampleSize = n()) %>% 
  left_join(select(statDat_wide,commonID,lat,long)) %>%
  distinct() %>%
  filter(!is.na(lat)& yearDiffDOC >1) %>%
  mutate(dataset = "all")
docBsMTimeDat

## Make maps of deltaDOC and deltaSecchi
## Map of DOC, only BsM
pDOCtimeDiff_BsM <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  #geom_point(data = ontario_ohnWaterbody_latlon, 
  #           aes(x=Longitude,y=Latitude),size=1,alpha=0.1) +
  geom_point(data=docBsMTimeDat, aes(y=lat,x=long,fill=diffDOCperYear), size=4, 
              shape = 21, stroke = 0.5, colour = "black", alpha=0.8) +
  scale_fill_gradient2(
    low = "navyblue",
    mid = "white",
    high = "brown",
    midpoint = 0)+
  ylim(c(42.5,57)) +
  labs(y="Latitude",x="Longitude",fill="Δ DOC (mg/L/year)") +
  theme_DOC()  +
  theme(legend.position = "none")
pDOCtimeDiff_BsM

pDOCTimeDiff_BsM_legTop = pDOCtimeDiff_BsM +
  theme(legend.position = "top")  #This is just to add legend to plot (custom)
pDOCTimeDiff_BsM_legRight = pDOCtimeDiff_BsM +
  theme(legend.position = "right")  #This is just to add legend to plot (custom)

# ggsave(filename = ("plotsTables/plots/doc_bsmMap_20260207.svg"),device = "svg", plot = pDOCtimeDiff_BsM,
#        width = 8, height = 6)
# ggsave(filename = ("plotsTables/plots/doc_bsmMap_legTop_202600207.svg"),device = "svg", plot = pDOCTimeDiff_BsM_legTop,
#        width = 8, height = 6)
# ggsave(filename = ("plotsTables/plots/doc_bsmMap_legRight_20260207.svg"),device = "svg", plot = pDOCTimeDiff_BsM_legRight,
#        width = 8, height = 6)

## Map of Secchi, only BsM
pSecchitimeDiff_BsM <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docBsMTimeDat, aes(y=lat,x=long,fill=diffSecchiperYear), 
             shape = 21, stroke = 0.5, colour = "black", alpha=0.8, size=4) +
  scale_fill_gradient2(
    high = "navyblue",            
    mid = "white",         
    low = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,57)) +
  labs(y="Latitude",x="Longitude",fill="Δ Secchi depth (m/year)") +
  theme_DOC() +
  theme(legend.position = "none")
pSecchitimeDiff_BsM

# ggsave(filename = ("plotsTables/plots/secchi_bsmMap_20260207.svg"),device = "svg", plot = pSecchitimeDiff_BsM,
#        width = 8, height = 6)


pSecchiTimeDiff_BsM_legTop = pSecchitimeDiff_BsM +
  theme(legend.position = "top")  #This is just to add legend to plot (custom)
# ggsave(filename = ("plotsTables/plots/secchi_bsmMap_legTop_20260207.svg"),device = "svg", plot = pSecchiTimeDiff_BsM_legTop,
#        width = 8, height = 6)

waterClarBsM_panelPlot <- pDOCTime_bsm+secchiBsmYearPlot+pDOCtimeDiff_BsM+pSecchitimeDiff_BsM + plot_layout(heights = c(1.5,2))

# Save as svg
# ggsave(filename = ("plotsTables/plots/waterClarBsM_panelPlot_20260207.svg"),device = "svg", plot = waterClarBsM_panelPlot,
#        width = 12, height = 12)

## Change in DOC and Secchi with environmental variables; ONLY BsM ####
# Model it all: full, biologically informed IC approach
docSpatial_catTrophicStatus <- glmmTMB(DOC ~ scaled_yearSample+
                                         scaled_maxDepth+
                                         lakeTrophicStatus+
                                         lat +
                                         scaled_yearSample:lakeTrophicStatus+
                                         scaled_yearSample:lat+
                                         (1|commonID), 
                                       data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_SI <- glmmTMB(DOC ~ scaled_yearSample+
                                            scaled_maxDepth+
                                            lakeTrophicStatus+
                                            lat +
                                            scaled_yearSample:lakeTrophicStatus+
                                            scaled_yearSample:lat+
                                            (scaled_yearSample|commonID), 
                                          data=mainDF_onlyBsM, family=Gamma(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_logN_I <- glmmTMB(DOC ~ scaled_yearSample+
                                                scaled_maxDepth+
                                                lakeTrophicStatus+
                                                lat +
                                                scaled_yearSample:lakeTrophicStatus+
                                                scaled_yearSample:lat+
                                                (1|commonID), 
                                              data=mainDF_onlyBsM, family=lognormal(link="log"), na.action = "na.fail")
docSpatial_catTrophicStatus_SI_logN_SI <- glmmTMB(DOC ~ scaled_yearSample+
                                                    scaled_maxDepth+
                                                    lakeTrophicStatus+
                                                    lat +
                                                    scaled_yearSample:lakeTrophicStatus+
                                                    scaled_yearSample:lat+
                                                    (scaled_yearSample|commonID), 
                                                  data=mainDF_onlyBsM, family=lognormal(link="log"), na.action = "na.fail")
AIC(docSpatial_catTrophicStatus,docSpatial_catTrophicStatus_SI,
    docSpatial_catTrophicStatus_logN_I,docSpatial_catTrophicStatus_SI_logN_SI)

#dredgeresults = dredge(docSpatial_catTrophicStatus)

bestmodel4=glmmTMB(DOC ~ scaled_yearSample+
                     scaled_maxDepth+
                     lakeTrophicStatus+
                     lat +
                     scaled_yearSample:lakeTrophicStatus+
                     (1|commonID), 
                   data=mainDF_onlyBsM, family=Gamma(link="log"))

# Investigate the results
summary(bestmodel4)
car::Anova(bestmodel4)
DHARMa::simulateResiduals(bestmodel4, plot=T)
hist(residuals(bestmodel4), breaks = 30, main = "Histogram of Residuals")
r.squaredGLMM(bestmodel4) 

# Predict **only use this for maxDepth**; can't parse out mean lakeTrophicStatus
mainDF_onlyBsM$predictedDOC_wEnvVars = predict(bestmodel4, type = "response", re.form = NA)
mainDF_onlyBsM$predictedDOC_wRandomEnvVars = predict(bestmodel4, type = "response",re.form = NULL)

pMaxDepthBsM <- ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=DOC, x=maxDepth, fill=DOC), alpha=0.5, size=2.5, shape=21)+
  geom_smooth(aes(x=maxDepth, y = DOC),
              colour = "blue",
              method="glm",
              formula = y~x,
              se = TRUE,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "DOC (mg/L)")+
  ylim(c(0,21)) +
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  theme_DOC()+
  theme(legend.position="none")
pMaxDepthBsM

# Plot partial regressions of interaction between year and trophic status
newdata = mainDF_onlyBsM
newdata$lat2 = newdata$lat
newdata$lat = mean(newdata$lat)
newdata$scaled_maxDepth = mean(newdata$scaled_maxDepth)

newdata$predictedDOC_wEnvVars = predict(bestmodel4,
                                        newdata = newdata,
                                        type = "response", re.form = NA)

# Need 95%CIs to plot; function to predict from bootstrapped models
bootFun <- function(model) {
  predict(model, newdata = newdata, type = "response", re.form = NA)
}

# Perform bootstrapping (adjust nsim for precision vs. speed)
bootResults <- bootMer(bestmodel4, FUN = bootFun, nsim = 1000, re.form = NA)

# Compute confidence intervals
newdata$predictedDOC_wEnvVars <- predict(bestmodel4, newdata = newdata, type = "response", re.form = NA)
newdata$lwr <- apply(bootResults$t, 2, quantile, probs = 0.025)  # 2.5% CI
newdata$upr <- apply(bootResults$t, 2, quantile, probs = 0.975)  # 97.5% CI

pSpatialDOCBsM <- ggplot(newdata)+
  geom_jitter(aes(y=DOC, x=yearSample, colour=lat2), alpha=0.3, size=2, width=0.22)+
  geom_ribbon(aes(x = yearSample, ymin = lwr, ymax = upr), fill = "#5E8FC6", alpha = 0.4) +  # Confidence interval
  geom_line(aes(x=yearSample, y = predictedDOC_wEnvVars), 
            colour = "blue", linewidth = 1) + 
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  ylim(c(0,21)) +
  labs(colour="Latitude", y="DOC (mg/L)", x="Year" )+
  theme_DOC()
pSpatialDOCBsM


# Summarize modelled differences in DOC over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData <- ggplot_build(pSpatialDOCBsM)
# Extract the data related to geom_smooth
lineData <- plotData$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumDOC_timeBsM <- lineData %>%
  group_by(PANEL) %>%  # PANEL corresponds to lakeTrophicStatus facet (1, Oligo; 2, Meso; 3, Eutro)
  summarise(minDOC = min(y),
            maxDOC = max(y),
            diffDOC = maxDOC - minDOC,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffDOC/yearDiff)*10)


## Now look at Secchi
# Use model formulation of final DOC model, sub in secchi depths and the whole dataset
secchiSpatial_catTrophicStatus <- glmmTMB(secchiDepth ~ scaled_yearSample+
                                            scaled_maxDepth+
                                            lakeTrophicStatus+
                                            lat +
                                            scaled_yearSample:lakeTrophicStatus+
                                            (1|commonID), 
                                          data=mainDF_onlyBsM, family=Gamma(link="log"))

#Investigate the results
summary(secchiSpatial_catTrophicStatus)
DHARMa::simulateResiduals(secchiSpatial_catTrophicStatus, plot=T)
hist(residuals(secchiSpatial_catTrophicStatus), breaks = 30, main = "Histogram of Residuals")
r.squaredGLMM(secchiSpatial_catTrophicStatus)  


# Plot some stuff
mainDF_onlyBsM$predictedSecchi_wEnvVars = predict(secchiSpatial_catTrophicStatus, type = "response", re.form = NA)
mainDF_onlyBsM$predictedSecchi_wRandomEnvVars = predict(secchiSpatial_catTrophicStatus, type = "response",re.form = NULL)

pMaxDepth_sec_BsM <- ggplot(mainDF_onlyBsM)+
  geom_point(aes(y=secchiDepth, x=maxDepth, fill = secchiDepth), alpha=0.5, size=2.5, shape = 21)+
  geom_smooth(aes(x=maxDepth, y = secchiDepth), 
              colour = "blue",
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  #ylim(c(0,21)) +
  theme(legend.position = "none")
pMaxDepth_sec_BsM

## Plot only partial regression containing trophic status * year interaction
newdata_sec = mainDF_onlyBsM
newdata_sec$lat2 = newdata_sec$lat
newdata_sec$lat = mean(newdata_sec$lat)
newdata_sec$scaled_maxDepth = mean(newdata_sec$scaled_maxDepth)

newdata_sec$predictedSecchi_wEnvVars = predict(secchiSpatial_catTrophicStatus,
                                               newdata = newdata_sec,
                                               type = "response", re.form = NA)

# Need 95%CIs to plot; function to predict from bootstrapped models
bootFun_sec <- function(model) {
  predict(model, newdata = newdata_sec, type = "response", re.form = NA)
}

# Perform bootstrapping (adjust nsim for precision vs. speed)
bootResults_sec <- bootMer(secchiSpatial_catTrophicStatus, FUN = bootFun_sec, nsim = 1000, re.form = NA)

# Compute confidence intervals
newdata_sec$predictedSecchi_wEnvVars <- predict(secchiSpatial_catTrophicStatus, newdata = newdata_sec, type = "response", re.form = NA)
newdata_sec$lwr <- apply(bootResults_sec$t, 2, quantile, probs = 0.025)  # 2.5% CI
newdata_sec$upr <- apply(bootResults_sec$t, 2, quantile, probs = 0.975)  # 97.5% CI


pSpatialSecBsM <- ggplot(newdata_sec)+
  geom_jitter(data = mainDF_onlyBsM, aes(y=secchiDepth, x=yearSample, colour=lat), alpha=0.3, size=2, width=0.22)+
  geom_line(aes(x=yearSample, y=predictedSecchi_wEnvVars),linewidth=1, colour="blue") +
  geom_ribbon(aes(x = yearSample, ymin = lwr, ymax = upr), fill = "#5E8FC6", alpha = 0.2) +  # Confidence interval
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)") +  
  labs(colour="Latitude")+
  # ylim(c(0,21)) +
  theme_DOC() 
pSpatialSecBsM

# Summarize modelled differences in Secchi over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_sec <- ggplot_build(pSpatialSecBsM)
# Extract the data related to geom_smooth
lineData_sec <- plotData_sec$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumSec_timeBsM <- lineData_sec %>%
  group_by(PANEL) %>%  # PANEL corresponds to lakeTrophicStatus facet (1, Oligo; 2, Meso; 3, Eutro)
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)

spatioTempBsM_panelPlot <- pSpatialDOCBsM + pMaxDepthBsM + plot_layout(ncol = 1)
spatioTempSecchiBsM_panelPlot <- pSpatialSecBsM + pMaxDepth_sec_BsM + plot_layout(ncol = 1)

# ggsave(filename = ("plotsTables/plots/spatioTempBsM_panelPlot_20260121.svg"),device = "svg", plot = spatioTempBsM_panelPlot,
#         width = 8, height = 6,bg="white")
# 
# ggsave(filename = ("plotsTables/plots/spatioTempSecchiBsM_panelPlot_20260121.svg"),device = "svg", plot = spatioTempSecchiBsM_panelPlot,
#        width = 8, height = 6,bg="white")


## Change in Secchi; BsM and AHI ####
## Now Secchi
M1_sec <- glmmTMB(secchiDepth ~ scaled_yearSample + (1|commonID),
                  data = mainDF,
                  family = Gamma(link = "log"))
M2_sec <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
                  data = mainDF,
                  family = Gamma(link = "log"))          
M3_sec = update(M1_sec, family = lognormal(link="log"))
M4_sec = update(M2_sec, family = lognormal(link="log"))
AIC(M1_sec,M2_sec,M3_sec,M4_sec)
# Gamma models are best. 

# Parameterize best model for DOC over time
bestmodel3_sec = M2_sec

# Diagnostics
DHARMa::simulateResiduals(bestmodel3_sec, plot=T)
hist(residuals(bestmodel3_sec), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumSec_allDat_noEnv <- summary(bestmodel3_sec)) #Secchi decreases by 0.37% per year; (3.61% decrease per decade); 0.12 m decrease per decade
exp(sumSec_allDat_noEnv$coefficients$cond[, "Estimate"]) 
r.squaredGLMM(bestmodel3_sec)          # pretty much all lake effect

mainDF$Secchi_fit_allData_noEnv = predict( bestmodel3_sec, type = "response", re.form = NA)
mainDF$Secchi_fitWRandom_allData_noEnv = predict( bestmodel3_sec, type = "response", re.form = NULL)  # To plot just predicted lines

#Plot it
pSecchi_allData_time <- ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = secchiDepth, fill=secchiDepth), alpha=0.7, size=2.5, shape = 21) +
  geom_smooth(aes(x=yearSample, y = secchiDepth),   # Use raw data, create a model - no partial effects as only one predictor (all marginal anyways)
              colour = "blue",
              method="glm",
              formula = y~x,
              se = TRUE,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214","#D2B48C","#F5DEB3")) +
  theme_DOC()+
  theme(legend.position="none")
pSecchi_allData_time

## Plot random slopes associated w. Secchi year model
#Add fitted values from the mixed model (conditional = include random effects)
mainDF$fit_cond_secchiYear <- predict(
  bestmodel3_sec,
  type = "response",   # get on original gamma scale
  re.form = NULL       # include random effects
)

#Remove instances where lakes were only sampled once to emphasize lakes w. repeated samples
plotDF_secYearAll <- mainDF %>%
  group_by(commonID) %>%
  filter(n() > 1) %>%
  ungroup()

#Plot lake-specific fitted lines
pSecchiBsMAndAHI_randomSlopes <- ggplot(plotDF_secYearAll, aes(x = yearSample, y = Secchi_fitWRandom_allData_noEnv)) +
  geom_point(aes(fill = Secchi_fitWRandom_allData_noEnv),
            alpha = 0.7, size = 2.5, shape = 21, width = 0.25) +
  geom_line(aes(y = fit_cond_secchiYear, colour = commonID), alpha=0.5) +   # <- one slope per lake
  scale_x_continuous("Year") +
  scale_y_continuous("Secchi depth (m)") +
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  scale_colour_viridis_d() +
  theme_DOC() +
  theme(legend.position = "none")
pSecchiBsMAndAHI_randomSlopes

# ggsave(filename = ("plotsTables/plots/secchiTime_bsmAndAHI_randomSlopes_20260121.svg"),device = "svg", plot = pSecchiBsMAndAHI_randomSlopes,
#        width = 8, height = 6)



# Summarize modelled differences in DOC over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_secAllYear <- ggplot_build(pSecchi_allData_time)
# Extract the data related to geom_smooth
lineData_secAllYear <- plotData_secAllYear$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumsecNoEnv_allYear <- lineData_secAllYear %>%
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)

## Change in Secchi with AHI and BsM
secTimeDat <- mainDF %>% 
  group_by(commonID) %>% 
  summarise(yearDiffSec = max(yearSample) - min(yearSample),
            initSecchi = Secchi_fitWRandom_allData_noEnv[which.min(yearSample)],
            finalSecchi = Secchi_fitWRandom_allData_noEnv[which.max(yearSample)],
            diffSecchi = finalSecchi-initSecchi,
            diffSecchiPerYear = diffSecchi/yearDiffSec,
            diffSecchiPerYear_cm = diffSecchiPerYear*100,
            sampleSize = n()) %>% 
  left_join(select(statDat_wide,commonID,lat,long)) %>% 
  distinct() %>% 
  filter(!is.na(lat)) %>% 
  mutate(dataset = "all")
secTimeDat


## Map of Secchi, including AHI
pSecchitimeDiff <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=secTimeDat, aes(y=lat,x=long,fill=diffSecchiPerYear_cm,), 
             shape = 21, stroke = 0.5, colour = "black",alpha=0.8, size=4) +
  scale_fill_gradient2(
    high = "navyblue",          
    mid = "white",         
    low = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,55)) +
  labs(y="Latitude",x="Longitude",fill="Δ Secchi depth (cm/year)") +
  theme_DOC() +
  theme(legend.position = "none")
pSecchitimeDiff

waterClar_panelPlot <- pSecchi_allData_time+pSecchitimeDiff + plot_layout(nrow=1)
# ggsave(filename = ("plotsTables/plots/waterClar_panelPlot_updated_20260207.svg"),device = "svg", plot = waterClar_panelPlot,
#        width = 12, height = 12)


waterClar_panelPlot_legTop = pSecchitimeDiff +
  theme(legend.position = "top")  #This is just to add legend to plot (custom)
# ggsave(filename = ("plotsTables/plots/secchi_allLakes_legTop_20260207.svg"),device = "svg", plot = waterClar_panelPlot_legTop,
#        width = 8, height = 6)
waterClar_panelPlot_legRight = pSecchitimeDiff +
  theme(legend.position = "right")  #This is just to add legend to plot (custom)
# ggsave(filename = ("plotsTables/plots/secchi_allLakes_legRight_20260207.svg"),device = "svg", plot = waterClar_panelPlot_legRight,
#        width = 8, height = 6)

## Using ALL data, see if BsM enviro results hold w. Secchi ####
# Expect similar direction and magnitude of effects
mainDF_secMod <- mainDF %>%
  group_by(commonID) %>%
  mutate(lakeTrophicStatus = case_when(
    TDP < 5 ~ "oligotrophic",
    TDP >= 5 & TDP < 10 ~ "mesotrophic",
    TDP >= 10 ~ "eutrophic")) %>%
  fill(lakeTrophicStatus, .direction = "downup")  # Fill NAs for ARU observations within each commonID

mainDF_secMod$lakeTrophicStatus <- factor(mainDF_secMod$lakeTrophicStatus, 
                                          levels = c("oligotrophic","mesotrophic","eutrophic"))  # just changing order for plotting

# Use model formulation of final DOC model, sub in secchi depths and the whole dataset
secchiSpatial_catTrophicStatus_allData <- glmmTMB(secchiDepth ~ scaled_yearSample+
                                                    scaled_maxDepth+
                                                    lakeTrophicStatus+
                                                    lat +
                                                    scaled_yearSample:lakeTrophicStatus+
                                                    (1|commonID), 
                                                  data=mainDF_secMod, family=Gamma(link="log"))


summary(secchiSpatial_catTrophicStatus_allData)
anova(secchiSpatial_catTrophicStatus_allData)

#plot some stuff

mainDF_secMod$predictedSecchi = predict(secchiSpatial_catTrophicStatus_allData, type = "response", re.form = NA)
mainDF_secMod$predictedSecchi_wRandom = predict(secchiSpatial_catTrophicStatus_allData, type = "response",re.form = NULL)

# Plot partial regressions of interaction between year and trophic status
newdata_allSec = mainDF_secMod
newdata_allSec$lat2 = newdata_allSec$lat
newdata_allSec$lat = mean(newdata_allSec$lat)
newdata_allSec$scaled_maxDepth = mean(newdata_allSec$scaled_maxDepth)

newdata_allSec$predictedSecchi = predict(secchiSpatial_catTrophicStatus_allData,
                                         newdata = newdata_allSec,
                                         type = "response", re.form = NA)

# Need 95%CIs to plot; function to predict from bootstrapped models; same as above w. only BsM data
bootFun_allSec <- function(model) {
  predict(model, newdata = newdata_allSec, type = "response", re.form = NA)
}

# Perform bootstrapping (adjust nsim for precision vs. speed)
bootResults_allSec <- bootMer(secchiSpatial_catTrophicStatus_allData, FUN = bootFun_allSec, nsim = 1000, re.form = NA)

# Compute confidence intervals
newdata_allSec$lwr <- apply(bootResults_allSec$t, 2, quantile, probs = 0.025)  # 2.5% CI
newdata_allSec$upr <- apply(bootResults_allSec$t, 2, quantile, probs = 0.975)  # 97.5% CI


pSpatialSecAll <- ggplot(data=newdata_allSec)+
  geom_point(aes(y=secchiDepth, x=yearSample, colour=lat2), alpha=0.3, size=2)+
  geom_line(aes(x=yearSample, y=predictedSecchi),linewidth=1, colour = "blue") +
  geom_ribbon(aes(x = yearSample, ymin = lwr, ymax = upr), fill = "#5E8FC6", alpha = 0.2) +  # Confidence interval
  scale_colour_viridis_c(option="inferno", direction = -1) +
  facet_wrap(~lakeTrophicStatus, labeller = labeller(lakeTrophicStatus = function(x) str_to_title(x)))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "Secchi depth (m)") +  
  labs(colour="Latitude")+
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pSpatialSecAll

pMaxDepth_secAll <- ggplot(mainDF_secMod)+
  geom_point(aes(y=secchiDepth, x=maxDepth, fill = secchiDepth), alpha=0.7, size=2.5, shape = 21)+
  geom_smooth(aes(x=maxDepth, y = secchiDepth), 
              colour="blue",
              method="glm",
              formula = y~x,
              se=TRUE, level = 0.95,  #Plot 95%CI
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Maximum depth (m)")+
  scale_y_continuous(name = "Secchi depth (m)")+
  scale_fill_gradientn(colours = c("#3D2B1F","#704214", "#D2B48C", "#F5DEB3")) +
  theme_DOC()+
  theme(legend.position = "none")
pMaxDepth_secAll

spatioTemp_panelPlotAll <- pSpatialSecAll + pMaxDepth_secAll + plot_layout(ncol = 1)

# ggsave(filename = ("plotsTables/plots/spatioTemp_panelPlotAll_20260121.svg"),device = "svg", plot = spatioTemp_panelPlotAll,
#        width = 8, height = 6,bg="white")



# Summarize modelled differences in Secchi over time period for each trophic group
# Extract data from plot above - want the min/max of the fitted lines
plotData_secAll <- ggplot_build(pSpatialSecAll)
# Extract the data related to geom_smooth
lineData_secAll <- plotData_secAll$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumSec_timeAll <- lineData_secAll %>%
  group_by(PANEL) %>%  # PANEL corresponds to lakeTrophicStatus facet (1, Oligo; 2, Meso; 3, Eutro)
  summarise(minSec = min(y),
            maxSec = max(y),
            diffSec = maxSec - minSec,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffSec/yearDiff)*10)




## Testing observed secchi slope; different than chance? ####
#This is in response to reviewer concern that a longer time series (that includes BsM and AHI)
# smooths over decadal-scal variability that may be detected as directional change in shorter
# time periods (like that of BsM only)

# Try adding a "period" term into the model, and test for an interaction
# If significant interaction, we then know something about period differences, if no interaction
# then it's a general trend over time

# Step 1: Fit model to the actual present period
present_mod <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample | commonID),
                       data = mainDF_onlyBsM,
                       family = Gamma(link = "log"))

observed_slope <- summary(present_mod)$coefficients$cond["scaled_yearSample", "Estimate"]

summary(present_mod)

# Step 2: Use all data, include an interaction term
interactionSec_mod <- glmmTMB(secchiDepth ~ scaled_yearSample*samplingProgram + (scaled_yearSample | commonID),
                              data = mainDF,
                              family = Gamma(link = "log"))

summary(interactionSec_mod) ## Significant interaction, time points are different
# Model summary indicate

# Predict this common slope out
min_scaled_year <- min(mainDF$scaled_yearSample)
max_scaled_year <- max(mainDF$scaled_yearSample)
predDat_commonSlope <- data.frame(scaled_yearSample = seq(min_scaled_year, max_scaled_year, length.out = 100))

predDat_commonSlope$predicted_secchiDepth <- predict(commonslope_Secmod,
                                                     newdata = predDat_commonSlope,
                                                     re.form = NA, # Ignore random effects for common slope
                                                     type = "response")



# Step 4; plot it all
pIntMod_secBsMandAHI <- ggplot(mainDF, aes(x = scaled_yearSample, y= secchiDepth, colour = samplingProgram)) +
  geom_point(aes(fill = samplingProgram), shape = 21, colour = "black", stroke = 0.5, alpha = 0.7) +
  geom_smooth(aes(x=scaled_yearSample, y = secchiDepth), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log'))) +
  scale_x_continuous(breaks = c(-40,-20,0, 20), labels = c(1960,1980,2000,2020)) +
  # Add the common slope line
  geom_line(data = predDat_commonSlope,
            aes(x = scaled_yearSample, y = predicted_secchiDepth),
            inherit.aes = FALSE, 
            linetype = "dashed", 
            color = "black",     
            size = 1) +          
  scale_fill_manual(values = c("#278B9AFF", "#D8AF39FF")) +
  scale_colour_manual(values = c("#278B9AFF", "#D8AF39FF")) +
  labs(x= "Year sample", y = "Secchi depth (m)") +
  guides(fill = "none", colour = "none") +
  theme(legend.position = "none") +
  theme_DOC()
pIntMod_secBsMandAHI

# ggsave(filename = ("plotsTables/plots/interactionModel_secchiBsMandAHI_20260201.png"),device = "png", plot = pIntMod_secBsMandAHI,
#        width = 8, height = 6,bg="white")

## SUPPLEMENT - PLOTTING RANDOM SLOPES FOR SECCHI MODELS #### 
# For models where random slopes + intercepts significant, visualizing lake-specific effects demonstrates 
# clear lakes have the largest capacity for change

## To show meaningful change, the mean and SD of all slopes is calculated for a population average
# (which is equivalent to the slope of the models fixed effect). Lines are then colour-coded to
# represent cases where decreases in water clarity are quite steep (lake-specific slope is less
# than the mean - 1SD), changes are not strong (slope is within 1 SD of the mean), and lakes where 
# clarity is actually increasing (lake-specific slope is > mean slope plus 1 SD). This is done
# for both the BsM data in isolation as well as the BsM + AHI dataset

#### BsM only model
#Add fitted values from the mixed model (conditional = include random effects)
mainDF_onlyBsM$fit_cond_secchiYear <- predict(
  secchiTimeModb,
  type = "response",   # get on original gamma scale
  re.form = NULL       # include random effects
)

# Calculate the "Population Slope" for each lake on the response scale
plotDF_secYearBsM <- mainDF_onlyBsM %>%
  group_by(commonID) %>%
  filter(n() > 1) %>%
  mutate(slope_response_scale = (last(fit_cond_secchiYear) - first(fit_cond_secchiYear)) / 
        (last(yearSample) - first(yearSample))) %>%
  ungroup()

# Calculate the average slope for all lakes on the same response scale
lake_stats <- plotDF_secYearBsM %>%
  select(commonID, slope_response_scale) %>%
  distinct() %>%
  summarise( mean_s = mean(slope_response_scale, na.rm = TRUE),
             sd_s = sd(slope_response_scale,na.rm = TRUE),
             se_s   = sd(slope_response_scale, na.rm = TRUE) / sqrt(n()))  #Also calculate SE..maybe plot

# Extract BsM values
mean_val <- lake_stats$mean_s
sd_val <- lake_stats$sd_s
se_val   <- lake_stats$se_s

# Create categories for plotting
plotDF_secYearBsM <- plotDF_secYearBsM %>%
  mutate(slope_cat = case_when(
      slope_response_scale < (mean_val - (sd_val)) ~ "Decrease (slope < mean - 1SD)",   # Steep decrease change
      slope_response_scale > (mean_val + (sd_val)) ~ "Increase (slope > mean + 1SD)",    # Steep negative change
      # slope_response_scale < (mean_val - (2*se_val)) ~ "Decrease (slope < mean - 2SE)",   # Steep decrease change
      # slope_response_scale > (mean_val + (2*se_val)) ~ "Increase (slope > mean + 2SE)",    # Steep negative change
      TRUE ~ "No change (within 1SD of mean)"))                                        # Flat / Average change

# Plot it
pSecchiBsM_randomSlopes <- ggplot(plotDF_secYearBsM, aes(x = yearSample, y = Secchi_fitWRandom)) +
  geom_point(alpha = 0.7, size = 2.5, shape = 21) +
    geom_line(aes(y = fit_cond_secchiYear, group = commonID, colour = slope_cat), 
              alpha = 0.5, lwd=1.1) +
  scale_x_continuous("Year") +
  scale_y_continuous("Secchi depth (m)") +
  #Define specific colors
  scale_colour_manual(values = c("Decrease (slope < mean - 1SD)" = "#D8AF39FF",
                                 "No change (within 1SD of mean)" = "#E75B64FF",
                                 "Increase (slope > mean + 1SD)"  = "#278B9AFF"),
                      name = "Slope Relative to Population Mean:") +
  # scale_colour_manual(values = c("Decrease (slope < mean - 2SE)" = "#D8AF39FF",
  #                                "No change" = "#E75B64FF",
  #                                "Increase (slope > mean + 2SE)"  = "#278B9AFF"),
  #                     name = "Slope Relative to Population Mean:") +
  theme_DOC() +
  facet_wrap(~slope_cat) +
  theme(legend.position = "none")
pSecchiBsM_randomSlopes

# ggsave(filename = ("plotsTables/plots/randomSlopes_secchiBsM_20260131.svg"),device = "svg", plot = pSecchiBsM_randomSlopes,
#        width = 12, height = 8,bg="white")

#### BsM and AHI data
#Add fitted values from the mixed model (conditional = include random effects)
mainDF$fit_cond_secchiYear <- predict(
  bestmodel3_sec,
  type = "response",   # get on original gamma scale
  re.form = NULL       # include random effects
)

# Calculate the "Population Slope" for each lake on the response scale
plotDF_secYearAll <- mainDF %>%
  group_by(commonID) %>%
  filter(n() > 1) %>% 
  mutate(slope_response_scale = (last(fit_cond_secchiYear) - first(fit_cond_secchiYear)) / 
           (last(yearSample) - first(yearSample))) %>%
  ungroup()

# Calculate the average slope for all lakes on the same response scale
lake_statsAll <- plotDF_secYearAll %>%
  select(commonID, slope_response_scale) %>%
  distinct() %>%
  summarise( mean_s = mean(slope_response_scale, na.rm = TRUE),
             sd_s = sd(slope_response_scale,na.rm = TRUE),
             se_s   = sd(slope_response_scale, na.rm = TRUE) / sqrt(n()))  #Also calculate SE..maybe plot

# Extract BsM values
mean_valAll <- lake_statsAll$mean_s
sd_valAll <- lake_statsAll$sd_s
se_valAll   <- lake_statsAll$se_s

# Create categories for plotting
plotDF_secYearAll<- plotDF_secYearAll %>%
  mutate(slope_cat = case_when(
    slope_response_scale < (mean_valAll - (sd_valAll)) ~ "Decrease (slope < mean - 1SD)",   # Steep decrease change
    slope_response_scale > (mean_valAll + (sd_valAll)) ~ "Increase (slope > mean + 1SD)",    # Steep negative change
    # slope_response_scale < (mean_valAll - (2*se_valAll)) ~ "Decrease (slope < mean - 2SE)",   # Steep decrease change
    # slope_response_scale > (mean_valAll + (2*se_valAll)) ~ "Increase (slope > mean + 2SE)",    # Steep negative change
    TRUE ~ "No change (within 1SD of mean)"))                                        # Flat / Average change

# Plot it
pSecchiAll_randomSlopes <- ggplot(plotDF_secYearAll, aes(x = yearSample, y = Secchi_fitWRandom_allData_noEnv)) +
  geom_point(alpha = 0.7, size = 2.5, shape = 21) +
  geom_line(aes(y = fit_cond_secchiYear, group = commonID, colour = slope_cat), 
            alpha = 0.5, lwd=1.1) +
  scale_x_continuous("Year") +
  scale_y_continuous("Secchi depth (m)") +
  #Define specific colors
  scale_colour_manual(values = c("Decrease (slope < mean - 1SD)" = "#D8AF39FF",
                                 "No change (within 1SD of mean)" = "#E75B64FF",
                                 "Increase (slope > mean + 1SD)"  = "#278B9AFF"),
                      name = "Slope Relative to Population Mean:") +
  # scale_colour_manual(values = c("Decrease (slope < mean - 2SE)" = "#D8AF39FF",
  #                                "No change" = "#E75B64FF",
  #                                "Increase (slope > mean + 2SE)"  = "#278B9AFF"),
  #                     name = "Slope Relative to Population Mean:") +
  theme_DOC() +
  facet_wrap(~slope_cat) +
  theme(legend.position = "none")
pSecchiAll_randomSlopes
# 
# ggsave(filename = ("plotsTables/plots/randomSlopes_secchiAll_20260131.svg"),device = "svg", plot = pSecchiAll_randomSlopes,
#        width = 12, height = 8,bg="white")

## DATA EXPLORATION; HISTOGRAMS ####

## Histogram of lake size and trophic statuses
histSize <- ggplot(ontario_ohnWaterbody_latlon) +
  #geom_histogram(aes(x=SYSTEM_CALCULATED_AREA))
  geom_histogram(aes(x=(SYSTEM_CALCULATED_AREA/10000))) +
  scale_x_log10(labels = scales::label_number()) +  # Removes scientific notation  theme_minimal()
  xlab("Surface area (ha)") +
  ylab("Frequency") +
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~geoRegion)
histSize

histSize <- ggplot(ontario_ohnWaterbody_latlon) +
  #geom_histogram(aes(x=SYSTEM_CALCULATED_AREA))
  geom_histogram(aes(x=(SYSTEM_CALCULATED_AREA/10000))) +
  scale_x_log10(labels = scales::label_number()) +  # Removes scientific notation  theme_minimal()
  xlab("Surface area (ha)") +
  ylab("Frequency") +
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~geoRegion)
histSize

mainDF_onlyBsM_wSA <- mainDF_onlyBsM %>% 
  left_join(select(statDat_wide,commonID,surfaceArea)) %>% 
  distinct()

histTrophicStatus_loc <- ggplot(mainDF_onlyBsM_wSA) +
  geom_histogram(aes(x=surfaceArea)) +
  facet_wrap(~lakeTrophicStatus) +
  scale_x_log10(labels = scales::label_number()) +  # Removes scientific notation  theme_minimal()
  theme_DOC() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
histTrophicStatus_loc

# pCor_TDPSec <- ggplot(mainDF_onlyBsM) +
#   geom_point(aes(x=TDP,y=secchiDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
#   geom_smooth(aes(x=(TDP),y=(secchiDepth)), method = "lm") +
#   scale_colour_viridis_c(option="inferno", direction = -1) +
#   scale_x_continuous(trans='log', labels = function(x) signif(x, 2)) +
#   scale_y_continuous(trans='log', labels = function(x) signif(x, 2)) +
#   ylab("Secchi (m)") +
#   xlab("Total dissolved phosphorus (ug/L)") +
#   # ylim(c(0,21))+
#   # xlim(c(0,21)) +
#   theme_DOC() 
# pCor_TDPSec

pCor_TDPSec <- ggplot(mainDF_onlyBsM) +
  geom_point(aes(x=yearSample,y=TDP), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_smooth(aes(x=(yearSample),y=(TDP)), method = "lm") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  #scale_x_continuous(trans='log', labels = function(x) signif(x, 2)) +
  #scale_y_continuous(trans='log', labels = function(x) signif(x, 2)) +
  ylab("Secchi (m)") +
  xlab("Total dissolved phosphorus (ug/L)") +
  # ylim(c(0,21))+
  # xlim(c(0,21)) +
  theme_DOC() 
pCor_TDPSec

## SUPPLEMENT - sensitivity analysis of models ####
# See if lakes with sparse data provide similar results to those with more

library(dotwhisker)
library(broom.mixed)

## Create new datasets of "sparse" and "robust" data; BsM data
bsmDat_sparse <- mainDF_onlyBsM %>% 
  group_by(commonID) %>% 
  filter(n() == 2) %>% 
  ungroup()
bsmDat_robust <- mainDF_onlyBsM %>% 
  group_by(commonID) %>% 
  filter(n() > 2) %>% 
  ungroup()

length(unique(bsmDat_sparse$commonID)) #330
length(unique(bsmDat_robust$commonID)) #113
length(unique(mainDF_onlyBsM$commonID)) #684



## DOC over time
docTimeBsM_originalModel <- glmmTMB(DOC ~ scaled_yearSample + (1|commonID),
                                    data = mainDF_onlyBsM,
                                    family = Gamma(link = "log"))
docTimeBsM_sparseModel <- update(docTimeBsM_originalModel,data=bsmDat_sparse)
docTimeBsM_robustModel <- update(docTimeBsM_originalModel,data=bsmDat_robust)

# Combine models into a named list
modelsDOC_bsm <- list(
  "Original (n = 684)" = docTimeBsM_originalModel,
  "Sparse (n = 330)" = docTimeBsM_sparseModel,
  "Robust (n = 113)" = docTimeBsM_robustModel
)

# Create a forest plot of the 'scaled_yearSample' coefficient
docBsM_sensAnal <- dwplot(modelsDOC_bsm, effects = "fixed") |>
  relabel_predictors(scaled_yearSample = "Year sample") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "a) Effect of Year on DOC",
       x = "Estimate (Log scale)",
       y = "") +
  scale_colour_discrete(name = "Model Subset") +
  theme_DOC()
docBsM_sensAnal

# ggsave(filename = ("plotsTables/plots/sensitivityAnalysis_DOCfromBsM_20260213.svg"),device = "svg", plot = docBsM_sensAnal,
#        width = 8, height = 6,bg="white")

## Secchi over time
secTimeBsM_originalModel <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
                                    data = mainDF_onlyBsM,
                                    family = Gamma(link = "log"))   
secTimeBsM_sparseModel <- update(secTimeBsM_originalModel,data=bsmDat_sparse)
secTimeBsM_robustModel <- update(secTimeBsM_originalModel,data=bsmDat_robust)

# Combine models into a named list
modelsSec_BsM <- list(
  "Original (n = 684)" = secTimeBsM_originalModel,
  "Sparse (n = 330)" = secTimeBsM_sparseModel,
  "Robust (n = 113)" = secTimeBsM_robustModel
)

# Create a forest plot of the 'scaled_yearSample' coefficient
secBsM_sensAnal <- dwplot(modelsSec_BsM, effects = "fixed") |>
  relabel_predictors(c(scaled_yearSample = "Year sample")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "b) Effect of Year on Secchi",
       x = "Estimate (Log Scale)",
       y = "") +
  scale_colour_discrete(name = "Model Subset") +
  theme_DOC()
secBsM_sensAnal

# ggsave(filename = ("plotsTables/plots/sensitivityAnalysis_secFromBsM_20260213.svg"),device = "svg", plot = secBsM_sensAnal,
#        width = 8, height = 6,bg="white")

## DOC with spatial variables
docTimeBsM_spatial_originalModel <- glmmTMB(DOC ~ scaled_yearSample+
                                    scaled_maxDepth+
                                    lakeTrophicStatus+
                                    lat +
                                    scaled_yearSample:lakeTrophicStatus+
                                    (1|commonID), 
                                    data=mainDF_onlyBsM, family=Gamma(link="log"))
docTimeBsM_spatial_sparseModel <- update(docTimeBsM_spatial_originalModel,data=bsmDat_sparse)
docTimeBsM_spatial_robustModel <- update(docTimeBsM_spatial_originalModel,data=bsmDat_robust)

# Combine models into a named list
modelsDOC_spatial_BsM <- list(
  "Original (n = 684)" = docTimeBsM_spatial_originalModel,
  "Sparse (n = 330)" = docTimeBsM_spatial_sparseModel,
  "Robust (n = 113)" = docTimeBsM_spatial_robustModel
)

# Create a forest plot of the 'scaled_yearSample' coefficient
docSpatialBsM_sensAnal <- dwplot(modelsDOC_spatial_BsM, effects = "fixed") |>
  relabel_predictors(scaled_yearSample = "Year sample",
                     scaled_maxDepth = "Max. depth",
                     lakeTrophicStatusmesotrophic = "Mesotrophic lakes",
                     lakeTrophicStatuseutrophic = "Eutrophic lakes",
                     lat = "Latitude",
                     `scaled_yearSample:lakeTrophicStatusmesotrophic` = "Year sample * Mesotrophic",
                     `scaled_yearSample:lakeTrophicStatuseutrophic` = "Year sample * Eutrophic") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "c) Effect of Year + Env. on DOC",
       x = "Estimate (Log Scale)",
       y = "") +
  scale_colour_discrete(name = "Model Subset") +
  theme_DOC()
docSpatialBsM_sensAnal

# ggsave(filename = ("plotsTables/plots/sensitivityAnalysis_DOCSpatialFromBsM_20260213.svg"),device = "svg", plot = docSpatialBsM_sensAnal,
#        width = 8, height = 6,bg="white")

## Secchi with spatial variables
secTimeBsM_spatial_originalModel <- glmmTMB(secchiDepth ~ scaled_yearSample+
                                              scaled_maxDepth+
                                              lakeTrophicStatus+
                                              lat +
                                              scaled_yearSample:lakeTrophicStatus+
                                              (1|commonID), 
                                            data=mainDF_onlyBsM, family=Gamma(link="log"))
secTimeBsM_spatial_sparseModel <- update(secTimeBsM_spatial_originalModel,data=bsmDat_sparse)
secTimeBsM_spatial_robustModel <- update(secTimeBsM_spatial_originalModel,data=bsmDat_robust)

# Combine models into a named list
modelsSec_spatial_BsM <- list(
  "Original (n = 684)" = secTimeBsM_spatial_originalModel,
  "Sparse (n = 330)" = secTimeBsM_spatial_sparseModel,
  "Robust (n = 113)" = secTimeBsM_spatial_robustModel
)

# Create a forest plot of the 'scaled_yearSample' coefficient
secSpatialBsM_sensAnal <- dwplot(modelsSec_spatial_BsM, effects = "fixed") |>
  relabel_predictors(scaled_yearSample = "Year sample",
                     scaled_maxDepth = "Max. depth",
                     lakeTrophicStatusmesotrophic = "Mesotrophic lakes",
                     lakeTrophicStatuseutrophic = "Eutrophic lakes",
                     lat = "Latitude",
                     `scaled_yearSample:lakeTrophicStatusmesotrophic` = "Year sample * Mesotrophic",
                     `scaled_yearSample:lakeTrophicStatuseutrophic` = "Year sample * Eutrophic") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "d) Effect of Year + Env. on Secchi",
       x = "Estimate (Log Scale)",
       y = "") +
  scale_colour_discrete(name = "Model Subset") +
  theme_DOC()
secSpatialBsM_sensAnal

# ggsave(filename = ("plotsTables/plots/sensitivityAnalysis_secSpatialFromBsM_20260213.svg"),device = "svg", plot = secSpatialBsM_sensAnal,
#        width = 8, height = 6,bg="white")


## Now look at Secchi ~ yearSample for all BsM + AHI
allDat_sparse <- mainDF_secMod %>% 
  group_by(commonID) %>% 
  filter(n() == 2) %>% 
  ungroup()
allDat_robust <- mainDF_secMod %>% 
  group_by(commonID) %>% 
  filter(n() > 2) %>% 
  ungroup()

length(unique(mainDF_secMod$commonID)) #684
length(unique(allDat_sparse$commonID)) #248
length(unique(allDat_robust$commonID)) #436


secTimeAllDat_originalModel <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample|commonID),
                  data = mainDF,
                  family = Gamma(link = "log"))
secTimeAllDat_sparseModel <- update(secTimeAllDat_originalModel,data=allDat_sparse)
secTimeAllDat_robustModel <- update(secTimeAllDat_originalModel,data=allDat_robust)

# Combine models into a named list
modelsSec_allDat <- list(
  "Original (n = 684)" = secTimeAllDat_originalModel,
  "Sparse (n = 248)" = secTimeAllDat_sparseModel,
  "Robust (n = 436)" = secTimeAllDat_robustModel
)

# Create a forest plot of the 'scaled_yearSample' coefficient
secAll_sensAnal <- dwplot(modelsSec_allDat, effects = "fixed") |>
  relabel_predictors(scaled_yearSample = "Year sample") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "a) Effect of Year on Secchi",
       x = "Estimate (Log Scale)",
       y = "") +
  scale_colour_discrete(name = "Model Subset") +
  theme_DOC()
secAll_sensAnal

# ggsave(filename = ("plotsTables/plots/sensitivityAnalysis_secFromAllData_20260213.svg"),device = "svg", plot = secAll_sensAnal,
#        width = 8, height = 6,bg="white")

## Now look at Secchi ~ spatial variables for all BsM + AHI
secSpatialAllData_originalModel <- glmmTMB(secchiDepth ~ scaled_yearSample+
                                                    scaled_maxDepth+
                                                    lakeTrophicStatus+
                                                    lat +
                                                    scaled_yearSample:lakeTrophicStatus+
                                                    (1|commonID), 
                                                  data=mainDF_secMod, family=Gamma(link="log"))
secSpatialAllDat_sparseModel <- update(secSpatialAllData_originalModel,data=allDat_sparse)
secSpatialAllDat_robustModel <- update(secSpatialAllData_originalModel,data=allDat_robust)

# Combine models into a named list
modelsSecSpatial_allDat <- list(
  "Original (n = 684)" = secSpatialAllData_originalModel,
  "Sparse (n = 248)" = secSpatialAllDat_sparseModel,
  "Robust (n = 436)" = secSpatialAllDat_robustModel
)

# Create a forest plot of the 'scaled_yearSample' coefficient
secAllSpatial_sensAnal <- dwplot(modelsSecSpatial_allDat, effects = "fixed") |>
  relabel_predictors(scaled_yearSample = "Year sample",
                     scaled_maxDepth = "Max. depth",
                     lakeTrophicStatusmesotrophic = "Mesotrophic lakes",
                     lakeTrophicStatuseutrophic = "Eutrophic lakes",
                     lat = "Latitude",
                     `scaled_yearSample:lakeTrophicStatusmesotrophic` = "Year sample * Mesotrophic",
                     `scaled_yearSample:lakeTrophicStatuseutrophic` = "Year sample * Eutrophic") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "b) Effect of Year + Env. on Secchi",
       x = "Estimate (Log Scale)",
       y = "Fixed Effect") +
  scale_colour_discrete(name = "Model Subset") +
  theme_DOC()
secAllSpatial_sensAnal

# ggsave(filename = ("plotsTables/plots/sensitivityAnalysis_secSpatialFromAllData_20260213.svg"),device = "svg", plot = secAllSpatial_sensAnal,
#        width = 8, height = 6,bg="white")




# NOT USED ANYMORE - Estimate DOC for BsM/AHI secchi interpretation ####
# Can we predict DOC w. lake variables? Secchi, depth, lat? ####
mainDF_BsM <- filter(mainDF, samplingProgram %in% "BsM")

# Fit the initial model - create full models with log-normal and gamma, look at AIC
clarDOCmod_gamma <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID), 
                            data = mainDF_BsM, family = Gamma(link = "log"), na.action = "na.fail")   # added lat in prediction because predictions w. only secchi and maxDepth under-predicted high DOC lakes. Only included as additive effect as there's a S/N gradient in DOC across the province
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_gamma)) #full model best
# clarDOCmod_LMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
#                             data = mainDF_BsM,  na.action = "na.fail")
# (clarDOCmod_LMM_sel <- dredge(clarDOCmod_LMM))
# clarDOCmod_logLMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
#                       data = mainDF_BsM, family = lognormal(link="log"),  na.action = "na.fail")
# (clarDOCmod_logLMM_sel <- dredge(clarDOCmod_logLMM))
AIC(clarDOCmod_gamma,clarDOCmod_LMM,clarDOCmod_logLMM)
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_gamma)) #full model best

# AICc lowest with gamma. Do that.

bestmodel2 <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth + lat + (1|commonID),
                      data = mainDF_BsM, family = Gamma(link = "log"))

# Diagnostics
DHARMa::simulateResiduals(bestmodel2, plot=T)
hist(residuals(bestmodel2), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Summary
summary(bestmodel2)
r2_secchiMod <- r.squaredGLMM(bestmodel2)        # R2 = 0.49 fixed effects, R2 = 0.92 w. random effects 
car::Anova(bestmodel2,type = "III")

# Partial variance explained 
# Fit reduced models
model_no_secchi <- update(bestmodel2, . ~ . - scaled_secchi)
model_no_depth <- update(bestmodel2, . ~ . - scaled_maxDepth)
model_no_lat <- update(bestmodel2, . ~ . - lat)
model_no_int <- update(bestmodel2, . ~ . - scaled_secchi:scaled_maxDepth)

# Calculate R² for reduced models
r2_no_secchi <- r.squaredGLMM(model_no_secchi)
r2_no_depth <- r.squaredGLMM(model_no_depth)
r2_no_lat <- r.squaredGLMM(model_no_lat)
r2_no_int <- r.squaredGLMM(model_no_int)

# Calculate partial R²
partial_R2_secchi <- r2_secchiMod[1] - r2_no_secchi[1]
partial_R2_depth <- r2_secchiMod[1] - r2_no_depth[1]
partial_R2_lat <- r2_secchiMod[1] - r2_no_lat[1]
partial_R2_int <- r2_secchiMod[1] - r2_no_int[1]


# Print results
cat("Partial R² for scaled_secchi: ", partial_R2_secchi, "\n")
cat("Partial R² for scaled_maxDepth: ", partial_R2_depth, "\n")
cat("Partial R² for lat: ", partial_R2_lat, "\n")
cat("Partial R² for int: ", partial_R2_int, "\n")

# Add predicted values to df
predictedVals <- predict(bestmodel2, type = "response", re.form=NA) #for discussion (should random effects be used or not)
#predictedVals <- predict(bestmodel2, type = "response", re.form=NULL) #for discussion (should random effects be used or not)
mainDF_BsM$predictedDOC <- predictedVals
cor(mainDF_BsM$DOC, mainDF_BsM$predictedDOC)  # obs-fitted have 0.75 correlation - good
cor.test(mainDF_BsM$DOC, mainDF_BsM$predictedDOC)
r.squaredGLMM(bestmodel2)
## Calculate RMSE
# Compute residuals
residualsDOCmod <- mainDF_BsM$DOC - predictedVals
# Compute RMSE
(rmseDOCMod <- sqrt(mean(residualsDOCmod^2)))

## Relationship between observed and predicted DOC for BSM data
pObsFit <- ggplot(mainDF_BsM) +
  geom_point(aes(x=DOC,y=predictedDOC, colour = lat, size=maxDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  ylab("Model Predicted DOC (mg/L)") +
  xlab("Measured DOC (mg/L)") +
  labs(colour="Latitude", size="Max. depth (m)")+
  ylim(c(0,21))+
  xlim(c(0,21)) +
  theme_minimal() +
  theme_DOC()
pObsFit

# ggsave(filename = ("plotsTables/plots/obsFittedDOC_20241118.png"),device = "png", plot = pObsFit,
#        width = 7, height = 6, units = "in", dpi = 450)

# If we do the same with secchi, do predictors have the same sign ####
# i.e., do we see water clarity changing in a complimentary way to DOC?
bestmodel2b <- glmmTMB(secchiDepth ~ scaled_maxDepth + lat + (1|commonID),
                       data = mainDF_BsM, family = Gamma(link = "log"))

# Diagnostics
DHARMa::simulateResiduals(bestmodel2b, plot=T)
hist(residuals(bestmodel2b), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Summary
summary(bestmodel2b)             # Beauty, response moving in same direction as DOC, and coefficients (i.e., slopes) are similar
r2_secchiMod <- r.squaredGLMM(bestmodel2b)        # R2 = 0.30 fixed effects, R2 = 0.77 w. random effects 
car::Anova(bestmodel2,type = "III")

# Partial variance explained 
# Fit reduced models
model_no_depth <- update(bestmodel2b, . ~ . - scaled_maxDepth)
model_no_lat <- update(bestmodel2b, . ~ . - lat)

# Calculate R² for reduced models
r2_no_depth <- r.squaredGLMM(model_no_depth)
r2_no_lat <- r.squaredGLMM(model_no_lat)

# Calculate partial R²
partial_R2_depth <- r2_secchiMod[1] - r2_no_depth[1]
partial_R2_lat <- r2_secchiMod[1] - r2_no_lat[1]


# Print results
cat("Partial R² for scaled_maxDepth: ", partial_R2_depth, "\n")
cat("Partial R² for lat: ", partial_R2_lat, "\n")

# Add predicted values to df
predictedVals <- predict(bestmodel2b, type = "response", re.form=NA) #for discussion (should random effects be used or not)
mainDF_BsM$predictedSecchi <- predictedVals
cor(mainDF_BsM$secchiDepth, mainDF_BsM$predictedSecchi)  # obs-fitted have 0.48 correlation - good
cor.test(mainDF_BsM$secchiDepth, mainDF_BsM$predictedSecchi)

# Can we predict DOC w. lake variables if we know trophic status? Secchi, depth, lat, trophic status? ####
mainDF_BsM_wTrophic <-  mainDF %>% 
  filter(samplingProgram == "BsM") %>% 
  mutate(lakeTrophicStatus = case_when(TDP < 5 ~ "oligotrophic",
                                       TDP >= 5 & TDP < 10 ~ "mesotrophic",
                                       TDP >= 10 ~ "eutrophic"),
         scaledTDP = as.vector(scale(TDP, scale = FALSE)))

# Fit the initial model - create full models with log-normal and gamma, look at AIC
clarDOCmod_wTrophic_gamma <- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth * lakeTrophicStatus + lat + elevation + (1|commonID), 
                                     data = mainDF_BsM_wTrophic, family = Gamma(link = "log"), na.action = "na.fail") 
#data = filter(mainDF_BsM_wTrophic, lakeTrophicStatus != "eutrophic"), family = Gamma(link = "log"), na.action = "na.fail")   # added lat in prediction because predictions w. only secchi and maxDepth under-predicted high DOC lakes. Only included as additive effect as there's a S/N gradient in DOC across the province
# clarDOCmod_wTrophic_LMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth * lakeTrophicStatus + lat + elevation + (1|commonID),
#                             data = mainDF_BsM_wTrophic,  na.action = "na.fail")
# (clarDOCmod_LMM_sel <- dredge(clarDOCmod_LMM))
# clarDOCmod_wTrophic_logLMM<- glmmTMB(DOC ~ scaled_secchi * scaled_maxDepth * lakeTrophicStatus + lat + elevation + (1|commonID),
#                       data = mainDF_BsM_wTrophic, family = lognormal(link="log"),  na.action = "na.fail")
# (clarDOCmod_logLMM_sel <- dredge(clarDOCmod_logLMM))
AIC(clarDOCmod_wTrophic_gamma,clarDOCmod_wTrophic_LMM,clarDOCmod_wTrophic_logLMM)
(clarDOCmod_gamma_sel <- dredge(clarDOCmod_wTrophic_gamma)) #full model best

# AICc lowest with gamma. Do that.

bestmodel_predTrophicStat <- glmmTMB(DOC ~ scaled_secchi + scaled_maxDepth + lat + 
                                       lakeTrophicStatus + elevation + lakeTrophicStatus:scaled_secchi+ (1|commonID),
                                     data = mainDF_BsM_wTrophic, family = Gamma(link = "log"))

# Diagnostics
DHARMa::simulateResiduals(bestmodel_predTrophicStat, plot=T)
hist(residuals(bestmodel_predTrophicStat), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Summary
summary(bestmodel_predTrophicStat)
r2_secchiTrophicMod <- r.squaredGLMM(bestmodel_predTrophicStat)        # R2 = 0.54 fixed effects, R2 = 0.92 w. random effects 
car::Anova(bestmodel_predTrophicStat,type = "III")

# Partial variance explained 
# Fit reduced models
model_no_secchi <- update(bestmodel_predTrophicStat, . ~ . - scaled_secchi)
model_no_depth <- update(bestmodel_predTrophicStat, . ~ . - scaled_maxDepth)
model_no_lat <- update(bestmodel_predTrophicStat, . ~ . - lat)
model_no_trophic <- update(bestmodel_predTrophicStat, . ~ . - lakeTrophicStatus)
model_no_elevation <- update(bestmodel_predTrophicStat, . ~ . - elevation)
model_no_interaction <- update(bestmodel_predTrophicStat, . ~ . - scaled_secchi:lakeTrophicStatus)


# Calculate R² for reduced models
r2_no_secchi <- r.squaredGLMM(model_no_secchi)
r2_no_depth <- r.squaredGLMM(model_no_depth)
r2_no_lat <- r.squaredGLMM(model_no_lat)
r2_no_trophic <- r.squaredGLMM(model_no_trophic)
r2_no_elevation <- r.squaredGLMM(model_no_elevation)
r2_no_interaction <- r.squaredGLMM(model_no_interaction)


# Calculate partial R²
partial_R2_secchi <- r2_secchiTrophicMod[1] - r2_no_secchi[1]
partial_R2_depth <- r2_secchiTrophicMod[1] - r2_no_depth[1]
partial_R2_lat <- r2_secchiTrophicMod[1] - r2_no_lat[1]
partial_R2_trophic <- r2_secchiTrophicMod[1] - r2_no_trophic[1]
partial_R2_elevation <- r2_secchiTrophicMod[1] - r2_no_elevation[1]
partial_R2_noInteraction <- r2_secchiTrophicMod[1] - r2_no_interaction[1]


# Print results
cat("Partial R² for scaled_secchi: ", partial_R2_secchi, "\n")
cat("Partial R² for scaled_maxDepth: ", partial_R2_depth, "\n")
cat("Partial R² for lat: ", partial_R2_lat, "\n")
cat("Partial R² for trophic: ", partial_R2_trophic, "\n")
cat("Partial R² for elevation: ", partial_R2_elevation, "\n")
cat("Partial R² for interaction: ", partial_R2_noInteraction, "\n")

# Add predicted values to df
predictedVals_trophicStat <- predict(bestmodel_predTrophicStat, type = "response", re.form=NA) #for discussion (should random effects be used or not)
#predictedVals <- predict(bestmodel_predTrophicStat, type = "response", re.form=NULL) #for discussion (should random effects be used or not)
mainDF_BsM_wTrophic$predictedDOC <- predictedVals_trophicStat
cor(mainDF_BsM_wTrophic$DOC, mainDF_BsM_wTrophic$predictedDOC)  # obs-fitted have 0.98 correlation - good
cor.test(mainDF_BsM_wTrophic$DOC, mainDF_BsM_wTrophic$predictedDOC)
r.squaredGLMM(bestmodel_predTrophicStat)
## Calculate RMSE
# Compute residuals
residualsTrophMod <- mainDF_BsM_wTrophic$DOC - predictedVals_trophicStat
# Compute RMSE
(rmseTrophicMod <- sqrt(mean(residualsTrophMod^2)))

## Relationship between observed and predicted DOC for BSM data
pObsFit_trophic <- ggplot(mainDF_BsM_wTrophic) +
  geom_point(aes(x=DOC,y=predictedDOC, colour = lat, size=maxDepth), alpha=0.7) +    # exp y-axis to get back on raw scale
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_viridis_c(option="inferno", direction = -1) +
  ylab("Model Predicted DOC (mg/L)") +
  xlab("Measured DOC (mg/L)") +
  labs(colour="Latitude", size="Max. depth (m)")+
  ylim(c(0,21))+
  xlim(c(0,21)) +
  facet_wrap(~lakeTrophicStatus)+
  theme_DOC()
pObsFit_trophic

# ggsave(filename = ("plotsTables/plots/obsFittedDOC_20241118.png"),device = "png", plot = pObsFit,
#        width = 7, height = 6, units = "in", dpi = 450)



# Estimate DOC in ARU data - using NO lake trophic status #### 
mainDF$rowID = 1:nrow(mainDF)
newdata = mainDF %>% filter(samplingProgram == "ARU")
newdata$ARU_DOC_predicted = predict(bestmodel2, newdata= newdata, type = "response", re.form = NA) #prediction only on fixed effects
#newdata$ARU_DOC_predicted_wRE = predict(bestmodel2, newdata= newdata, type = "response") #

updatedDOC = NA
#updatedDOC_RE = NA
for(i in 1:nrow(mainDF)){
  if(mainDF$samplingProgram[i] == "ARU"){ 
    updatedDOC[i] = newdata$ARU_DOC_predicted[which(newdata$rowID == i)]
    #   updatedDOC_RE[i] = newdata$ARU_DOC_predicted_wRE[which(newdata$rowID == i)]
    
  } else {
    updatedDOC[i] = mainDF$DOC[i]
    # updatedDOC_RE[i] = mainDF$DOC[i]
  }
}

mainDF$updatedDOC = updatedDOC

# Using IC approach, assemble different models to look for changes in DOC over time
M1 <- glmmTMB(updatedDOC ~ scaled_yearSample + (1|commonID),
              data = mainDF,
              family = Gamma(link = "log"))
M2 <- glmmTMB(updatedDOC ~ scaled_yearSample + (scaled_yearSample|commonID),
              data = mainDF,
              family = Gamma(link = "log"))          
M3 = update(M1, family = lognormal(link="log"))
M4 = update(M2, family = lognormal(link="log"))
AIC(M1,M2,M3,M4)
# Gamma models are best. 

# Parameterize best model for DOC over time
bestmodel3 = M2

# Diagnostics
DHARMa::simulateResiduals(bestmodel3, plot=T)
hist(residuals(bestmodel3), breaks = 30, main = "Histogram of Residuals")
# Diagnostics look fine

# Model summaries
(sumDOC <- summary(bestmodel3)) 
exp(sumDOC$coefficients$cond[, "Estimate"]) #DOC increases by 0.08% per year; 1.05 mg/L increase over 62 years (4.9% increase over time); 0.3mg/L per decade
r.squaredGLMM(bestmodel3)          # pretty much all lake effect

mainDF$DOC_fit_allData_noEnv = predict( bestmodel3, type = "response", re.form = NA)
mainDF$DOC_fitWRandom_allData_noEnv = predict( bestmodel3, type = "response", re.form = NULL)  # To plot just predicted lines

# Plot it
pDOC_allData_time <- ggplot(mainDF)+
  geom_point(aes(x=yearSample, y = updatedDOC, fill=updatedDOC), alpha=0.7, size=2.5, shape = 21) +
  geom_smooth(aes(x=yearSample, y = DOC_fit_allData_noEnv), 
              method="glm",
              formula = y~x,
              method.args = list(family = Gamma(link = 'log')))+
  scale_x_continuous(name="Year")+
  scale_y_continuous(name = "DOC (mg/L)")+
  scale_fill_gradientn(colours = c("#F5DEB3","#D2B48C","#704214","#3D2B1F")) +
  theme_DOC()+
  theme(legend.position="none")
pDOC_allData_time

# Summarize modelled differences in DOC over time period 
# Extract data from plot above - want the min/max of the fitted lines
plotData_DOC <- ggplot_build(pDOC_allData_time)
# Extract the data related to geom_smooth
lineData_DOC <- plotData_DOC$data[[2]]  # The second layer corresponds to geom_smooth

# Group by lakeTrophicStatus and calculate the min/max of the fitted values
sumDOC_allYear <- lineData_DOC %>%
  summarise(minDOC = min(y),
            maxDOC = max(y),
            diffDOC = maxDOC - minDOC,
            yearDiff = max(x)-min(x),
            diffPerDecade = (diffDOC/yearDiff)*10)

## Change in DOC and Secchi with AHI and BsM
docTimeDat <- mainDF %>% 
  group_by(commonID) %>% 
  summarise(initDOC = DOC_fitWRandom_allData_noEnv[which.min(yearSample)],
            finalDOC = DOC_fitWRandom_allData_noEnv[which.max(yearSample)],
            diffDOC = finalDOC-initDOC,
            yearDiffDOC = max(yearSample) - min(yearSample),
            initSecchi = Secchi_fitWRandom_allData_noEnv[which.min(yearSample)],
            finalSecchi = Secchi_fitWRandom_allData_noEnv[which.max(yearSample)],
            diffSecchi = finalSecchi-initSecchi,
            sampleSize = n()) %>% 
  left_join(select(statDat_wide,commonID,lat,long)) %>% 
  distinct() %>% 
  filter(!is.na(lat)) %>% 
  mutate(dataset = "all")
docTimeDat

## Map of DOC, including estimated DOC from AHI
pDOCtimeDiff <- ggplot() +
  geom_sf(data = ontario_shapefile, fill = NA, color = "black") +
  geom_point(data=docTimeDat, aes(y=lat,x=long,fill=diffDOC,size=yearDiffDOC), 
             shape = 21, stroke = 0.5, colour = "black", alpha=0.7) +
  scale_fill_gradient2(
    low = "navyblue",          
    mid = "white",         
    high = "brown",     
    midpoint = 0)+  
  ylim(c(42.5,55)) +
  labs(y="Latitude",x="Longitude",fill="Δ DOC (mg/L)",size="Sample period (years)") +
  theme_DOC() +
  theme(legend.position = "right")
pDOCtimeDiff

## NOT USED: Bootstrap Secchi results - 15 yr slope legit or not? ####

##** Below doesn't answer reviewer question as the dataset is not continuous, 
##* and for lakes with repeated measurements, does not have enough unique
##* 'starting points' to assess 15yr periods

## Bootstrap years to see if observed slope of Secchi depth from 15yr time series is stronger
# then by random chance (i.e., if years were shuffled). Randomly draw 15 year segments (i.e.,
# different 'start' times) and create a distribution of slope values

## Filter for time ranges
# 15-yr BsM data range
contempBsM_years <- sort(unique(mainDF_onlyBsM$yearSample))
contempBsM_period_length <- length(unique(contempBsM_years))
# BsM and AHI data range
all_years <- sort(unique(mainDF$yearSample))
all_period_length <- length(unique(all_years))

n_boot <- 999

## Go back to Mikes Lecture code -- figure out how to shuffle years **** ##

# Step 1: Fit model to the actual present period
present_mod <- glmmTMB(secchiDepth ~ scaled_yearSample + (scaled_yearSample | commonID),
                       data = mainDF_onlyBsM,
                       family = Gamma(link = "log"))

observed_slope <- summary(present_mod)$coefficients$cond["scaled_yearSample", "Estimate"]

# Step 2: Bootstrapping: shuffle years
boot_slopes <- numeric(n_boot)

set.seed(42)  # reproducibility
for (i in 1:n_boot) {
  # Randomly choose a starting year such that you can get a full 15-year period
  possible_start_years <- 1962:1972
  #possible_start_years <- all_years[all_years <= max(all_years) - 14]
  start_year <- sample(possible_start_years, 1)
  selected_years <- start_year:(start_year + 14)
  
  # Filter your data to the selected 15-year period
  boot_data <- mainDF %>%
    filter(yearSample %in% selected_years) %>%
    mutate(realYear = yearSample,
           yearSample_fake = scale(yearSample))  # scale the actual years
  
  # Run bootstrap model - has to just be fixed-effect model; not enough data where lakeIDs > 1 because AHI data only had 1 lake visit
  mod <- glmmTMB(secchiDepth ~ yearSample_fake,
                 data = boot_data,
                 family = Gamma(link = "log"))
  
  # collect bootstrap slopes
  boot_slopes[i] <- summary(mod)$coefficients$cond["yearSample_fake", "Estimate"]
}

boot_slopes <- na.omit(boot_slopes)
# Step 3: Plot results
df_boot <- data.frame(slope = boot_slopes)
ggplot(df_boot, aes(x = slope)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = observed_slope, color = "red", size = 1.2, linetype = "dashed") +
  labs(title = "Bootstrapped Slopes vs Observed Slope",
       x = "Slope (scaled yearSample)", y = "Frequency") +
  theme_minimal()

# Optional: empirical p-value
empirical_p <- mean(abs(boot_slopes) >= abs(observed_slope))
cat("Empirical p-value:", empirical_p, "\n")


# Step 2b: Bootstrapping: shuffle years randomly (not 15 consecutive years)
boot_slopes <- numeric(n_boot)

set.seed(42)  # reproducibility
for (i in 1:n_boot) {
  # Randomly choose a starting year such that you can get a full 15-year period
  #possible_start_years <- 1962:1972
  possible_start_years <- all_years[all_years <= max(all_years) - 16]
  #selected_years <- sample(possible_start_years,15)  # use this one if don't want any BsM
  selected_years <- sample(all_years,15)   # use this one if want to use all years (including BsM)
  
  # Filter your data to the selected 15-year period
  boot_data <- mainDF %>%
    filter(yearSample %in% selected_years) %>%
    mutate(realYear = yearSample,
           yearSample_fake = scale(yearSample))  # scale the actual years
  
  # Run bootstrap model - has to just be fixed-effect model; not enough data where lakeIDs > 1 because AHI data only had 1 lake visit
  mod <- glmmTMB(secchiDepth ~ yearSample_fake,
                 data = boot_data,
                 family = Gamma(link = "log"))
  
  # collect bootstrap slopes
  boot_slopes[i] <- summary(mod)$coefficients$cond["yearSample_fake", "Estimate"]
}

boot_slopes <- na.omit(boot_slopes)
# Step 3: Plot results
df_boot <- data.frame(slope = boot_slopes)
ggplot(df_boot, aes(x = slope)) +
  geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
  geom_vline(xintercept = observed_slope, color = "red", size = 1.2, linetype = "dashed") +
  labs(title = "Bootstrapped Slopes vs Observed Slope",
       x = "Slope (scaled yearSample)", y = "Frequency") +
  theme_minimal()

# Optional: empirical p-value
empirical_p <- mean(abs(boot_slopes) >= abs(observed_slope))
cat("Empirical p-value:", empirical_p, "\n")
