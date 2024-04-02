require(tidyverse)
require(GGally)
require(glmmTMB)
require(ggeffects)
require(sf)
require(here)
require(basemapR)

statDat_wide <- read_csv(here("data","BsM_WideSampleModelData_20240311.csv"))
data = statDat_wide

data_sub = data %>% select(`COLTR (TCU)`, DOC, lat, long, meanDepth, TDP, FEUT, TKN, `SSO4UR (mg/L)`, secchiDepth ) %>%
filter(TKN <2000)  

names(data_sub)[8] <- "SSO"

ggplot(data_sub, aes(y = DOC, x = meanDepth))+
  geom_point()+
  geom_smooth( method="gam")

#pairplot
ggpairs(data_sub)

data_sub = data_sub %>% select(-"TKN")
data_sub <-na.omit(data_sub)
M1 = glmmTMB(DOC ~ meanDepth + TDP + FEUT + SSO + lat,
        family = Gamma(link="log"),
        data=data_sub)

summary(M1)


##spatial patterns of residuals

data_sub$resid = residuals(M1)
sf_points <- st_as_sf(data_sub, coords = c("long", "lat"), crs = 4326)

#read in some basemaps

bm = st_read("data\\Province\\Province.shp")
bm<-st_transform(bm, crs = 4326)


# Plot the heatmap
ggplot()+
  geom_sf(data=bm, fill = "grey85")+
  geom_sf(data=filter(sf_points, resid > -15),
          aes(fill=resid),
          color = "black",
          alpha = 0.75, shape = 21,
          size = 2)+
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "brown", 
    midpoint = 0
  )+
  theme_minimal()

1##do it as a heatmap


# Define the extent of the grid
bbox <- st_as_sfc(st_bbox(sf_points))

# Create a grid for heatmap
grid <- st_make_grid(bbox, square = FALSE,
                     cellsize = c(1)) %>%
  st_sf() %>% 
  mutate(ID = row_number())
  
# spatially join grid to points, so that each point is assigned the grid ID into which it falls
pointsID <- st_join(sf_points, grid)
  
pointsID <- pointsID %>% 
  as.data.frame() %>% 
  group_by(ID) %>% 
  summarize(avg_resid = mean(resid),
            count = n())  

# join aggregated values back to your grid
grid<- left_join(grid, pointsID, by = "ID")

bbox2 = st_transform(bbox, crs=4326)


# Plot the heatmap
ggplot()+
  geom_sf(data=bm, fill = "grey")+
  geom_sf(data = grid,
          aes(fill = avg_resid),  color = "transparent") +
  scale_fill_gradient2(low = "blue", high = "brown", mid= "white", 
                       midpoint = 0, na.value = "transparent",
                       name = "Residual\nDOC") +
  theme_minimal()

##
E=ggpredict(M1)
E1 = ggpredict(M1, "meanDepth")
plot(E1)
E2 = ggpredict(M1, "TDP")
plot(E2)


plot(resid(M1)~predict(M1))

