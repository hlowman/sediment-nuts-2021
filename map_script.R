# Figure 1 Map Creation
# Heili Lowman
# 5/25/21

# The following script will create the map to be used in the 
# sediment nutrient recycling manuscript.

# Load packages
library(tidyverse) 
library(patchwork) 
library(ggmap)
library(ggrepel)
library(sf)
library(USAboundaries)
library(USAboundariesData)
library(ggspatial)

# Load data
map_data <- read_csv("data_raw/sites.csv")

# Create data sf object
map_sf <- st_as_sf(map_data,
                   coords = c("Long", "Lat"),
                   remove = F,
                   crs = 4326) %>% # WGS 84 projection
  mutate(site_f = factor(Site, levels = c("REFU", "GOLP", "GOLB", "ABUR", "MOHK", "MICR"))) %>%
  mutate(n_x = c(0.0,0.0,0.0,-0.02,0.01,0.0)) %>%
  mutate(n_y = c(-0.02,-0.02,-0.02,-0.02,-0.02,-0.02))
  
# Base plot to see how things are looking
plot(map_sf$geometry)

# create base terrain map tile
# create bounding box
lats <- c(34.35, 34.5)
lons <- c(-120.15, -119.6)
bb <- make_bbox(lon = lons, lat = lats, f = 0.05)

sb_basemap <- get_stamenmap(bb,
                      zoom = 12,
                      maptype = 'terrain-background')

ggmap(sb_basemap)

# create california state map for inset
states <- us_states()
ca <- states %>%
  dplyr::filter(name %in% c('California'))
ca <- st_transform(ca, 4326)

camap <- ggplot(ca) +
  geom_sf() +
  annotate("rect", xmin = -120.2, xmax = -119.4, 
           ymin = 34.2, ymax = 34.6,
           alpha = 0, color = 'black') + # adds box for zoom in definition
  theme_classic() + 
  theme(text=element_text(family="Times New Roman", size=14),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) + # removes ALL axes and gridlines
  coord_sf(crs = st_crs(4326))

camap

# create full color map for displaying Goleta Pier sampling site
fullmap <- ggmap(sb_basemap) + # base google maps tile
  geom_point(data = map_sf, aes(x = Long, y = Lat),
             size = 3,
             shape = 16,
             inherit.aes = FALSE) + # adds points
  geom_text_repel(data = map_sf, 
                  aes(x = Long, 
                      y = Lat, 
                      label = site_f),
                  nudge_x = map_sf$n_x,
                  nudge_y = map_sf$n_y,
                  segment.size = 0.2,
                  size = 5) +
  ggspatial::annotation_north_arrow(location = "tr", 
                                    height = unit(1, "cm"),
                                    width = unit(1, "cm")) + # adds compass due north
  ggspatial::annotation_scale() + # adds scale
  geom_text(x = -120, y = 34.38, label = "Santa Barbara Channel", color = "gray40", size = 4, fontface = "italic") +
  labs(x = "Longitude (WGS84)",
       y = "Latitude") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=14)) +
  coord_sf(crs = st_crs(4326))

fullmap

# Now use patchwork to knit them together
side <- camap + fullmap + plot_annotation(tag_levels = 'A')
side # side-by-side plot

# Export map to desktop.
# ggsave(("Figure_1.tiff"),
#        path = "figures",
#        width = 25,
#        height = 10,
#        units = "cm"
#        )

# inside <- fullmap + inset_element(camap, left = 0.75, bottom = 0.05, right = 0.95, top = 0.35)
# inside # inset plot

# End of script.
