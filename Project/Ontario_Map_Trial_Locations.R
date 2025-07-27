# Install if not already installed
install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata", "ggspatial", "dplyr"))

# Load packages
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(dplyr)

# Load Canada provinces map (as sf object)
canada <- ne_states(country = "Canada", returnclass = "sf")

# Filter for Ontario
ontario <- canada %>% filter(name == "Ontario")

# Your dataset with latitude and longitude
# Replace with your actual file
df <- read.csv("Environmental_Summary.csv")

# Make sure your data has columns: Location, latitude, longitude
head(df)

# Plot Ontario map with points
ggplot() +
  geom_sf(data = ontario, fill = "antiquewhite", color = "black") +
  geom_point(data = df, aes(x = Longitude, y = Latitude), color = "red", size = 2) +
  coord_sf(xlim = c(-95, -74), ylim = c(41, 57), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         style = north_arrow_fancy_orienteering) +
  labs(title = "Trial Locations in Ontario",
       x = "Longitude", y = "Latitude") +
  theme_minimal()
ggsave("ontario_map.jpg", width = 10, height = 8, dpi = 300)
