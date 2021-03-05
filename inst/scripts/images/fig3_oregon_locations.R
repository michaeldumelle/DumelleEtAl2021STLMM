# Figure 3 Preliminaries -----------------------------------------
# Oregon Locations Image

# set to TRUE if you don't want to write csv's
write <- FALSE

# load the package
library(DumelleEtAl2021STLMM)

# load ggplot2
library(ggplot2)

# load the data
data("or_data")

# subset into training and test
or_train <- subset(or_data, TYPE == "TRAIN")
or_test <- subset(or_data, TYPE == "TEST")

# load us states data
us_states <- map_data("state")
oregon <- us_states[us_states$region == "oregon", ]
# project their coordinates
proj_coords <- LLtoUTM(cm = -121.5158, lat = oregon$lat, lon = oregon$long)
# store as a data frame for plotting
map_oregon <- data.frame(LONGITUDE = proj_coords$xy[, "x"], LATITUDE = proj_coords$xy[, "y"])
# scale accounting for likely different datum
map_oregon$LONGITUDE <- map_oregon$LONGITUDE - 8

# Oregon Locations Image (Figure 3) -------------------------------------------------
or_locations <- ggplot() +
  geom_path(map_oregon, mapping = aes(x = LONGITUDE, y = LATITUDE), col = "grey", size = 6) +
  geom_point(data = or_test, mapping = aes(x = LONGITUDE, y = LATITUDE, color = "Test Data"), pch = 4, size = 5) +
  geom_point(data = or_train, mapping = aes(x = LONGITUDE, y = LATITUDE, color = "Training Data"), pch = 19, size = 6) +
  scale_color_manual(values = c("black", "black"), name = "Spatial Locations") +
  guides(colour = guide_legend(override.aes = list(shape = c(19, 4)))) +
  expand_limits(y = c(0, 600)) +
  theme(
    legend.position = c(0.44, 0.88),
    legend.text = element_text(size = 23),
    legend.title = element_text(size = 23, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.key.size = unit(2.5, "line"),
    panel.background = element_blank()
  )

if (write) {
  ggsave(plot = or_locations, width = 9, height = 7, units = "in", filename = "inst/images/or_locations.jpeg", dpi = 1200)
}
