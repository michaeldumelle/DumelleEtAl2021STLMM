library(ggplot2)
write <- TRUE
data("or_train")
or_train <- unique(or_train[, c("LATITUDE", "LONGITUDE", "STATION")])
or_train$set <- "train"
data("or_test")
or_test <- unique(or_test[, c("LATITUDE", "LONGITUDE", "STATION")])
or_test$set <- "test"
or <- rbind(or_train, or_test)
us_states <- map_data("state")
oregon <- us_states[us_states$region == "oregon", ]
proj_coords <- LLtoUTM(cm = -121.5158, lat = oregon$lat, lon = oregon$long)
map_oregon <- data.frame(LONGITUDE = proj_coords$xy[, "x"], LATITUDE = proj_coords$xy[, "y"])
map_oregon$LONGITUDE <- map_oregon$LONGITUDE - 8


or_locations <- ggplot() +
  geom_path(map_oregon, mapping = aes(x = LONGITUDE, y = LATITUDE), col = "grey", size = 6) +
  geom_point(data = or_test, mapping = aes(x = LONGITUDE, y = LATITUDE, color = "Test Data"), pch = 4, size = 5) +
  geom_point(data = or_train, mapping = aes(x = LONGITUDE, y = LATITUDE, color = "Training Data"), pch = 19, size = 6) +
  scale_color_manual(values = c("black", "black"), name = "Spatial Locations") +
  guides(colour = guide_legend(override.aes = list(shape = c(19, 4)))) +
  expand_limits(y = c(0, 600)) +
  theme(legend.position = c(0.44, 0.88),
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 23, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.key.size = unit(2.5,"line"),
        panel.background = element_blank())

if (write) {
  library(readr)
  pathsave <- "C:/Users/mdumelle/Documents/publications/DumelleEtAl2020STLMM/submission/SpatialStatistics/resubmission/GitHub/DumelleEtAl2021STLMM/inst/plots"
  ggsave(plot = or_locations, width = 9, height = 7, units = "in", filename = paste0(pathsave, "/or_locations.jpeg"), dpi = 1200)
}



