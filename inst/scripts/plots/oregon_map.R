library(ggplot2)
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


