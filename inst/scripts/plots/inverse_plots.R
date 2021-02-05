library(ggplot2) # load ggplot2
path <- "C:/Users/mdumelle/Documents/publications/DumelleEtAl2020STLMM/submission/SpatialStatistics/resubmission/GitHub/DumelleEtAl2021STLMM/inst/scripts/inverses"
inverse_avg <- read.csv(paste0(path, "/inverse_avg.csv"))
inverse_ratios <- read.csv(paste0(path, "/inverse_ratios.csv"))
stempsv_avg <- read.csv(paste0(path, "/stempsv_avg.csv"))
stempsv_avg$algorithm <- "stempsv"

write <- TRUE

if (write) {
  pathsave <- "C:/Users/mdumelle/Documents/publications/DumelleEtAl2020STLMM/submission/SpatialStatistics/resubmission/GitHub/DumelleEtAl2021STLMM/inst/plots"
}

# raw inverse image
inverse_raw <- ggplot(inverse_avg, mapping = aes(x = n_st, y = seconds, shape = algorithm)) +
  geom_line(lty = "solid", alpha = 0.1, size = 3) +
  geom_point(size = 6) +
  labs(x = "Sample Size (Thousands)", y = "Inversion Seconds") +
  scale_x_continuous(
    limits = c(0, 15000),
    breaks = c(0, 5000, 10000, 15000),
    labels = c(0, 5, 10, 15)
  ) +
  scale_y_continuous(
    limits = c(0, 2000),
    breaks = seq(0, 2000, length.out = 5),
    labels = seq(0, 2000, length.out = 5)
  ) +
  scale_shape_manual(
    name = "Algorithm",
    values = c(8, 19, 17, 15),
    breaks = c("cholesky", "productsum", "sum_with_error", "product"),
    labels = c("CHOL", "PS", "SWE", "P")
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 30, face = "bold"),
    axis.title = element_text(face = "bold", size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(face = "bold", size = 30),
    legend.title = element_text(face = "bold", size = 30),
    legend.position = c(0.2, 0.8)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 7)))

if (write) {
  ggsave(plot = inverse_raw, width = 9, height = 7, units = "in", filename = paste0(pathsave, "/inverse_raw.jpeg"), dpi = 1200)
}


# zoomed inverse image
inverse_zoom <- ggplot(inverse_avg[inverse_avg$algorithm != "cholesky", ], mapping = aes(x = n_st, y = seconds, shape = algorithm)) +
  geom_line(lty = "solid", alpha = 0.1, size = 3) +
  geom_point(size = 6) +
  coord_cartesian(xlim = c(12000, 15000), ylim = c(0, 90)) +
  labs(x = "Sample Size (Thousands)", y = "Inversion Seconds") +
  scale_x_continuous(
    limits = c(0, 15000),
    breaks = c(seq(12000, 15000, by = 1000)),
    labels = c(seq(12, 15, by = 1))
  ) +
  scale_y_continuous(
    limits = c(0, 2000),
    breaks = seq(0, 80, length.out = 5),
    labels = seq(0, 80, length.out = 5)
  ) +
  scale_shape_manual(
    name = "Algorithm",
    values = c(19, 17, 15),
    breaks = c("productsum", "sum_with_error", "product"),
    labels = c("PS", "SWE", "P")
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 30, face = "bold"),
    axis.title = element_text(face = "bold", size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(face = "bold", size = 30),
    legend.title = element_text(face = "bold", size = 30),
    legend.position = c(0.25, 0.85)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 7)))

inverse_zoom
if (write) {
  ggsave(plot = inverse_zoom, width = 9, height = 7, units = "in", filename = paste0(pathsave, "/inverse_zoom.jpeg"), dpi = 1200)
}

# ratio inverse image
inverse_ratio <- ggplot(inverse_ratios, mapping = aes(x = n_st, y = ratio, shape = algorithm)) +
  geom_line(lty = "solid", alpha = 0.1, size = 3) +
  geom_point(size = 6) +
  labs(x = "Sample Size (Thousands)", y = "Inversion Ratios") +
  scale_x_continuous(
    limits = c(0, 15000),
    breaks = c(seq(0, 15000, by = 5000)),
    labels = c(seq(0, 15, by = 5))
  ) +
  scale_y_continuous(
    limits = c(0, 150),
    breaks = seq(0, 150, length.out = 6),
    labels = seq(0, 150, length.out = 6)
  ) +
  scale_shape_manual(
    name = "Algorithm",
    values = c(15, 17, 19),
    breaks = c("pratio", "sweratio", "psratio"),
    labels = c("CHOL / P", "CHOL / SWE", "CHOL / PS")
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 30, face = "bold"),
    axis.title = element_text(face = "bold", size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(face = "bold", size = 30),
    legend.title = element_text(face = "bold", size = 30),
    legend.position = c(0.25, 0.85)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 7)))

inverse_ratio
if (write) {
  ggsave(plot = inverse_ratio, width = 9, height = 7, units = "in", filename = paste0(pathsave, "/inverse_ratio.jpeg"), dpi = 1200)
}

# semivariogram image
stempsv_raw <- ggplot(stempsv_avg, mapping = aes(x = n_st, y = seconds)) +
  geom_line(lty = "solid", alpha = 0.1, size = 3) +
  geom_point(size = 6, shape = 18) +
  labs(x = "Sample Size (Thousands)", y = "Semivariogram Seconds") +
  scale_x_continuous(
    limits = c(0, 15000),
    breaks = c(seq(0, 15000, by = 5000)),
    labels = c(seq(0, 15, by = 5))
  ) +
  scale_y_continuous(
    limits = c(0, 50),
    breaks = seq(0, 50, length.out = 6),
    labels = seq(0, 50, length.out = 6)
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 30, face = "bold"),
    axis.title = element_text(face = "bold", size = 30),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
  )

stempsv_raw
if (write) {
  ggsave(plot = stempsv_raw, width = 9, height = 7, units = "in", filename = paste0(pathsave, "/stempsv_raw.jpeg"), dpi = 1200)
}
