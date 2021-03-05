# Figure 1 and 4 (and supplementary 1) Preliminaries --------------------------
# STLMM Limiting Behavior Images

# set to TRUE if you don't want to write csv's
write <- FALSE

# load ggplot2 and latex2exp
library(ggplot2)
library(latex2exp)

# Figure 1 (and supplementary 1) ----------------------------------------------

# setting parameter values
s_de <- 4
s_ie <- 4
t_de <- 6
t_ie <- 7
st_de <- 1
st_ie <- 2
total_var <- sum(s_de, s_ie, t_de, t_ie, st_de, st_ie)
s_range <- 2
t_range <- 2
rangetol <- 1 / 4
s_rangetol <- s_range + s_range * rangetol
t_rangetol <- t_range + t_range * rangetol

brspace <- 0.15
brlen <- 0.08

# epsilon tolerance
epstol <- 1e-5
axistext_size <- 25
legendtext_size <- 25
annotate_size <- 11

# setting unique distance values
h_s_seq <- seq(0, s_rangetol, length.out = 200)
h_s_seq <- append(h_s_seq, epstol, s_range / 2)
h_t_seq <- seq(0, t_rangetol, length.out = 200)
h_t_seq <- append(h_t_seq, epstol, t_range / 2)
h_s <- rep(h_s_seq, times = length(h_t_seq))
h_t <- rep(h_t_seq, each = length(h_s_seq))

# make covparam object
covparams <- make_covparam_object(
  s_de = s_de,
  s_ie = s_ie,
  t_de = t_de,
  t_ie = t_ie,
  st_de = st_de,
  st_ie = st_ie,
  s_range = s_range,
  t_range = t_range,
  stcov = "productsum"
)

# make stcovariance
sigma <- make_stcovariance(covparams, h_s, h_t, "spherical", "spherical")

# make semivariogram
gamma <- make_stsemivariogram(covparams, h_s, h_t, "spherical", "spherical")

# make data frame
data <- data.frame(
  h_s = h_s,
  h_t = h_t,
  sigma = sigma,
  gamma = gamma
)

# spatial plotting subset
s_plot_pos <- data %>%
  dplyr::filter((h_t %in% c(0, epstol, max(h_t))) & (h_s > 0))

s_plot_zero <- data %>%
  dplyr::filter((h_t %in% c(0, epstol, max(h_t))) & (h_s == 0))

# temporal plotting subset
t_plot_pos <- data %>%
  dplyr::filter((h_s %in% c(0, epstol, max(h_s))) & (h_t > 0))

t_plot_zero <- data %>%
  dplyr::filter((h_s %in% c(0, epstol, max(h_s))) & (h_t == 0))


# making plots
cov_tempx <- ggplot(t_plot_pos) +
  geom_line(mapping = aes(x = h_t, y = sigma, linetype = as.factor(h_s)), size = 1.5) +
  geom_point(t_plot_zero, mapping = aes(x = h_t, y = sigma), size = 2.5) +
  labs(x = TeX("Temporal Distance $(h_t)$"), y = "Covariance") +
  scale_linetype_discrete(name = TeX("Spatial Distance $(h_s)$"), labels = c(bquote(h[s] == 0), bquote(h[s] == 0^"+"), bquote(h[s] == infinity))) +
  scale_x_continuous(breaks = c(0, t_rangetol), labels = c(TeX("$h_t = 0$"), TeX("$h_t = \\infty$"))) +
  scale_y_continuous(breaks = c(0, total_var), labels = c(0, TeX("$\\sigma^2$"))) +
  theme(
    axis.title.x = element_text(face = "bold", size = axistext_size, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = axistext_size, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = legendtext_size),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.75, 0.8),
    legend.text = element_text(size = legendtext_size),
    legend.key.size = unit(3, "line"),
    axis.text.x = element_text(size = axistext_size, face = "bold", colour = "black"),
    axis.text.y = element_text(size = axistext_size, face = "bold", colour = "black"),
    legend.key = element_rect(fill = NA)
  ) +
  expand_limits(x = c(-3.2 * rangetol, t_rangetol + 4 * rangetol)) +
  # spatial independent
  annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = s_de + s_ie, yend = s_de + s_ie, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = s_de, yend = s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + brlen, y = s_de + s_ie, yend = s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + 2 * brlen, y = s_de + s_ie / 2, yend = s_de + s_ie / 2, size = 1.5) +
  annotate("text", x = t_rangetol + brspace + 5 * brlen, y = s_de + s_ie / 2, label = TeX("$\\sigma^2_{\\delta}$"), size = annotate_size, color = "black") +
  # spatial dependent
  annotate("segment", x = t_rangetol + brspace + 2 * brlen, xend = t_rangetol + brspace + brlen + 2 * brlen, y = s_de, yend = s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + 2 * brlen, xend = t_rangetol + brspace + brlen + 2 * brlen, y = 0, yend = 0, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen + 2 * brlen, xend = t_rangetol + brspace + brlen + 2 * brlen, y = s_de, yend = 0, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen + 2 * brlen, xend = t_rangetol + brspace + 2 * brlen + 2 * brlen, y = s_de / 2, yend = s_de / 2, size = 1.5) +
  annotate("text", x = t_rangetol + brspace + 5 * brlen + 2 * brlen, y = s_de / 2, label = TeX("$\\sigma^2_{\\gamma}$"), size = annotate_size, color = "black") +
  # temporal  independent
  annotate("segment", x = -brspace, xend = -(brspace + brlen), y = t_de + t_ie, yend = t_de + t_ie, size = 1.5) +
  annotate("segment", x = -brspace, xend = -(brspace + brlen), y = t_de, yend = t_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = -(brspace + brlen), y = t_de + t_ie, yend = t_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = 0 - (brspace + 2 * brlen), y = t_de + t_ie / 2, yend = t_de + t_ie / 2, size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen), y = t_de + t_ie / 2, label = TeX("$\\sigma^2_{\\tau}$"), size = annotate_size, color = "black") +
  # temporal dependent
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = t_de, yend = t_de, size = 1.5) +
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = 0, yend = 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = t_de, yend = 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + 2 * brlen + 2 * brlen), y = t_de / 2, yend = t_de / 2, size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen + 2 * brlen), y = t_de / 2, label = TeX("$\\sigma^2_{\\eta}$"), size = annotate_size, color = "black")


if (write) {
  ggsave(plot = cov_tempx, width = 9, height = 7, units = "in", filename = "inst/images/cov_tempx.jpeg", dpi = 1200)
}



cov_spx <- ggplot(s_plot_pos) +
  geom_line(mapping = aes(x = h_s, y = sigma, linetype = as.factor(h_t)), size = 1.5) +
  geom_point(s_plot_zero, mapping = aes(x = h_s, y = sigma), size = 2.5) +
  labs(x = TeX("Spatial Distance $(h_s)$"), y = "Covariance") +
  scale_linetype_discrete(name = TeX("Temporal Distance $(h_t)$"), labels = c(bquote(h[t] == 0), bquote(h[t] == 0^"+"), bquote(h[t] == infinity))) +
  scale_x_continuous(breaks = c(0, s_rangetol), labels = c(TeX("$h_s = 0$"), TeX("$h_s = \\infty$"))) +
  scale_y_continuous(breaks = c(0, total_var), labels = c(0, TeX("$\\sigma^2$"))) +
  theme(
    axis.title.x = element_text(face = "bold", size = axistext_size, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = axistext_size, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = legendtext_size),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.75, 0.8),
    legend.text = element_text(size = legendtext_size),
    legend.key.size = unit(3, "line"),
    axis.text.x = element_text(size = axistext_size, face = "bold", colour = "black"),
    axis.text.y = element_text(size = axistext_size, face = "bold", colour = "black"),
    legend.key = element_rect(fill = NA)
  ) +
  expand_limits(x = c(-3.2 * rangetol, s_rangetol + 4 * rangetol)) +
  # spatial independent
  annotate("segment", x = -(brspace), xend = -(brspace + brlen), y = s_de + s_ie, yend = s_de + s_ie, size = 1.5) +
  annotate("segment", x = -(brspace), xend = -(brspace + brlen), y = s_de, yend = s_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = -(brspace + brlen), y = s_de + s_ie, yend = s_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = -(brspace + 2 * brlen), y = s_de + s_ie / 2, yend = s_de + s_ie / 2, size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen), y = s_de + s_ie / 2, label = TeX("$\\sigma^2_{\\delta}$"), size = annotate_size, color = "black") +
  # spatial dependent
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = s_de, yend = s_de, size = 1.5) +
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = 0, yend = 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = s_de, yend = 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + 2 * brlen + 2 * brlen), y = s_de / 2, yend = s_de / 2, size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen + 2 * brlen), y = s_de / 2, label = TeX("$\\sigma^2_{\\gamma}$"), size = annotate_size, color = "black") +
  # temporal independent
  annotate("segment", x = s_rangetol + brspace, xend = s_rangetol + (brspace + brlen), y = t_de + t_ie, yend = t_de + t_ie, size = 1.5) +
  annotate("segment", x = s_rangetol + brspace, xend = s_rangetol + (brspace + brlen), y = t_de, yend = t_de, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen), xend = s_rangetol + (brspace + brlen), y = t_de + t_ie, yend = t_de, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen), xend = s_rangetol + (brspace + 2 * brlen), y = t_de + t_ie / 2, yend = t_de + t_ie / 2, size = 1.5) +
  annotate("text", x = s_rangetol + (brspace + 5 * brlen), y = t_de + t_ie / 2, label = TeX("$\\sigma^2_{\\tau}$"), size = annotate_size, color = "black") +
  # temporal dependent
  annotate("segment", x = s_rangetol + (brspace + 2 * brlen), xend = s_rangetol + (brspace + brlen + 2 * brlen), y = t_de, yend = t_de, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + 2 * brlen), xend = s_rangetol + (brspace + brlen + 2 * brlen), y = 0, yend = 0, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen + 2 * brlen), xend = s_rangetol + (brspace + brlen + 2 * brlen), y = t_de, yend = 0, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen + 2 * brlen), xend = s_rangetol + (brspace + 2 * brlen + 2 * brlen), y = t_de / 2, yend = t_de / 2, size = 1.5) +
  annotate("text", x = s_rangetol + (brspace + 5 * brlen + 2 * brlen), y = t_de / 2, label = TeX("$\\sigma^2_{\\eta}$"), size = annotate_size, color = "black")

if (write) {
  ggsave(plot = cov_spx, width = 9, height = 7, units = "in", filename = "inst/images/cov_spx.jpeg", dpi = 1200)
}


sv_tempx <- ggplot(t_plot_pos) +
  geom_line(mapping = aes(x = h_t, y = gamma, linetype = as.factor(h_s)), size = 1.5) +
  geom_point(t_plot_zero, mapping = aes(x = h_t, y = gamma), size = 2.5) +
  labs(x = TeX("Temporal Distance $(h_t)$"), y = "Semivariance") +
  scale_linetype_discrete(name = TeX("Spatial Distance $(h_s)$"), labels = c(bquote(h[s] == 0), bquote(h[s] == 0^"+"), bquote(h[s] == infinity))) +
  scale_x_continuous(breaks = c(0, t_rangetol), labels = c(TeX("$h_t = 0$"), TeX("$h_t = \\infty$"))) +
  scale_y_continuous(breaks = c(0, total_var), labels = c(0, TeX("$\\sigma^2$"))) +
  theme(
    axis.title.x = element_text(face = "bold", size = axistext_size, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = axistext_size, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = legendtext_size),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.75, 0.25),
    legend.text = element_text(size = legendtext_size),
    legend.key.size = unit(3, "line"),
    axis.text.x = element_text(size = axistext_size, face = "bold", colour = "black"),
    axis.text.y = element_text(size = axistext_size, face = "bold", colour = "black"),
    legend.key = element_rect(fill = NA)
  ) +
  expand_limits(x = c(-3.2 * rangetol, t_rangetol + 4 * rangetol), y = c(-20 * rangetol, NA)) +
  # spatial independent
  annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = total_var - (s_de + s_ie), yend = total_var - (s_de + s_ie), size = 1.5) +
  annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = total_var - s_de, yend = total_var - s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + brlen, y = total_var - (s_de + s_ie), yend = total_var - s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + 2 * brlen, y = total_var - (s_de + s_ie / 2), yend = total_var - (s_de + s_ie / 2), size = 1.5) +
  annotate("text", x = t_rangetol + brspace + 5 * brlen, y = total_var - (s_de + s_ie / 2), label = TeX("$\\sigma^2_{\\delta}$"), size = annotate_size, color = "black") +
  # spatial dependent
  annotate("segment", x = t_rangetol + brspace + 2 * brlen, xend = t_rangetol + brspace + brlen + 2 * brlen, y = total_var - s_de, yend = total_var - s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + 2 * brlen, xend = t_rangetol + brspace + brlen + 2 * brlen, y = total_var - 0, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen + 2 * brlen, xend = t_rangetol + brspace + brlen + 2 * brlen, y = total_var - s_de, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = t_rangetol + brspace + brlen + 2 * brlen, xend = t_rangetol + brspace + 2 * brlen + 2 * brlen, y = total_var - s_de / 2, yend = total_var - s_de / 2, size = 1.5) +
  annotate("text", x = t_rangetol + brspace + 5 * brlen + 2 * brlen, y = total_var - s_de / 2, label = TeX("$\\sigma^2_{\\gamma}$"), size = annotate_size, color = "black") +
  # temporal  independent
  annotate("segment", x = -brspace, xend = -(brspace + brlen), y = total_var - (t_de + t_ie), yend = total_var - (t_de + t_ie), size = 1.5) +
  annotate("segment", x = -brspace, xend = -(brspace + brlen), y = total_var - t_de, yend = total_var - t_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = -(brspace + brlen), y = total_var - (t_de + t_ie), yend = total_var - t_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = 0 - (brspace + 2 * brlen), y = total_var - (t_de + t_ie / 2), yend = total_var - (t_de + t_ie / 2), size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen), y = total_var - (t_de + t_ie / 2), label = TeX("$\\sigma^2_{\\tau}$"), size = annotate_size, color = "black") +
  # temporal dependent
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = total_var - t_de, yend = total_var - t_de, size = 1.5) +
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = total_var - 0, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = total_var - t_de, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + 2 * brlen + 2 * brlen), y = total_var - t_de / 2, yend = total_var - t_de / 2, size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen + 2 * brlen), y = total_var - t_de / 2, label = TeX("$\\sigma^2_{\\eta}$"), size = annotate_size, color = "black")


if (write) {
  ggsave(plot = sv_tempx, width = 9, height = 7, units = "in", filename = "inst/images/sv_tempx.jpeg", dpi = 1200)
}



sv_spx <- ggplot(s_plot_pos) +
  geom_line(mapping = aes(x = h_s, y = gamma, linetype = as.factor(h_t)), size = 1.5) +
  geom_point(s_plot_zero, mapping = aes(x = h_s, y = gamma), size = 2.5) +
  labs(x = TeX("Spatial Distance $(h_s)$"), y = "Semivariance") +
  scale_linetype_discrete(name = TeX("Temporal Distance $(h_t)$"), labels = c(bquote(h[t] == 0), bquote(h[t] == 0^"+"), bquote(h[t] == infinity))) +
  scale_x_continuous(breaks = c(0, s_rangetol), labels = c(TeX("$h_s = 0$"), TeX("$h_s = \\infty$"))) +
  scale_y_continuous(breaks = c(0, total_var), labels = c(0, TeX("$\\sigma^2$"))) +
  theme(
    axis.title.x = element_text(face = "bold", size = axistext_size, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = axistext_size, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = legendtext_size),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.75, 0.25),
    legend.text = element_text(size = legendtext_size),
    legend.key.size = unit(3, "line"),
    axis.text.x = element_text(size = axistext_size, face = "bold", colour = "black"),
    axis.text.y = element_text(size = axistext_size, face = "bold", colour = "black"),
    legend.key = element_rect(fill = NA)
  ) +
  expand_limits(x = c(-3.2 * rangetol, s_rangetol + 4 * rangetol), y = c(-20 * rangetol, NA)) +
  # spatial independent
  annotate("segment", x = -(brspace), xend = -(brspace + brlen), y = total_var - (s_de + s_ie), yend = total_var - (s_de + s_ie), size = 1.5) +
  annotate("segment", x = -(brspace), xend = -(brspace + brlen), y = total_var - s_de, yend = total_var - s_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = -(brspace + brlen), y = total_var - (s_de + s_ie), yend = total_var - s_de, size = 1.5) +
  annotate("segment", x = -(brspace + brlen), xend = -(brspace + 2 * brlen), y = total_var - (s_de + s_ie / 2), yend = total_var - (s_de + s_ie / 2), size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen), y = total_var - (s_de + s_ie / 2), label = TeX("$\\sigma^2_{\\delta}$"), size = annotate_size, color = "black") +
  # spatial dependent
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = total_var - s_de, yend = total_var - s_de, size = 1.5) +
  annotate("segment", x = -(brspace + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = total_var - 0, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + brlen + 2 * brlen), y = total_var - s_de, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = -(brspace + brlen + 2 * brlen), xend = -(brspace + 2 * brlen + 2 * brlen), y = total_var - s_de / 2, yend = total_var - s_de / 2, size = 1.5) +
  annotate("text", x = -(brspace + 5 * brlen + 2 * brlen), y = total_var - s_de / 2, label = TeX("$\\sigma^2_{\\gamma}$"), size = annotate_size, color = "black") +
  # temporal independent
  annotate("segment", x = s_rangetol + brspace, xend = s_rangetol + (brspace + brlen), y = total_var - (t_de + t_ie), yend = total_var - (t_de + t_ie), size = 1.5) +
  annotate("segment", x = s_rangetol + brspace, xend = s_rangetol + (brspace + brlen), y = total_var - t_de, yend = total_var - t_de, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen), xend = s_rangetol + (brspace + brlen), y = total_var - (t_de + t_ie), yend = total_var - t_de, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen), xend = s_rangetol + (brspace + 2 * brlen), y = total_var - (t_de + t_ie / 2), yend = total_var - (t_de + t_ie / 2), size = 1.5) +
  annotate("text", x = s_rangetol + (brspace + 5 * brlen), y = total_var - (t_de + t_ie / 2), label = TeX("$\\sigma^2_{\\tau}$"), size = annotate_size, color = "black") +
  # temporal dependent
  annotate("segment", x = s_rangetol + (brspace + 2 * brlen), xend = s_rangetol + (brspace + brlen + 2 * brlen), y = total_var - t_de, yend = total_var - t_de, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + 2 * brlen), xend = s_rangetol + (brspace + brlen + 2 * brlen), y = total_var - 0, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen + 2 * brlen), xend = s_rangetol + (brspace + brlen + 2 * brlen), y = total_var - t_de, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = s_rangetol + (brspace + brlen + 2 * brlen), xend = s_rangetol + (brspace + 2 * brlen + 2 * brlen), y = total_var - t_de / 2, yend = total_var - t_de / 2, size = 1.5) +
  annotate("text", x = s_rangetol + (brspace + 5 * brlen + 2 * brlen), y = total_var - t_de / 2, label = TeX("$\\sigma^2_{\\eta}$"), size = annotate_size, color = "black")

if (write) {
  ggsave(plot = sv_spx, width = 9, height = 7, units = "in", filename = "inst/images/sv_spx.jpeg", dpi = 1200)
}

# Figure 4  -------------------------------------------------------------------
# Fitted REML Semivariogram


covparams <- read.csv("inst/output/dataanalysis/covparams.csv")
ps_reml <- covparams %>%
  dplyr::filter(estmethod == "reml" & stcov == "productsum") %>%
  select(covparam, value) %>%
  pivot_wider(names_from = covparam)
s_de <- ps_reml$s_de
s_ie <- ps_reml$s_ie
t_de <- ps_reml$t_de
t_ie <- ps_reml$t_ie
st_de <- ps_reml$st_de
st_ie <- ps_reml$st_ie
total_var <- sum(s_de, s_ie, t_de, t_ie, st_de, st_ie)
s_range <- ps_reml$s_range
t_range <- ps_reml$t_range
rangetol <- 1 / 4
s_rangetol <- s_range + s_range * rangetol
# t_rangetol <- t_range + t_range * rangetol
t_rangetol <- t_range * 2


brspace <- 0.075
brlen <- 0.15

# epsilon tolerance
epstol <- 1e-5
axistext_size <- 24
legendtext_size <- 24
annotate_size <- 9

s_epstol <- c(10, 750, 1500, Inf)

# setting unique distance values
h_s_seq <- seq(0, s_rangetol, length.out = 200)
h_s_seq <- append(h_s_seq, s_epstol, s_range / 2)
h_t_seq <- seq(0, t_rangetol, length.out = 200)
h_t_seq <- append(h_t_seq, t_range / 2)
h_s <- rep(h_s_seq, times = length(h_t_seq))
h_t <- rep(h_t_seq, each = length(h_s_seq))

# make covparam object
covparams <- make_covparam_object(
  s_de = s_de,
  s_ie = s_ie,
  t_de = t_de,
  t_ie = t_ie,
  st_de = st_de,
  st_ie = st_ie,
  s_range = s_range,
  t_range = t_range,
  stcov = "productsum"
)

# make stcovariance
sigma <- make_stcovariance(covparams, h_s, h_t, "exponential", "exponential")

# make semivariogram
gamma <- make_stsemivariogram(covparams, h_s, h_t, "exponential", "exponential")

# make data frame
data <- data.frame(
  h_s = h_s,
  h_t = h_t,
  sigma = sigma,
  gamma = gamma
)


# temporal plotting subset
t_plot_pos <- data %>%
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t > 0.2))

t_plot_zero <- data %>%
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t == 0))

plot_s_de <- total_var - t_plot_pos$gamma[t_plot_pos$h_s == 10 & t_plot_pos$h_t == max(t_plot_pos$h_t)]
plot_s_ie <- t_plot_pos$gamma[t_plot_pos$h_s == 10 & t_plot_pos$h_t == max(t_plot_pos$h_t)] -
  t_plot_pos$gamma[t_plot_pos$h_s == 0 & t_plot_pos$h_t == max(t_plot_pos$h_t)]

sv_tempx_ps_reml <- ggplot(t_plot_pos) +
  geom_line(mapping = aes(x = h_t, y = gamma, linetype = as.factor(h_s)), size = 1.5) +
  geom_point(t_plot_zero, mapping = aes(x = h_t, y = gamma), size = 4) +
  labs(x = TeX("Temporal Distance $(h_t)$"), y = "Semivariance") +
  scale_linetype_discrete(name = TeX("Spatial Distance $(h_s)$ in km"), labels = c(bquote(h[s] == 0), bquote(h[s] == 0^"+"), bquote(h[s] == 750), bquote(h[s] == 1500), bquote(h[s] == infinity))) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, TeX("$\\infty$"))) +
  scale_y_continuous(breaks = c(0, 30, 60, 90, 120), labels = c(0, 30, 60, 90, 120)) +
  theme(
    axis.title.x = element_text(face = "bold", size = axistext_size, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = axistext_size, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = legendtext_size),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.5, 0.5),
    legend.text = element_text(size = legendtext_size),
    legend.key.size = unit(3, "line"),
    axis.text.x = element_text(size = axistext_size, face = "bold", colour = "black"),
    axis.text.y = element_text(size = axistext_size, face = "bold", colour = "black"),
    legend.key = element_rect(fill = NA),
    legend.key.height = unit(1, "cm")
  ) +
  expand_limits(x = c(-1 * rangetol, t_rangetol + 5.5 * rangetol), y = c(0, 120)) +
  # spatial independent
  # annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = total_var - (plot_s_de + plot_s_ie), yend =  total_var - (plot_s_de + plot_s_ie), size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = total_var - plot_s_de, yend =  total_var - plot_s_de, size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + brlen, y = total_var - (plot_s_de + plot_s_ie), yend =  total_var - plot_s_de, size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + 2 * brlen, y = total_var - (plot_s_de + plot_s_ie / 2), yend =  total_var - (plot_s_de + plot_s_ie / 2), size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace + 2 * brlen, xend = t_rangetol + brspace + 2 * brlen, y = total_var - (plot_s_de + plot_s_ie / 2), yend =  total_var - (plot_s_de + plot_s_ie / 2) - 40 * brlen, size = 1.5) +
  # annotate("text", x = t_rangetol + brspace + 2.75 * brlen, y = total_var - (plot_s_de + plot_s_ie / 2) - 80 * brlen, label = TeX("$\\sigma^2_{\\gamma}$"), size = annotate_size, color = "black") +
  # spatial dependent
  annotate("segment", x = t_rangetol + 4 * brspace + 2 * brlen, xend = t_rangetol + 4 * brspace + brlen + 2 * brlen, y = total_var - plot_s_de, yend = total_var - plot_s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + 4 * brspace + 2 * brlen, xend = t_rangetol + 4 * brspace + brlen + 2 * brlen, y = total_var - 0, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = t_rangetol + 4 * brspace + brlen + 2 * brlen, xend = t_rangetol + 4 * brspace + brlen + 2 * brlen, y = total_var - plot_s_de, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = t_rangetol + 4 * brspace + brlen + 2 * brlen, xend = t_rangetol + 4 * brspace + 2.75 * brlen + 2 * brlen, y = total_var - plot_s_de / 2, yend = total_var - plot_s_de / 2, size = 1.5) +
  annotate("text", x = t_rangetol + brspace + 5 * brlen + 5 * brlen, y = total_var - plot_s_de / 2, label = TeX("$\\sigma^2_{\\delta}$"), size = annotate_size, color = "black")

if (write) {
  ggsave(plot = sv_tempx_ps_reml, width = 9, height = 7, units = "in", filename = "inst/images/sv_tempx_ps_reml.jpeg", dpi = 1200)
}
