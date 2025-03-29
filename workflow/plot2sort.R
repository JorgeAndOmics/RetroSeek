# ----------------------------
# PLOT GENERATION SCRIPT
# ----------------------------

# --------------------------------
# 1. LOAD LIBRARIES & SET OPTIONS
# --------------------------------
options(warn = -1)  # Suppress warnings

suppressMessages({  # Load required libraries silently
  library(argparse)     # Modern argument parser
  library(arrow)
  library(tidyverse)
  library(ggsci)
  library(ggalluvial)
  library(ggdist)
})

# --------------------------
# 2. HELPER COLOR PALETTE
# --------------------------
futurama_unlimited_palette <- function(input_colour_number = 12, output_colour_number) {
  planet_express <- pal_futurama("planetexpress")(input_colour_number)
  output_colour <- colorRampPalette(planet_express)(output_colour_number)
  return(output_colour)
}

# --------------------------
# 3. PARSE ARGUMENTS (argparse)
# --------------------------
parser <- ArgumentParser(description = "Generate plots from enERVate Parquet results")

parser$add_argument("--input", required = TRUE, help = "Input directory with *_main.parquet and *_accessory.parquet files")
parser$add_argument("--output", required = TRUE, help = "Directory to save output plots")
parser$add_argument("--threshold", required = TRUE, type = "double", help = "Filter threshold as proportion (e.g. 0.01)")

args <- parser$parse_args()

input.dir <- args$input
output.dir <- args$output
plot_filter_threshold <- args$threshold

print("Generating plots...")

# --------------------------
# 4. READ PARQUET FILES
# --------------------------
main_files <- list.files(input.dir, pattern = "_main\\.parquet$", full.names = TRUE)
accessory_files <- list.files(input.dir, pattern = "_accessory\\.parquet$", full.names = TRUE)

if (length(main_files) == 0) stop("No main parquet files found.")
if (length(accessory_files) == 0) stop("No accessory parquet files found.")

all.main <- purrr::map_dfr(main_files, arrow::read_parquet)
all.accessory <- purrr::map_dfr(accessory_files, arrow::read_parquet)
all.full <- bind_rows(all.main, all.accessory)

# --------------------------
# 5. QUANTILE CALCULATIONS
# --------------------------
q_stats <- function(df) {
  list(
    mean = mean(df$mean_bitscore),
    q1 = quantile(df$mean_bitscore, 0.25),
    median = quantile(df$mean_bitscore, 0.5),
    q3 = quantile(df$mean_bitscore, 0.75)
  )
}

bit.main <- q_stats(all.main)
bit.accessory <- q_stats(all.accessory)

# --------------------------
# 6. COUNT PROBES PER GROUP
# --------------------------
group_count <- function(df) {
  df %>%
    group_by(species_name, virus, probe, family, abbreviation) %>%
    summarise(count = n(), .groups = "keep") %>%
    ungroup()
}

all.counted_probe <- group_count(all.full)
main.counted_probe <- group_count(all.main)
accessory.counted_probe <- group_count(all.accessory)

# --------------------------
# 7. DENSITY PLOTS
# --------------------------
density_bitscore_plot <- function(data, q1, median, q3) {
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))
  ggplot(data, aes(x = mean_bitscore, fill = probe, colour = probe, weight = mean_identity)) +
    geom_density(alpha = 0.4, adjust = 3) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    labs(x = "HSP Bitscore", y = "Density") +
    geom_vline(xintercept = c(q1, median, q3), linetype = "dashed", linewidth = 0.2) +
    annotate("text", x = q1 + 20, y = 1e-4, label = "Q1") +
    annotate("text", x = median + 20, y = 1e-4, label = "Q2") +
    annotate("text", x = q3 + 20, y = 1e-4, label = "Q3") +
    scale_x_continuous(breaks = seq(0, max(data$mean_bitscore, na.rm = TRUE), by = 100)) +
    theme_minimal() +
    theme(
      text = element_text(face = "bold"),
      axis.title.x = element_text(size = 15, margin = margin(b = 15)),
      axis.title.y = element_text(size = 15, margin = margin(l = 10)),
      axis.text = element_text(size = 12)
    )
}

all.main.density.plot <- density_bitscore_plot(all.main, bit.main$q1, bit.main$median, bit.main$q3)
all.accessory.density.plot <- density_bitscore_plot(all.accessory, bit.accessory$q1, bit.accessory$median, bit.accessory$q3)

# --------------------------
# 8. BAR PLOTS
# --------------------------
bar_plot <- function(data) {
  data <- data %>%
    group_by(species_name, family) %>%
    summarise(count = sum(count), .groups = "drop")
  
  ggplot(data) +
    aes(x = species_name, y = count, fill = family, alpha = count) +
    geom_col(color = "black", linewidth = 0.2) +
    scale_fill_futurama() +
    scale_alpha_continuous(range = c(0.5, 1)) +
    theme_void() +
    labs(fill = "Group", y = "Count") +
    theme(
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11)
    ) +
    guides(alpha = "none")
}

all.full.bar.plot <- bar_plot(all.counted_probe)
all.main.bar.plot <- bar_plot(main.counted_probe)
all.accessory.bar.plot <- bar_plot(accessory.counted_probe)

# --------------------------
# 9. RAINCLOUD PLOTS
# --------------------------
raincloud_bitscore_plot <- function(data) {
  manual_colours <- futurama_unlimited_palette(3, length(unique(data$probe)))
  ggplot(data, aes(y = probe, x = mean_bitscore, fill = probe, color = probe)) +
    ggdist::stat_halfeye(adjust = 0.5, justification = 0, alpha = 0.5, .width = c(0.5, 0.95)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    ggdist::stat_dots(side = "right", dotsize = 0.1, alpha = 0.01, binwidth = 0.2) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    theme_minimal() +
    labs(x = "HSP Bitscore", y = "Probe") +
    theme(
      text = element_text(face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
}

main.full.raincloud.bitscore_probe.plot <- raincloud_bitscore_plot(all.main)
accessory.full.raincloud.bitscore_probe.plot <- raincloud_bitscore_plot(all.accessory)

# --------------------------
# 10. BALLOON PLOTS
# --------------------------
balloon_virus_species_plot <- function(data) {
  manual_colours <- futurama_unlimited_palette(5, length(unique(data$probe)))
  ggplot(data, aes(x = abbreviation, y = factor(species_name, levels = rev(unique(species_name))))) +
    geom_point(aes(color = probe, fill = probe, alpha = count, size = count^3), shape = 15) +
    facet_wrap(~ probe) +
    scale_color_manual(values = manual_colours) +
    scale_fill_manual(values = manual_colours) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_minimal() +
    labs(x = "Virus", y = "Species") +
    theme(
      text = element_text(face = "bold"),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    ) +
    guides(size = "none")
}

main.full.balloon.virus_species.plot <- balloon_virus_species_plot(main.counted_probe)
accessory.full.balloon.virus_species.plot <- balloon_virus_species_plot(accessory.counted_probe)

# --------------------------
# 11. EXPORT PLOTS
# --------------------------
save_plot <- function(name, plot, w = 15, h = 12) {
  ggsave(filename = file.path(output.dir, name), plot = plot, width = w, height = h, dpi = 300)
}

save_plot("main_density.png", all.main.density.plot)
save_plot("accessory_density.png", all.accessory.density.plot)
save_plot("main_raincloud.png", main.full.raincloud.bitscore_probe.plot)
save_plot("accessory_raincloud.png", accessory.full.raincloud.bitscore_probe.plot)
save_plot("full_bar.png", all.full.bar.plot)
save_plot("main_bar.png", all.main.bar.plot)
save_plot("accessory_bar.png", all.accessory.bar.plot)
save_plot("main_balloon.png", main.full.balloon.virus_species.plot, h = 15)
save_plot("accessory_balloon.png", accessory.full.balloon.virus_species.plot, h = 15)
